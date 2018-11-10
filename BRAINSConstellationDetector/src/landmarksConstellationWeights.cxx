/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*
*/

// ==============================================================================================

// 1- BCD is run through all the training data set.
//    We use an inputTrainingList the same as what we use for the BRAINSConstellationModeler.
//    For each training volume, just the "outputLandmarksInInputSpace" is needed.

// 2- The results of the BCD is compared with the training landmarks files. For each landmark point, the average and the
// STD (standard deviation) of error is calculated.

// 3- Resultant STD values are used to compute the wights such that the weight of a landmark point, which has the lowest
// std, is mapped to ONE, and the weights of the other landmarks are calculated correspondingly. The landmark with the
// highest STD has the lowest weight that is close to zero.

// 4- Resultant weight map is written to the disk as a CSV file with the extension of "wts".

// ==============================================================================================

// I N C L U D E S ////////////////////////////////////////////////////////////

#include "BRAINSConstellationDetectorPrimary.h"
#include "landmarksConstellationCommon.h"
#include "landmarksConstellationTrainingDefinitionIO.h"
#include "landmarkIO.h"
#include "Slicer3LandmarkWeightIO.h"
#include "landmarksConstellationWeightsCLP.h"

// D E F I N E S //////////////////////////////////////////////////////////////

using LandmarksDistanceMapType = std::map<std::string, std::vector<double> >;
using LandmarksValueMapType = std::map<std::string, double>;
using LandmarksMapTypeVec = std::vector<std::map<std::string, LandmarkPointType> >;

// F U N C T I O N S //////////////////////////////////////////////////////////

bool value_comparer(const std::pair<std::string, float> & i,
                    const std::pair<std::string, float> & j)
{
  return i.second < j.second;
}

double get_max_value(LandmarksValueMapType const & m)  // returns the max value of a map
{
  return std::max_element(m.begin(), m.end(), value_comparer)->second;
}

double get_min_value(LandmarksValueMapType const & m)  // returns the min value of a map
{
  return std::min_element(m.begin(), m.end(), value_comparer)->second;
}

// M A I N ////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  unsigned int             numNamedLandmarks = 0;
  double                   d0, d1, d2, dist;
  LandmarksDistanceMapType LandmarksDistanceMap;
  LandmarksValueMapType    LandmarksAverageMap; // for average
  LandmarksValueMapType    LandmarksSTDMap;     // for standard deviation
  LandmarksValueMapType    LandmarksWeightMap;  // for weights
  LandmarksMapTypeVec      LandmarksMapVector;

  if( ( inputTrainingList.compare("") == 0 )
      || ( outputWeightsList.compare("") == 0 )
      || ( inputTemplateModel.compare("") == 0 )
      || ( llsModel.compare("") == 0 ) )
    {
    std::cerr << "To run the program please specify the training list filename and the "
              << "output weight list filename."
              << "\nAlso, the \"LLSModel\" and the \"templateModel\" of the target landmarks list should be defined."
              << std::endl;
    std::cerr << "Type " << argv[0] << " -h for more help." << std::endl;
    return EXIT_FAILURE;
    }

  landmarksConstellationTrainingDefinitionIO mDef(inputTrainingList);
  LandmarksMapVector.clear();
  // For each training case, it loads the landmarks file into the LandmarksMapVector
  for( unsigned int currentDataset = 0; currentDataset < mDef.GetNumDataSets(); currentDataset++ )
    {
    std::cout << "PROCESSING:" <<  mDef[currentDataset].GetImageFilename() << std::endl;

    SImageType::Pointer volOrig = itkUtil::ReadImage<SImageType>( mDef[currentDataset].GetImageFilename() );
    if( volOrig.IsNull() )
      {
      printf( "\nCould not open image %s, aborting ...\n\n", mDef[currentDataset].GetImageFilename().c_str() );
      return EXIT_FAILURE;
      }

    LandmarksMapVector.push_back( mDef[currentDataset] );
    }

  // running the BCD on each training case, and loading the output landmarks into the LandmarkMapVector
  BRAINSConstellationDetectorPrimary BCD;
  BCD.SetInputTemplateModel( inputTemplateModel );
  BCD.SetLLSModel( llsModel );
  for( unsigned int currentDataset = 0; currentDataset < mDef.GetNumDataSets(); currentDataset++ )
    {
    std::cout << "\n====================================================================================" << std::endl;
    std::cout << "RUNNING BRAINSConstellationDetector ON: " <<  mDef[currentDataset].GetImageFilename() << std::endl;
    BCD.SetInputVolume( mDef[currentDataset].GetImageFilename() );
    BCD.Compute();
    LandmarksMapVector.push_back( BCD.GetOutputLandmarksInInputSpace() );
    }

  // Putting all the landmark names in a vector of strings and computing the number of landmarks
  std::vector<std::string> LandmarksNames;
  for( LandmarksMapType::const_iterator it = LandmarksMapVector[0].begin(); it != LandmarksMapVector[0].end(); ++it )
    {
    if( ( it->first ).compare("") != 0 )
      {
      LandmarksNames.push_back(it->first);
      numNamedLandmarks++;
      }
    }

  unsigned int k = mDef.GetNumDataSets();
  // So the LandmarkMapVector contains 2*k landmark files

  std::cout << "\n====================================================================================" << std::endl;
  std::cout << "Number of training cases = " << k << std::endl;
  std::cout << "Number of landmarks = " << numNamedLandmarks << std::endl;
  /*
  // TEST PRINT FOR TRACKING
  for (unsigned int i=0; i<(2*k); i++)
  {
      std::cout << "\n\nValues of landmarks(" << i << "): " << std::endl;
      for( LandmarksMapType::const_iterator it = LandmarksMapVector[i].begin(); it != LandmarksMapVector[i].end(); ++it )
      {
          std::cout << it->first << " : (" << it->second[0] << ", " << it->second[1] << ", " << it->second[2] << ")" << std::endl;
      }
      std::cout << "==============================================================" << std::endl;
  }
  */
  // Computing the average coordinate for each landmark
  for( unsigned int j = 0; j < numNamedLandmarks; j++ )
    {
    std::string name = LandmarksNames[j];
    for( unsigned int i = 0; i < k; i++ )
      {
      d0 = pow( LandmarksMapVector[i][name][0] - LandmarksMapVector[i + k][name][0], 2);
      d1 = pow( LandmarksMapVector[i][name][1] - LandmarksMapVector[i + k][name][1], 2);
      d2 = pow( LandmarksMapVector[i][name][2] - LandmarksMapVector[i + k][name][2], 2);

      dist = sqrt(d0 + d1 + d2);

      LandmarksDistanceMap[name].push_back(dist);
      }
    }
  // Computing the Average of the distances for each landmark
  for( unsigned int j = 0; j < numNamedLandmarks; j++ )
    {
    std::string name = LandmarksNames[j];
    double      sum = 0;
    for( unsigned int i = 0; i < k; i++ )
      {
      sum += LandmarksDistanceMap[name][i];
      }
    LandmarksAverageMap[name] = sum / k;
    }
  // Computing the Variance and the Standard Deviation of the distances for each landmark
  for( unsigned int j = 0; j < numNamedLandmarks; j++ )
    {
    std::string name = LandmarksNames[j];
    double      sum = 0;
    for( unsigned int i = 0; i < k; i++ )
      {
      sum += pow( LandmarksDistanceMap[name][i] - LandmarksAverageMap[name], 2);
      }
    LandmarksSTDMap[name] = sqrt(sum / k);
    }

  /*
  // TEST PRINT FOR TRACKING
  std::cout << "\n=============Standard Deviation of Error for Each Landmarks====================" << std::endl;
for (unsigned int j=0; j<numNamedLandmarks; j++)
{
  std::string name = LandmarksNames[j];
  std::cout << "-" << name << ": " << LandmarksSTDMap[name] << std::endl;
}
std::cout << "==============================================================\n" << std::endl;
  */

  // Computing LandmarksWeightMap
  // the weight of a landmark point, which has the lowest std, should be mapped to ONE, and weights of the other
  // landmarks should be calculated correspondingly
  double Margin = 0.05;
  double minValue = get_min_value( LandmarksSTDMap );
  double maxValue = get_max_value( LandmarksSTDMap );
  double UpperBound = Margin + minValue + maxValue;
  /*
  // TEST PRINT FOR TRACKING
  std::cout << "minValue = " << minValue << std::endl;
  std::cout << "maxValue = " << maxValue << std::endl;
  */
  for( unsigned int j = 0; j < numNamedLandmarks; j++ )
    {
    std::string name = LandmarksNames[j];
    LandmarksWeightMap[name] = ( UpperBound - LandmarksSTDMap[name] ) / ( Margin + maxValue);
    }

  // Writing the weights into the file
  WriteITKtoSlicer3LmkWts( outputWeightsList, LandmarksWeightMap );
  std::cout << "The output landmark weights list file is written." << std::endl;

  return EXIT_SUCCESS;
}
