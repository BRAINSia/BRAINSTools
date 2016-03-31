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
 * Author: Ali Ghayoor
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

 /*
  This program gets several fcsv file each one contains several landmarks with the same name but slightly different coordinates.
  For EACH landmark we compute the average coordination.

  The final output of this program is a new landmark fcsv file which contains the average coordinate of each landmark.

  Usage:
  $/GenerateAverageLmkFile \
    --inputLandmarkFiles lmk1.fcsv,lmk2.fcsv,...,lmkn.fcsv \
    --outputLandmarkFile outputAveLmk.fcsv
 */

#include "itkImage.h"
#include <cmath>
#include "Slicer3LandmarkIO.h"

#include "GenerateAverageLmkFileCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::vector<std::string> inputFileNames;
  if( inputLandmarkFiles.size() > 0 )
    {
    inputFileNames = inputLandmarkFiles;
    }
  else
    {
    std::cerr << "ERROR: No input file name is defined" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int k = inputFileNames.size(); // The number of input landmark files
  if( k==0 )
    {
    std::cerr << "ERROR: Number of input landmark files is zero!" << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Computing the average file for " << k << " input landmark files..." << std::endl;

  unsigned int numNamedLandmarks = 0;

  std::map<std::string, LandmarkPointType> LandmarksAverageMap;

  typedef std::vector<std::map<std::string, LandmarkPointType> > LandmarksMapTypeVec;
  LandmarksMapTypeVec LandmarksMapVector;

  // LandmarksMapType is as "std::map<std::string, LandmarkPointType>" which means a map between landmarks and their
  // coordinates.
  // For each input landmark file this LandmarksMapType is computed and is set in a vector: "LandmarksMapTypeVec"

  LandmarksMapType temp;
  for( unsigned int i = 0; i < k; i++ )
    {
    temp = ReadSlicer3toITKLmk( inputFileNames[i] );
    LandmarksMapVector.push_back(temp);
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
  // Computing the average coordinate for each landmark
  for( unsigned int j = 0; j < numNamedLandmarks; j++ )
    {
    double x_ave = 0.0;
    double y_ave = 0.0;
    double z_ave = 0.0;
    std::string name = LandmarksNames[j];
    for( unsigned int i = 0; i < k; i++ )
      {
      x_ave += LandmarksMapVector[i][name][0];
      y_ave += LandmarksMapVector[i][name][1];
      z_ave += LandmarksMapVector[i][name][2];
      }
    x_ave = x_ave / static_cast<double>(k);
    y_ave = y_ave / static_cast<double>(k);
    z_ave = z_ave / static_cast<double>(k);

    LandmarksAverageMap[name][0] = x_ave;
    LandmarksAverageMap[name][1] = y_ave;
    LandmarksAverageMap[name][2] = z_ave;
    }

  WriteITKtoSlicer3Lmk( outputLandmarkFile, LandmarksAverageMap );

  return EXIT_SUCCESS;
}
