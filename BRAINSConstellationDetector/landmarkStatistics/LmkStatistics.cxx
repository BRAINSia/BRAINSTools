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

// This program gets a row of the fcsv files. The first argument is the number of pair landmark files for example
// (e0326_10335.fcsv , k0326_10335.fcsv)
// After that the fcsv file names come consequently.
// For instance:
//              LMKstatis 10 E_lmk1.fcsv E_lmk2.fcsv ... E_lmk10.fcsv K_lmk1.fcsv K_lmk2.fcsv ... K_lmk10.fcsv

// *** FOR EACH PAIR FILES: ****
// We have two fcsv file for each volume case (e.g. we have 20 cases, and there are two landmark files (.fcsv)
// correspondence to each case).
// Each fcsv file contains several landmarks (e.g. 30 landmarks as AC, PC, ...)
// Therefore, for each landmark point in each case, we have two correspondent coordinate in the two correspondent fcsv
// files.
// For example for the PC point:

/*
         lmk1.fcsv   lmk2.fcsv
case 1:     PC11         PC12
case 2:     PC21         PC22
 .           .            .
 .           .            .
 .           .            .
case n:     PCn1         PCn2

*/

// First, for EACH landmark we compute the distance between its coordinates in each case,
// Then, the average, variance and standard deviation of these distances are calculated over all n cases.

#include "Slicer3LandmarkIO.h"
#include "itkImage.h"
#include "math.h"
#include <cmath>


int
main(int argc, char * argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage : " << argv[0] << " <#of landmark pairs> <lmk-pair1> <lmk-pair2> ... <lmk-pairN>" << std::endl;
    return EXIT_FAILURE;
  }
  const unsigned int k = std::stoi(argv[1]); // Number of landmark pairs
  // So the number of input landmark files = 2*k

  if ((argc % 2 != 0) || (static_cast<unsigned int>(argc) != (2 * (k + 1))))
  {
    std::cout << " First argument indicates the number of landmark pairs for comparison.\n"
              << "Two landmarks files are needed corresponding to each data case.\n"
              << "Therefore, you must insert an even number of landmark files at the command line.\n"
              << std::endl;

    return EXIT_FAILURE;
  }

  unsigned int numNamedLandmarks = 0;
  double       d0 = NAN;

  double       d1 = NAN;

  double       d2 = NAN;

  double       dist = NAN;
  using LandmarksDistanceMapType = std::map<std::string, std::vector<double>>;
  LandmarksDistanceMapType      LandmarksDistanceMap;
  std::map<std::string, double> LandmarksAverageMap;  // for average
  std::map<std::string, double> LandmarksVarianceMap; // for variance
  std::map<std::string, double> LandmarksSTDMap;      // for standard deviation

  using LandmarksMapTypeVec = std::vector<std::map<std::string, LandmarkPointType>>;
  LandmarksMapTypeVec LandmarksMapVector;

  // LandmarksMapType is as "std::map<std::string, LandmarkPointType>" which means a map between landmarks and their
  // coordinates.
  // For each input landmark file this LandmarksMapType is computed and is set in a vector: "LandmarksMapTypeVec"

  LandmarksMapType temp;
  for (unsigned int i = 0; i < (2 * k); i++)
  {
    temp = ReadSlicer3toITKLmk(argv[i + 2]);
    LandmarksMapVector.push_back(temp);
  }

  // Putting all the landmark names in a vector of strings and computing the number of landmarks
  std::vector<std::string> LandmarksNames;
  for (LandmarksMapType::const_iterator it = LandmarksMapVector[0].begin(); it != LandmarksMapVector[0].end(); ++it)
  {
    if ((it->first).compare("") != 0)
    {
      LandmarksNames.push_back(it->first);
      numNamedLandmarks++;
    }
  }

  std::cout << "\n\nNumber of landmarks = " << numNamedLandmarks << std::endl;

  // Computing the average coordinate for each landmark
  std::cout << "\n=============Distance Values For Each Landmarks====================" << std::endl;
  for (unsigned int j = 0; j < numNamedLandmarks; j++)
  {
    std::string name = LandmarksNames[j];
    for (unsigned int i = 0; i < k; i++)
    {
      d0 = pow(LandmarksMapVector[i][name][0] - LandmarksMapVector[i + k][name][0], 2);
      d1 = pow(LandmarksMapVector[i][name][1] - LandmarksMapVector[i + k][name][1], 2);
      d2 = pow(LandmarksMapVector[i][name][2] - LandmarksMapVector[i + k][name][2], 2);

      dist = sqrt(d0 + d1 + d2);
      // dist = sqrt(d0); // distance in left/right direction
      // dist = sqrt(d1); // distance in anterior/posterior direction
      // dist = sqrt(d2); // distance in superior/inferior direction

      LandmarksDistanceMap[name].push_back(dist);
    }

    std::cout << "-" << name << std::endl;
    for (unsigned int l = 0; l < k; l++)
    {
      std::cout << "(" << l + 1 << "):" << LandmarksDistanceMap[name][l] << ", ";
    }
    std::cout << "\n==============================================================" << std::endl;
  }

  // Computing the Average of the distances for each landmark
  std::cout << "\n=============Average Error for Each Landmarks====================" << std::endl;
  for (unsigned int j = 0; j < numNamedLandmarks; j++)
  {
    std::string name = LandmarksNames[j];
    double      sum = 0;
    for (unsigned int i = 0; i < k; i++)
    {
      sum += LandmarksDistanceMap[name][i];
    }
    LandmarksAverageMap[name] = sum / k;
    std::cout << "-" << name << ": " << LandmarksAverageMap[name] << std::endl;
  }

  // Computing the Variance of the distances for each landmark
  std::cout << "\n=============Variance of Error for Each Landmarks====================" << std::endl;
  for (unsigned int j = 0; j < numNamedLandmarks; j++)
  {
    std::string name = LandmarksNames[j];
    double      sum = 0;
    for (unsigned int i = 0; i < k; i++)
    {
      sum += pow(LandmarksDistanceMap[name][i] - LandmarksAverageMap[name], 2);
    }
    LandmarksVarianceMap[name] = sum / k;
    LandmarksSTDMap[name] = sqrt(sum / k);
    std::cout << "-" << name << ": " << LandmarksVarianceMap[name] << std::endl;
  }

  // Computing the Standard Deviation of the distances for each landmark
  std::cout << "\n=============Standard Deviation of Error for Each Landmarks====================" << std::endl;
  for (unsigned int j = 0; j < numNamedLandmarks; j++)
  {
    std::string name = LandmarksNames[j];
    std::cout << "-" << name << ": " << LandmarksSTDMap[name] << std::endl;
  }
  std::cout << "==============================================================\n" << std::endl;

  return EXIT_SUCCESS;
}
