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
// Author: Hans J. Johnson, Ali Ghayoor

#include "LandmarksCompareCLP.h"
#include "Slicer3LandmarkIO.h"
#include <cstdlib>
#include <BRAINSCommonLib.h>

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // inputLandmarkFile2 can contain several baselines, and this program is modified
  // to work with different baselines
  std::vector<std::string> inputLandmarkFile2Names;
  inputLandmarkFile2Names = inputLandmarkFile2;

  const unsigned int numBaselines = inputLandmarkFile2Names.size(); // The number of landmark baseline files
  if (numBaselines == 0)
  {
    std::cerr << "ERROR: No input baseline file is defined!" << std::endl;
    return EXIT_FAILURE;
  }

  bool testIterationIsPassed = false;

  for (unsigned int l = 0; l < numBaselines; ++l)
  {
    std::cout << "\nCompare the input landmark file with the baseline files number: " << l + 1 << std::endl;

    // load corresponding landmarks in EMSP aligned space from file if possible
    const LandmarksMapType       landmarks1 = ReadSlicer3toITKLmk(inputLandmarkFile1);
    const LandmarksMapType       landmarks2 = ReadSlicer3toITKLmk(inputLandmarkFile2Names[l]);
    const LandmarksWeightMapType weight_values = ReadLandmarkWeights(weights);

    if (landmarks1.empty())
    {
      std::cout << "ERROR: " << inputLandmarkFile1 << " is empty" << std::endl;
      return EXIT_FAILURE;
    }
    if (landmarks2.empty())
    {
      std::cout << "ERROR: " << inputLandmarkFile2Names[l] << " is empty" << std::endl;
      return EXIT_FAILURE;
    }
    if (landmarks1.size() != landmarks2.size())
    {
      std::cout << "ERROR: number of landmarks differ." << std::endl;
      return EXIT_FAILURE;
    }
    if (weight_values.empty())
    {
      std::cout << "ERROR: Missing weights" << std::endl;
      return EXIT_FAILURE;
    }
    auto lmk1iter = landmarks1.begin();
    bool                             allSame = true;
    while (lmk1iter != landmarks1.end())
    {
      const LandmarksMapType::const_iterator       lmk2iter = landmarks2.find(lmk1iter->first);
      const LandmarksWeightMapType::const_iterator wtsiter = weight_values.find(lmk1iter->first);
      const double lmk_tolerance = (wtsiter == weight_values.end()) ? 5.0 : wtsiter->second;
      if (wtsiter == weight_values.end())
      {
        std::cout << "Missing weight in file for landmark: " << lmk1iter->first << std::endl;
      }

      if (lmk2iter == landmarks2.end())
      {
        std::cout << "Missing landmark in second file: " << lmk1iter->first << std::endl;
        allSame = false;
        continue;
      }
      else
      {
        bool         thisLmkOK = true;
        const double error_term = lmk1iter->second.EuclideanDistanceTo(lmk2iter->second);

        if (error_term > lmk_tolerance && error_term > tolerance)
        {
          std::cout << "FAILED: lmk " << lmk1iter->first << "\t\t differ    ||" << lmk1iter->second << " - "
                    << lmk2iter->second << "|| = " << error_term << "\t > " << lmk_tolerance << " wts " << lmk_tolerance
                    << " tol: " << tolerance << std::endl;
          allSame = false;
          thisLmkOK = false;
        }
        if (thisLmkOK)
        {
          std::cout << "  PASSED: lmk " << lmk1iter->first << "\t\t are same ||" << lmk1iter->second << " - "
                    << lmk2iter->second << "|| = " << error_term << "\t > " << lmk_tolerance << " wts " << lmk_tolerance
                    << " tol: " << tolerance << std::endl;
        }
      }
      ++lmk1iter;
    }
    if (allSame)
    {
      std::cout << "The input landmark file is identical to the baseline file: " << l + 1 << "!" << std::endl;
      testIterationIsPassed = true;
      break;
    }
    else
    {
      std::cout << "WARINING: The input landmark file is too different than the baseline file: " << l + 1 << "!"
                << std::endl;
    }
  }

  if (testIterationIsPassed)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    std::cout << "ERROR: The landmark files are too different!" << std::endl;
    return EXIT_FAILURE;
  }
}
