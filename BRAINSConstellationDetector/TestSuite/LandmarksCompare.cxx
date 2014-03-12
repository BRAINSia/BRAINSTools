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
// Author: Hans J. Johnson

#include "LandmarksCompareCLP.h"
#include "Slicer3LandmarkIO.h"
#include <stdlib.h>
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // load corresponding landmarks in EMSP aligned space from file if possible
  const LandmarksMapType landmarks1 = ReadSlicer3toITKLmk( inputLandmarkFile1 );
  const LandmarksMapType landmarks2 = ReadSlicer3toITKLmk( inputLandmarkFile2 );

  if( landmarks1.empty() )
    {
    std::cout << "ERROR: " << inputLandmarkFile1 << " is empty" << std::endl;
    return EXIT_FAILURE;
    }
  if( landmarks2.empty() )
    {
    std::cout << "ERROR: " << inputLandmarkFile2 << " is empty" << std::endl;
    return EXIT_FAILURE;
    }
  if( landmarks1.size() != landmarks2.size() )
    {
    std::cout << "ERROR: number of landmarks differ." << std::endl;
    return EXIT_FAILURE;
    }
  LandmarksMapType::const_iterator lmk1iter = landmarks1.begin();
  bool allSame = true;
  while( lmk1iter != landmarks1.end() )
    {
    const LandmarksMapType::const_iterator lmk2iter = landmarks2.find(lmk1iter->first);
    if ( lmk2iter == landmarks2.end() )
      {
      std::cout << "Missing landmark in second file" << lmk1iter->first << std::endl;
      allSame = false;
      continue;
      }
    else
      {
      bool thisLmkOK = true;
      for( unsigned int i = 0 ; i < 3 ; ++i )
        {
        const double error_term = vcl_abs(lmk1iter->second[i] - lmk2iter->second[i]);
        if ( error_term > tolerance )
          {
          std::cout << "\nFAIL: lmk" << lmk1iter->first << "[" << i << "] differ by greater than tolerance" << std::endl;
          std::cout << "FAIL: | "<< lmk1iter->second[i] << " - " << lmk2iter->second[i] << " | = " << error_term << " is greater than " << tolerance << std::endl;
          allSame = false;
          thisLmkOK = false;
          }
        }
      if (thisLmkOK)
        {
        std::cout << "PASS:  lmk" << lmk1iter->first << std::endl;
        }
      }
    ++lmk1iter;
    }
  if( allSame )
    {
    std::cout << "The landmark files are identical!" << std::endl;
    return EXIT_SUCCESS;
    }
  else
    {
    std::cout << "ERROR: The landmark files are too different!" << std::endl;
    return EXIT_FAILURE;
    }
}
