// Author: Hans J. Johnson

#include "LandmarksCompareCLP.h"
#include "Slicer3LandmarkIO.h"
#include <stdlib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

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
      for( unsigned int i = 0 ; i < 3 ; ++i )
        {
        if ( vcl_abs(lmk1iter->second[i] - lmk2iter->second[i]) > tolerance )
          {
          std::cout << "lmks differ by greater than tolerance" << std::endl;
          std::cout << "| "<< lmk1iter->second[i]
                    << " - " << lmk2iter->second[i]
                    << " | > " << tolerance << std::endl;
          allSame = false;
          }
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
