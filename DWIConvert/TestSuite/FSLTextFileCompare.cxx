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
#include "DWIConvertUtils.h"

int main(int argc, char *argv[])
{
  if( argc < 5 )
    {
    std::cerr << argv[0]
              << " Usage: FSLTestCompare bvecfile1 bvecfile2 bvalfile1 bvalfile2"
              << std::endl;
    return EXIT_FAILURE;
    }
  std::string bvecfile1(argv[1]);
  std::string bvecfile2(argv[2]);
  std::string bvalfile1(argv[3]);
  std::string bvalfile2(argv[4]);

  unsigned int                      bVecCount1;
  DWIMetaDataDictionaryValidator::GradientTableType bvecs1;
  if( ReadBVecs(bvecs1, bVecCount1, bvecfile1, false) != EXIT_SUCCESS )
    {
    std::cerr << "Can't read " << bvecfile1 << std::endl;
    return EXIT_FAILURE;
    }
  unsigned int                      bVecCount2;
  DWIMetaDataDictionaryValidator::GradientTableType bvecs2;
  if( ReadBVecs(bvecs2, bVecCount2, bvecfile2, false) != EXIT_SUCCESS )
    {
    std::cerr << "Can't read " << bvecfile2 << std::endl;
    return EXIT_FAILURE;
    }
  if( bVecCount1 != bVecCount2 )
    {
    std::cerr << "Size mismatch: " << bvecfile1 << " (" << bVecCount1
              << ") " << bvecfile2 << " (" << bVecCount2 << ")"
              << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int        bValCount1;
  std::vector<double> bvals1;
  if( ReadBVals(bvals1, bValCount1, bvalfile1) != EXIT_SUCCESS )
    {
    std::cerr << "Can't read " << bvalfile1 << std::endl;
    return EXIT_FAILURE;
    }
  unsigned int        bValCount2;
  std::vector<double> bvals2;
  if( ReadBVals(bvals2, bValCount2, bvalfile2) != EXIT_SUCCESS )
    {
    std::cerr << "Can't read " << bvalfile2 << std::endl;
    return EXIT_FAILURE;
    }

  if( bValCount1 != bValCount2 )
    {
    std::cerr << "Size mismatch: " << bvalfile1 << " (" << bValCount1
              << ") " << bvalfile2 << " (" << bValCount2 << ")"
              << std::endl;
    return EXIT_FAILURE;
    }

  if( !CloseEnough(bvecs1, bvecs2, 1000.0) )
    {
    PrintVec(bvecs1);
    PrintVec(bvecs2);
    return EXIT_FAILURE;
    }

  if( !CloseEnough(bvals1, bvals2, 1000.0) )
    {
    PrintVec(bvals1);
    PrintVec(bvals2);
    return EXIT_FAILURE;
    }
}
