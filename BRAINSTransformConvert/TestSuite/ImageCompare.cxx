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
#include "itkIO.h"
#include "itkImageRegionIterator.h"
#include <iostream>
int
main(int argc, char * *argv)
{
  if( argc < 3 )
    {
    std::cerr << "Usage: ImageCompare imagea imageb" << std::endl;
    exit(1);
    }

  std::string input1Name(argv[1]), input2Name(argv[2]);

  using ImageType = itk::Image<short, 3>;

  ImageType::Pointer image1, image2;

  try
    {
    image1 = itkUtil::ReadImage<ImageType>(input1Name);
    }
  catch( ... )
    {
    std::cerr << "Error reading " << input1Name << std::endl;
    exit(1);
    }
  try
    {
    image2 = itkUtil::ReadImage<ImageType>(input2Name);
    }
  catch( ... )
    {
    std::cerr << "Error reading " << input2Name << std::endl;
    exit(1);
    }
  if( image1.IsNull() )
    {
    std::cerr << "Error reading " << input1Name << std::endl;
    exit(1);
    }
  if( image2.IsNull() )
    {
    std::cerr << "Error reading " << input2Name << std::endl;
    exit(1);
    }

  using ImageIteratorType = itk::ImageRegionIterator<ImageType>;
  ImageIteratorType it1(image1, image1->GetLargestPossibleRegion() );
  ImageIteratorType it2(image2, image2->GetLargestPossibleRegion() );

  while( !it1.IsAtEnd() && !it2.IsAtEnd() )
    {
    short diff = it1.Value() - it2.Value();
    if( diff < 0 )
      {
      diff = -diff;
      }
    if( diff > 1 )
      {
      std::cerr << "Mismatch between " << input1Name
                << " and " << input2Name << std::endl;
      break;
      }
    ++it1; ++it2;
    }

  if( !it1.IsAtEnd() && !it2.IsAtEnd() )
    {
    exit(1);
    }
  exit(0);
}
