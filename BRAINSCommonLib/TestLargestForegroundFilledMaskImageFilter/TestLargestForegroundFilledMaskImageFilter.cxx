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
#include "itkLargestForegroundFilledMaskImageFilter.h"
#include "itkIO.h"
#include <string>

int
main(int argc, char * *argv)
{
  if( argc < 3 )
    {
    std::cerr << "TestLargestForegrounFilledMaskImageFilter <input image> <output image>"
              << std::endl;
    return EXIT_FAILURE;
    }
  using ImageType = itk::Image<float, 3>;
  using FilterType = itk::LargestForegroundFilledMaskImageFilter<ImageType>;
  std::string         inputname(argv[1]);
  std::string         outputname(argv[2]);
  ImageType::Pointer  image = itkUtil::ReadImage<ImageType>(inputname);
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->Update();
  ImageType::Pointer outputImage = filter->GetOutput();
  itkUtil::WriteImage<ImageType>(outputImage, outputname);
  return EXIT_SUCCESS;
}
