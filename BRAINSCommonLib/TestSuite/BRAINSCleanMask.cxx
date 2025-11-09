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
#include "CleanBrainLabelMap.h"
#include "itkIO.h"

int
main(int argc, char ** argv)
{
  if (argc < 3)
  {
    std::cerr << "Usage: BRAINSCleanMask inputLabelMap outputLabelMap" << std::endl;
    return 1;
  }
  using ImageType = itk::Image<unsigned char, 3>;

  const std::string inputName(argv[1]);

  const std::string outputName(argv[2]);

  const ImageType::Pointer input = itkUtil::ReadImage<ImageType>(inputName);
  const ImageType::Pointer output = CleanBrainLabelMap<ImageType, ImageType>(input);
  itkUtil::WriteImage<ImageType>(output, outputName);
}
