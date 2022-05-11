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
#include "BRAINSCleanMaskCLP.h"
#include "BRAINSThreadControl.h"

int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if (inputVolume == "")
  {
    std::cerr << "No input volume name given" << std::endl;
    return EXIT_FAILURE;
  }
  if (outputVolume == "")
  {
    std::cerr << "No output volume name given" << std::endl;
    return EXIT_FAILURE;
  }

  using ImageType = itk::Image<unsigned char, 3>;

  ImageType::Pointer input;
  try
  {
    input = itkUtil::ReadImage<ImageType>(inputVolume);
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "error reading " << inputVolume << std::endl << e << std::endl;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cerr << "Unable to open " << inputVolume << std::endl;
    return EXIT_FAILURE;
  }
  ImageType::Pointer output;
  try
  {
    output = CleanBrainLabelMap<ImageType, ImageType>(input);
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cerr << "Error during processing of " << inputVolume << std::endl;
    return EXIT_FAILURE;
  }

  try
  {
    itkUtil::WriteImage<ImageType>(output, outputVolume);
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "error writing " << inputVolume << std::endl << e << std::endl;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cerr << "Unable to write " << outputVolume << std::endl;
    return EXIT_FAILURE;
  }
  return 0;
}
