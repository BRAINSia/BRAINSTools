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
#include <iostream>
#include <fstream>
#include <sstream>
#include <itkThresholdImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>

#include <BRAINSCommonLib.h>
#include "GenerateCsfClippedFromClassifiedImageCLP.h"

int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  using PixelType = float;
  constexpr unsigned int Dim = 3;
  using ImageType = itk::Image<PixelType, Dim>;

  using ImageReaderType = itk::ImageFileReader<ImageType>;

  ImageReaderType::Pointer classified_imageReader = ImageReaderType::New();

  classified_imageReader->SetFileName(inputCassifiedVolume);

  using ThresholdFilterType = itk::ThresholdImageFilter<ImageType>;

  ThresholdFilterType::Pointer classified_gradientFilter = ThresholdFilterType::New();

  classified_gradientFilter->SetInput(classified_imageReader->GetOutput());

  classified_gradientFilter->SetLower(50);
  classified_gradientFilter->Update();

  // writer setting
  std::cout << "Writing output ... " << std::endl;
  using WriterType = itk::ImageFileWriter<ImageType>;

  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageType, ImageType>;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);

  rescaler->SetInput(classified_gradientFilter->GetOutput());
  writer->SetFileName(outputVolume);
  writer->SetInput(rescaler->GetOutput());
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & exp)
  {
    std::cerr << "ExceptionObject with writer" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }
  return 0;
}
