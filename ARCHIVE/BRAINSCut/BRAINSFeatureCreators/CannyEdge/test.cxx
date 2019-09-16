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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include "itkMaskImageFilter.h"
#include "CannyEdgeCLP.h"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if (inputVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --inputVolume Required! " << std::endl;
  }
  if (outputVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --outputVolume Required! " << std::endl;
  }
  if (violated)
  {
    return EXIT_FAILURE;
  }

  constexpr unsigned int dimension = 3;
  using PixelType = float;
  using InputImage = itk::Image<float, dimension>;
  using OutputImage = itk::Image<unsigned char, dimension>;

  itk::ImageFileReader<InputImage>::Pointer input = itk::ImageFileReader<InputImage>::New();
  input->SetFileName(inputVolume);

  // Set up filter
  itk::CannyEdgeDetectionImageFilter<InputImage, InputImage>::Pointer filter =
    itk::CannyEdgeDetectionImageFilter<InputImage, InputImage>::New();
  itk::SimpleFilterWatcher watcher(filter);
  filter->SetInput(input->GetOutput());
  filter->SetUpperThreshold(30);
  filter->SetLowerThreshold(10);
  filter->SetThreshold(30);
  filter->SetVariance(1.0f);
  filter->SetMaximumError(.01f);

  itk::RescaleIntensityImageFilter<InputImage, OutputImage>::Pointer rescale =
    itk::RescaleIntensityImageFilter<InputImage, OutputImage>::New();
  rescale->SetInput(filter->GetOutput());
  rescale->SetOutputMinimum(0);
  rescale->SetOutputMaximum(255);

  try
  {
    // Generate test image
    itk::ImageFileWriter<OutputImage>::Pointer writer;
    writer = itk::ImageFileWriter<OutputImage>::New();
    writer->UseCompressionOn();
    writer->SetInput(rescale->GetOutput());
    writer->SetFileName(outputVolume);
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    (&err)->Print(std::cerr);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
