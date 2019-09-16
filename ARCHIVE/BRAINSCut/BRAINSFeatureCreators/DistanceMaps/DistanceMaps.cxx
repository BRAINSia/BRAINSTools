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
#include "DistanceMapsCLP.h"
#include <BRAINSCommonLib.h>

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if (inputLabelVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --inputLabelVolume Required! " << std::endl;
  }
  if (inputMaskVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --inputMaskVolume Required! " << std::endl;
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

  using PixelType = float;
  // using PixelType = unsigned long;
  constexpr unsigned int Dimension = 3;

  using LabelImageType = itk::Image<PixelType, Dimension>;
  using MaskImageType = itk::Image<char, Dimension>;
  using ReaderType = itk::ImageFileReader<LabelImageType>;
  using MaskReaderType = itk::ImageFileReader<MaskImageType>;

  ReaderType::Pointer     labelReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  labelReader->SetFileName(inputLabelVolume.c_str());
  maskReader->SetFileName(inputMaskVolume.c_str());

  using MaskFilterType = itk::MaskImageFilter<LabelImageType, MaskImageType, LabelImageType>;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>;
  BinaryThresholdFilterType::Pointer binaryFilter = BinaryThresholdFilterType::New();

  using DistanceMapFilterType = itk::SignedMaurerDistanceMapImageFilter<LabelImageType, LabelImageType>;
  DistanceMapFilterType::Pointer distanceMapFilter = DistanceMapFilterType::New();
  distanceMapFilter->SetInsideIsPositive(true); // Makes all distances positive

  try
  {
    maskFilter->SetInput1(labelReader->GetOutput());
    maskFilter->SetInput2(maskReader->GetOutput());
    maskFilter->SetOutsideValue(0);
    // binaryFilter->SetInput(labelReader->GetOutput());
    binaryFilter->SetInput(maskFilter->GetOutput());
    binaryFilter->SetLowerThreshold(inputTissueLabel);
    binaryFilter->SetUpperThreshold(inputTissueLabel);
    binaryFilter->Update();
    distanceMapFilter->SetInput(binaryFilter->GetOutput());
    distanceMapFilter->UseImageSpacingOff();
    distanceMapFilter->SquaredDistanceOff();
    distanceMapFilter->Update();
    maskReader->Update();
  }
  catch (itk::ExceptionObject & excep)
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    throw;
  }

  using ImageWriterType = itk::ImageFileWriter<DistanceMapFilterType::OutputImageType>;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName(outputVolume);
  imageWriter->SetInput(distanceMapFilter->GetOutput());
  imageWriter->Update();

  return EXIT_SUCCESS;
}
