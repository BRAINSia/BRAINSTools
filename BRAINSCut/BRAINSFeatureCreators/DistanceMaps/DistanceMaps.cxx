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

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if( inputLabelVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputLabelVolume Required! "  << std::endl;
    }
  if( inputMaskVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputMaskVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef float PixelType;
  // typedef unsigned long       PixelType;
  constexpr unsigned int Dimension = 3;

  typedef itk::Image<PixelType,  Dimension>    LabelImageType;
  typedef itk::Image<char,  Dimension>         MaskImageType;
  typedef itk::ImageFileReader<LabelImageType> ReaderType;
  typedef itk::ImageFileReader<MaskImageType>  MaskReaderType;

  ReaderType::Pointer     labelReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  labelReader->SetFileName( inputLabelVolume.c_str() );
  maskReader->SetFileName( inputMaskVolume.c_str() );

  typedef itk::MaskImageFilter<LabelImageType, MaskImageType, LabelImageType> MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> BinaryThresholdFilterType;
  BinaryThresholdFilterType::Pointer binaryFilter = BinaryThresholdFilterType::New();

  typedef itk::SignedMaurerDistanceMapImageFilter<LabelImageType, LabelImageType> DistanceMapFilterType;
  DistanceMapFilterType::Pointer distanceMapFilter = DistanceMapFilterType::New();
  distanceMapFilter->SetInsideIsPositive(true); // Makes all distances positive

  try
    {
    maskFilter->SetInput1( labelReader->GetOutput() );
    maskFilter->SetInput2( maskReader->GetOutput() );
    maskFilter->SetOutsideValue(0);
    // binaryFilter->SetInput(labelReader->GetOutput());
    binaryFilter->SetInput( maskFilter->GetOutput() );
    binaryFilter->SetLowerThreshold(inputTissueLabel);
    binaryFilter->SetUpperThreshold(inputTissueLabel);
    binaryFilter->Update();
    distanceMapFilter->SetInput( binaryFilter->GetOutput() );
    distanceMapFilter->UseImageSpacingOff();
    distanceMapFilter->SquaredDistanceOff();
    distanceMapFilter->Update();
    maskReader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    throw;
    }

  typedef itk::ImageFileWriter<DistanceMapFilterType::OutputImageType> ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName(outputVolume);
  imageWriter->SetInput( distanceMapFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
