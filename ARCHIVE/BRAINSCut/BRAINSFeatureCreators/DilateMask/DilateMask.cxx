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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkBinaryThresholdImageFilter.h"

#include "DilateMaskCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  constexpr unsigned int Dimension = 3;

  using InputPixelType = float;
  using OutputPixelType = unsigned char;

  using InputImageType = itk::Image<InputPixelType,  Dimension>;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  using ReaderType = itk::ImageFileReader<InputImageType>;
  using WriterType = itk::ImageFileWriter<OutputImageType>;

  using ThresholdFilterType = itk::BinaryThresholdImageFilter<InputImageType, OutputImageType>;

  using StructuringElementType = itk::BinaryBallStructuringElement<
      InputPixelType,
      Dimension>;

  using DilateFilterType = itk::BinaryDilateImageFilter<
      OutputImageType,
      OutputImageType,
      StructuringElementType>;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writingFilter = WriterType::New();

  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();

  DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
  StructuringElementType    structuringElement;
  structuringElement.SetRadius( sizeStructuralElement);  // 3x3 structuring element
  structuringElement.CreateStructuringElement();
  binaryDilate->SetKernel( structuringElement );

  reader->SetFileName( inputBinaryVolume );

  writingFilter->SetFileName( outputVolume );

  thresholder->SetInput( reader->GetOutput() );

  InputPixelType background =   0;
  InputPixelType foreground = 1;

  thresholder->SetOutsideValue( background );
  thresholder->SetInsideValue(  foreground );

  thresholder->SetLowerThreshold( lowerThreshold );
  // thresholder->SetUpperThreshold( upperThreshold );

  binaryDilate->SetInput( thresholder->GetOutput() );

  binaryDilate->SetDilateValue( foreground );

  writingFilter->SetInput( binaryDilate->GetOutput() );
  writingFilter->Update();

  return EXIT_SUCCESS;
}
