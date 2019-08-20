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
#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMeanImageFilter.h"
#include "ErodeImageCLP.h"
#include <BRAINSCommonLib.h>

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if ( inputVolume.size() == 0 )
  {
    violated = true;
    std::cout << "  --inputVolume Required! " << std::endl;
  }
  if ( inputMaskVolume.size() == 0 )
  {
    violated = true;
    std::cout << "  --inputMaskVolume Required! " << std::endl;
  }
  if ( inputMaskVolume.size() == 0 )
  {
    violated = true;
    std::cout << "  --inputMaskVolume Required! " << std::endl;
  }
  if ( violated )
  {
    return EXIT_FAILURE;
  }

  using PixelType = float;
  // using PixelType = unsigned long;
  constexpr unsigned int Dimension = 3;

  using ImageType = itk::Image< PixelType, Dimension >;
  using MaskImageType = itk::Image< char, Dimension >;
  using ReaderType = itk::ImageFileReader< ImageType >;
  using MaskReaderType = itk::ImageFileReader< MaskImageType >;

  ReaderType::Pointer     imageReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  imageReader->SetFileName( inputVolume.c_str() );
  maskReader->SetFileName( inputMaskVolume.c_str() );

  using MaskFilterType = itk::MaskImageFilter< ImageType, MaskImageType, ImageType >;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  using StructuringElementType = itk::BinaryBallStructuringElement< PixelType, Dimension >;

  using ErodeFilterType = itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementType >;
  ErodeFilterType::Pointer grayscaleErodeFilter = ErodeFilterType::New();

  StructuringElementType structuringElement;
  structuringElement.SetRadius( inputRadius );
  structuringElement.CreateStructuringElement();
  grayscaleErodeFilter->SetKernel( structuringElement );

  try
  {
    maskFilter->SetInput1( imageReader->GetOutput() );
    maskFilter->SetInput2( maskReader->GetOutput() );
    maskFilter->SetOutsideValue( 0 );

    grayscaleErodeFilter->SetInput( maskFilter->GetOutput() );
    grayscaleErodeFilter->Update();
  }

  catch ( itk::ExceptionObject & excep )
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    throw;
  }

  using ImageWriterType = itk::ImageFileWriter< ImageType >;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName( outputVolume );
  imageWriter->SetInput( grayscaleErodeFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
