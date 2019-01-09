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
#include "itkGrayscaleDilateImageFilter.h"
// #include "itkPluginFilterWatcher.h"
#include "itkMaskImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkConstrainedValueDifferenceImageFilter.h"
#include "itkOrientImageFilter.h"
#include "itkImageDuplicator.h"
#include "FlippedDifferenceCLP.h"
#include <BRAINSCommonLib.h>

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
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

  //  using PixelType = signed short;
  using PixelType = float;
  constexpr unsigned int Dimension = 3;

  using ImageType = itk::Image<PixelType,  Dimension>;
  using ReaderType = itk::ImageFileReader<ImageType>;
  using MaskImageType = itk::Image<char,  Dimension>;
  using MaskReaderType = itk::ImageFileReader<MaskImageType>;

  ReaderType::Pointer     imageReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  imageReader->SetFileName( inputVolume.c_str() );
  maskReader->SetFileName( inputMaskVolume.c_str() );

  using MaskFilterType = itk::MaskImageFilter<ImageType, MaskImageType, ImageType>;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  using DifferenceFilterType = itk::ConstrainedValueDifferenceImageFilter<ImageType, ImageType, ImageType>;
  DifferenceFilterType::Pointer differenceFilter = DifferenceFilterType::New();

  using ImageDuplicatorType = itk::ImageDuplicator<ImageType>;
  ImageDuplicatorType::Pointer duplicateImageFilter = ImageDuplicatorType::New();

  try
    {
    maskReader->Update();

    maskFilter->SetInput1( imageReader->GetOutput() );
    maskFilter->SetInput2( maskReader->GetOutput() );
    maskFilter->SetOutsideValue(0);
    maskFilter->Update();

    duplicateImageFilter->SetInputImage( maskFilter->GetOutput() );
    duplicateImageFilter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << excep << std::endl;
    return EXIT_FAILURE;
    }

  using ImageRegionIteratorType = itk::ImageRegionIterator<ImageType>;
  ImageRegionIteratorType imgItr( maskFilter->GetOutput(), maskFilter->GetOutput()->GetRequestedRegion() );
  unsigned long           xMax = maskFilter->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  for( imgItr.GoToBegin(); !imgItr.IsAtEnd(); ++imgItr )
    {
    ImageType::IndexType idx = imgItr.GetIndex();
    ImageType::IndexType idxOpposite = imgItr.GetIndex();
    idxOpposite[0] = xMax - idx[0];
    if( maskReader->GetOutput()->GetPixel(idx) != 0 )
      {
      duplicateImageFilter->GetOutput()->SetPixel( idx, maskFilter->GetOutput()->GetPixel(idxOpposite) );
      }
    }

  try
    {
    differenceFilter->SetInput1( duplicateImageFilter->GetOutput() );
    differenceFilter->SetInput2( maskFilter->GetOutput() );
    differenceFilter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << excep << std::endl;
    return EXIT_FAILURE;
    }

  using FileWriterType = itk::ImageFileWriter<ImageType>;
  FileWriterType::Pointer imageWriter = FileWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetInput( differenceFilter->GetOutput() );
  imageWriter->SetFileName( outputVolume.c_str() );
  imageWriter->Modified();
  imageWriter->Update();

  return EXIT_SUCCESS;
}
