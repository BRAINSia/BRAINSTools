/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: SobelEdgeDetectorDetectionImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-16 21:52:48 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkSobelEdgeDetectionImageFilter.h"

#include "SobelEdgeDetectorCLP.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  typedef float CharPixelType;            //  IO
  typedef float RealPixelType;            //  Operations
  const unsigned int Dimension = 3;

  typedef itk::Image<CharPixelType, Dimension> CharImageType;
  typedef itk::Image<RealPixelType, Dimension> RealImageType;

  typedef itk::ImageFileReader<CharImageType> ReaderType;
  typedef itk::ImageFileWriter<CharImageType> WriterType;

  typedef itk::CastImageFilter<CharImageType, RealImageType> CastToRealFilterType;

  typedef itk::RescaleIntensityImageFilter<RealImageType, CharImageType> RescaleFilter;

  typedef itk::SobelEdgeDetectionImageFilter<RealImageType, RealImageType> CannyFilter;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
  RescaleFilter::Pointer        rescale = RescaleFilter::New();

  CannyFilter::Pointer sobelFilter = CannyFilter::New();

  reader->SetFileName(inputVolume);
  writer->SetFileName(outputVolume);

  // The output of an edge filter is 0 or 1
  rescale->SetOutputMinimum(0);
  rescale->SetOutputMaximum(255);

  toReal->SetInput( reader->GetOutput() );

  sobelFilter->SetInput( toReal->GetOutput() );

  rescale->SetInput( sobelFilter->GetOutput() );
  writer->SetInput( rescale->GetOutput() );
  writer->UseCompressionOn();

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
