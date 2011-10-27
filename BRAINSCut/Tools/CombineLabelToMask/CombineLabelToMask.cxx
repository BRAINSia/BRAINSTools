/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    BinaryThresholdImageFilter.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

     Author Eunyoung Regina Kim
     This is modified from itk Binary Threshold Image Filter.
=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkBinaryThresholdImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "CombineLabelToMaskCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef  double        InputPixelType;
  typedef  unsigned char OutputPixelType;

  const unsigned char Dim = 3;

  typedef itk::Image<InputPixelType,  Dim> InputImageType;
  typedef itk::Image<OutputPixelType, Dim> OutputImageType;

  typedef itk::BinaryThresholdImageFilter<
      InputImageType, OutputImageType>  FilterType;

  typedef itk::ImageFileReader<InputImageType>  ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  reader->SetFileName( inputVolume );

  filter->SetInput( reader->GetOutput() );

  filter->SetOutsideValue( outsideValue );
  filter->SetInsideValue(  insideValue  );

  filter->SetLowerThreshold( lowerThreshold );
  if( upperThreshold > 0 )
    {
    filter->SetUpperThreshold( upperThreshold );
    }

  filter->Update();

  writer->SetFileName( outputVolume );
  writer->Update();

  return EXIT_SUCCESS;
}
