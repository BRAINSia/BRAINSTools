/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    MathematicalMorphologyBinaryFilters.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
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

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  const unsigned int Dimension = 3;

  typedef float         InputPixelType;
  typedef unsigned char OutputPixelType;

  typedef itk::Image<InputPixelType,  Dimension> InputImageType;
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

  typedef itk::ImageFileReader<InputImageType>  ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;

  typedef itk::BinaryThresholdImageFilter<InputImageType, OutputImageType> ThresholdFilterType;

  typedef itk::BinaryBallStructuringElement<
      InputPixelType,
      Dimension>             StructuringElementType;

  typedef itk::BinaryDilateImageFilter<
      OutputImageType,
      OutputImageType,
      StructuringElementType>  DilateFilterType;

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
