/*=========================================================================
 * Author : Eun Young Kim
 *
 * This is modifed from the ITK Exmaple as appeared below.
=========================================================================*/

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    GradientAnisotropicDiffusionImageFilter.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include "GradientAnisotropicDiffusionImageFilterCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // PARSE_ARG brings following:
  //  input/outputVolume, timeStep, conductance, numberOfIterations

  typedef    float InputPixelType;
  typedef    float OutputPixelType;
  const int Dimension = 3;

  typedef itk::Image<InputPixelType,  Dimension> InputImageType;
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

  typedef itk::ImageFileReader<InputImageType> ReaderType;

  typedef itk::GradientAnisotropicDiffusionImageFilter<
      InputImageType, OutputImageType>  FilterType;
  FilterType::Pointer filter = FilterType::New();

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume );

  filter->SetInput( reader->GetOutput() );

  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );

  filter->Update();

  //  The output of the filter is rescaled here and then sent to a writer.
  typedef unsigned char                         WritePixelType;
  typedef itk::Image<WritePixelType, Dimension> WriteImageType;
  typedef itk::RescaleIntensityImageFilter<
      OutputImageType, WriteImageType> RescaleFilterType;

  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );

  typedef itk::ImageFileWriter<WriteImageType> WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume );

  rescaler->SetInput( filter->GetOutput() );
  writer->SetInput( rescaler->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
