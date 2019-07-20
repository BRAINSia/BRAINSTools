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
#if defined( _MSC_VER )
#  pragma warning( disable : 4786 )
#endif

#ifdef __BORLANDC__
#  define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include "GradientAnisotropicDiffusionImageFilterCLP.h"
#include <BRAINSCommonLib.h>

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // PARSE_ARG brings following:
  //  input/outputVolume, timeStep, conductance, numberOfIterations

  using InputPixelType = float;
  using OutputPixelType = float;
  constexpr int Dimension = 3;

  using InputImageType = itk::Image< InputPixelType, Dimension >;
  using OutputImageType = itk::Image< OutputPixelType, Dimension >;

  using ReaderType = itk::ImageFileReader< InputImageType >;

  using FilterType = itk::GradientAnisotropicDiffusionImageFilter< InputImageType, OutputImageType >;
  FilterType::Pointer filter = FilterType::New();

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume );

  filter->SetInput( reader->GetOutput() );

  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );

  filter->Update();

  //  The output of the filter is rescaled here and then sent to a writer.
  using WritePixelType = unsigned char;
  using WriteImageType = itk::Image< WritePixelType, Dimension >;
  using RescaleFilterType = itk::RescaleIntensityImageFilter< OutputImageType, WriteImageType >;

  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 255 );

  using WriterType = itk::ImageFileWriter< WriteImageType >;

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume );

  rescaler->SetInput( filter->GetOutput() );
  writer->SetInput( rescaler->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
