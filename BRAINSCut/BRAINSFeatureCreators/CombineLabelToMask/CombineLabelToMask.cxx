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

#include "itkBinaryThresholdImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "CombineLabelToMaskCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

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
