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
/*
 * \author: Ali Ghayoor
 * at SINAPSE Lab,
 * The University of Iowa 2016
 */

#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDivideImageFilter.h>

#include "GenerateMaxGradientImage.h"

#include "GenerateEdgeMapImageCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const unsigned int Dim = 3;

  typedef itk::Image<float, Dim>                  FloatImageType;
  typedef FloatImageType::Pointer                 FloatImagePointer;
  typedef std::vector<FloatImagePointer>          InputImageList;
  typedef itk::ImageFileReader<FloatImageType>    ImageReaderType;
  typedef itk::Image<unsigned char, Dim>          CharImageType;
  typedef itk::ImageFileReader<CharImageType>     MaskReaderType;

  if( outputEdgeMap.compare( "" ) == 0 &&
      outputMaximumGradientImage.compare( "" ) == 0 )
    {
    std::cerr << "ERROR: No output image is specified!" << std::endl;
    return EXIT_FAILURE;
    }

  std::vector<std::string> inputMRFileNames;
  if( inputMRVolumes.size() > 0 )
    {
    inputMRFileNames = inputMRVolumes;
    }
  else
    {
    std::cerr << "ERROR: At least one MR image modality is needed "
              << "to generate the maximum gradient image and edgemap image."
              << std::endl;
    return EXIT_FAILURE;
    }
  const unsigned int numberOfMRImages = inputMRFileNames.size(); // number of modality images

  // Read the input MR modalities and set them in a vector of images
  typedef ImageReaderType::Pointer             LocalReaderPointer;

  InputImageList inputMRImageModalitiesList;
  for( unsigned int i = 0; i < numberOfMRImages; i++ )
    {
    std::cout << "Reading image: " << inputMRFileNames[i] << std::endl;
    LocalReaderPointer imgreader = ImageReaderType::New();
    imgreader->SetFileName( inputMRFileNames[i].c_str() );
    try
      {
      imgreader->Update();
      }
    catch( ... )
      {
      std::cerr << "ERROR:  Could not read image " << inputMRFileNames[i] << "." << std::endl;
      return EXIT_FAILURE;
      }
    inputMRImageModalitiesList.push_back( imgreader->GetOutput() );
    }

  CharImageType::Pointer mask;
  if( inputMask.compare( "" ) != 0 )
    {
    MaskReaderType::Pointer maskreader = MaskReaderType::New();
    maskreader->SetFileName( inputMask );
    try
      {
      maskreader->Update();
      }
    catch( ... )
      {
      std::cerr << "ERROR:  Could not read mask " << inputMask << "." << std::endl;
      return EXIT_FAILURE;
      }
    mask = maskreader->GetOutput();
    }

  // Create maximum gradient image using the input structural MR images.
  CharImageType::Pointer MGI = GenerateMaxGradientImage<
    FloatImageType,CharImageType,CharImageType>(inputMRImageModalitiesList,
                                                lowerPercentileMatching,
                                                upperPercentileMatching,
                                                minimumOutputRange,
                                                maximumOutputRange,
                                                mask);

  if( outputMaximumGradientImage.compare( "" ) != 0 )
    {
    typedef itk::ImageFileWriter<CharImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->UseCompressionOn();
    writer->SetFileName( outputMaximumGradientImage );
    writer->SetInput( MGI );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & exp )
      {
      std::cerr << "ExceptionObject with writer" << std::endl;
      std::cerr << exp << std::endl;
      return EXIT_FAILURE;
      }
    }

  // EdgeMap is created as the inverse of Maximum Gradient Image
  typedef itk::DivideImageFilter <FloatImageType, CharImageType, FloatImageType> DivideImageFilterType;

  DivideImageFilterType::Pointer divideImageFilter = DivideImageFilterType::New ();
  divideImageFilter->SetInput1( 1.0 );
  divideImageFilter->SetInput2( MGI );
  divideImageFilter->Update();
  FloatImageType::Pointer edgeMap = divideImageFilter->GetOutput();

  if( outputEdgeMap.compare( "" ) != 0 )
    {
    typedef itk::ImageFileWriter<FloatImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->UseCompressionOn();
    writer->SetFileName( outputEdgeMap );
    writer->SetInput( edgeMap );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & exp )
      {
      std::cerr << "ExceptionObject with writer" << std::endl;
      std::cerr << exp << std::endl;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}
