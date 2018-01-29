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

#include "itkNoiseImageFilter.h"

#include "TextureFromNoiseImageFilterCLP.h"
#include <BRAINSCommonLib.h>

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  typedef float PixelType;
  constexpr unsigned int Dimension = 3;

  typedef itk::Image<PixelType,  Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>   ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();

  imageReader->SetFileName( inputVolume.c_str() );

  typedef itk::NoiseImageFilter<ImageType, ImageType> NoiseImageFilterType;
  NoiseImageFilterType::Pointer noiseFilter = NoiseImageFilterType::New();

  try
    {
    noiseFilter->SetInput( imageReader->GetOutput() );
    noiseFilter->SetRadius( inputRadius);
    noiseFilter->Update();
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    throw;
    }

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName(outputVolume);
  imageWriter->SetInput( noiseFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
