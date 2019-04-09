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
// #include "Utilities.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include <itkIntensityWindowingImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "BRAINSThreadControl.h"
#include "GenerateBrainClippedImageCLP.h"
#include <BRAINSCommonLib.h>

int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  using PixelType = float;
  constexpr unsigned int Dim = 3;
  using ImageType = itk::Image<PixelType, Dim>;

  using ImageReaderType = itk::ImageFileReader<ImageType>;

  ImageReaderType::Pointer imgReader = ImageReaderType::New();
  ImageReaderType::Pointer mskReader = ImageReaderType::New();

  imgReader->SetFileName(inputImg);
  mskReader->SetFileName(inputMsk);

  using ImageMultiplyFilterType = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>;
  ImageMultiplyFilterType::Pointer imgMultiplyFilter =
    ImageMultiplyFilterType::New();

  imgMultiplyFilter->SetInput1( imgReader->GetOutput() );
  imgMultiplyFilter->SetInput2( mskReader->GetOutput() );

  // writer setting
  std::cout << "Writing output ... " << std::endl;
  using WriterType = itk::ImageFileWriter<ImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  std::cout << "* origin of input   ::" << imgReader->GetOutput()->GetOrigin()
            << std::endl
            << "* origin of output  ::" << imgMultiplyFilter->GetOutput()->GetOrigin()
            << std::endl;
  writer->SetFileName(outputFileName);
  writer->SetInput( imgMultiplyFilter->GetOutput() );
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

  return 0;
}
