#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <itkIntensityWindowingImageFilter.h>

#include <itkMultiplyImageFilter.h>
#include "BRAINSThreadControl.h"

#include <GenerateBrainClippedImageCLP.h>
int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);
  typedef float PixelType;
  const unsigned int Dim = 3;
  typedef  itk::Image<PixelType, Dim> ImageType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  ImageReaderType::Pointer imgReader = ImageReaderType::New();
  ImageReaderType::Pointer mskReader = ImageReaderType::New();

  imgReader->SetFileName(inputImg);
  mskReader->SetFileName(inputMsk);

  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType>
    ImageMultiplyFilterType;
  ImageMultiplyFilterType::Pointer imgMultiplyFilter =
    ImageMultiplyFilterType::New();

  imgMultiplyFilter->SetInput1( imgReader->GetOutput() );
  imgMultiplyFilter->SetInput2( mskReader->GetOutput() );

  // writer setting
  std::cout << "Writing output ... " << std::endl;
  typedef itk::ImageFileWriter<ImageType> WriterType;
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
