#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <itkThresholdImageFilter.h>
// Iterator
#include <itkImageRegionConstIterator.h>

int
main(int argc, char * *argv)
{
  // PARSE_ARGS;
  if( argc < 3 )
    {
    std::cout << " Usage:: Need classified Image to calculate " << std::endl
              << argv[0] << "   "
              << " [classified Image Name] [output name] "
              << std::endl;
    return 0;
    }

  typedef float PixelType;
  const unsigned int Dim = 3;
  typedef  itk::Image<PixelType, Dim> ImageType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  ImageReaderType::Pointer classified_imageReader = ImageReaderType::New();

  classified_imageReader->SetFileName(argv[1]);

  typedef itk::ThresholdImageFilter<ImageType> ThresholdFilterType;

  ThresholdFilterType::Pointer classified_gradientFilter = ThresholdFilterType::New();

  classified_gradientFilter->SetInput( classified_imageReader->GetOutput() );

  classified_gradientFilter->SetLower(50);
  classified_gradientFilter->Update();

  // writer setting
  std::cout << "Writing output ... " << std::endl;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  const char *outputFileName = argv[2];
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);

  rescaler->SetInput( classified_gradientFilter->GetOutput() );
  writer->SetFileName(outputFileName);
  writer->SetInput( rescaler->GetOutput() );
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
