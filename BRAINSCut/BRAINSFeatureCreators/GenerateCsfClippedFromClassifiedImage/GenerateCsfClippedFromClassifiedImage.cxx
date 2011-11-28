#include <iostream>
#include <fstream>
#include <sstream>
#include <itkThresholdImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>

#include "GenerateCsfClippedFromClassifiedImageCLP.h"

int
main(int argc, char * *argv)
{
  PARSE_ARGS;

  typedef float PixelType;
  const unsigned int Dim = 3;
  typedef  itk::Image<PixelType, Dim> ImageType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  ImageReaderType::Pointer classified_imageReader = ImageReaderType::New();

  classified_imageReader->SetFileName( inputCassifiedVolume );

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
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);

  rescaler->SetInput( classified_gradientFilter->GetOutput() );
  writer->SetFileName( outputVolume );
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
