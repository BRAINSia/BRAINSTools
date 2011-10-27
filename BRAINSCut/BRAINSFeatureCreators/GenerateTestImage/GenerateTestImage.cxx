#include <iostream>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "GenerateTestImageCLP.h"
int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  typedef double PixelType;
  const unsigned int Dim = 3;
  typedef  itk::Image<PixelType, Dim> ImageType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
  inputImageReader->SetFileName(inputVolume);

  // Resclaer
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType>
    RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();

  rescaler->SetInput( inputImageReader->GetOutput() );
  rescaler->SetOutputMinimum( lowerBoundOfOutputVolume);
  rescaler->SetOutputMaximum( upperBoundOfOutputVolume);

  try
    {
    rescaler->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }
  // Get Rescaler Output

  ImageType::Pointer rescaledImage = rescaler->GetOutput();

  // Set first/last value to be extreme
  typedef itk::ImageRegionIterator<ImageType> IteratorType;
  IteratorType inputIt( rescaledImage,  rescaledImage->GetLargestPossibleRegion() );

  inputIt.GoToBegin();
  inputIt.Set( 0 );

  inputIt.GoToEnd();
  --inputIt;
  inputIt.Set( 30000 );

  // write
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rescaledImage );
  writer->SetFileName( outputVolume );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }
}
