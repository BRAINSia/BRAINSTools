#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkNoiseImageFilter.h"

#include "TextureFromNoiseImageFilterCLP.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  typedef float PixelType;
  const unsigned int Dimension = 3;

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
    throw excep;
    }

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName(outputVolume);
  imageWriter->SetInput( noiseFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
