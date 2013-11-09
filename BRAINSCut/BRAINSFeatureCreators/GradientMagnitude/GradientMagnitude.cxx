#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMyGradientMagnitudeImageFilter.h"

#include "GradientMagnitudeCLP.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef float PixelType;
  // typedef unsigned long       PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType,  Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>   ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();

  imageReader->SetFileName( inputVolume.c_str() );

  typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientFilterType;
  GradientFilterType::Pointer gradientFilter = GradientFilterType::New();

  try
    {
    std::cout << __LINE__ << "::" << __FILE__ << std::endl;
    gradientFilter->SetInput( imageReader->GetOutput() );
    gradientFilter->Update();
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
  imageWriter->SetInput( gradientFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
