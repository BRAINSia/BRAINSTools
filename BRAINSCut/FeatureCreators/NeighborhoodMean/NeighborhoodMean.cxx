#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"
#include "itkMaskImageFilter.h"
#include "itkMeanImageFilter.h"
#include "NeighborhoodMeanCLP.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;

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
    exit(1);
    }

  typedef float PixelType;
  // typedef unsigned long       PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType,  Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>   ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();

  imageReader->SetFileName( inputVolume.c_str() );

  typedef itk::MeanImageFilter<ImageType, ImageType> MeanFilterType;
  MeanFilterType::Pointer meanFilter = MeanFilterType::New();

  MeanFilterType::InputSizeType radius;
  radius[0] = inputRadius;
  radius[1] = inputRadius;
  radius[2] = inputRadius;

  try
    {
    meanFilter->SetInput( imageReader->GetOutput() );
    meanFilter->SetRadius(radius);
    meanFilter->Update();
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName(outputVolume);
  imageWriter->SetInput( meanFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
