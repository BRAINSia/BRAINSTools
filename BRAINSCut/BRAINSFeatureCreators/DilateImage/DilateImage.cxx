#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkMaskImageFilter.h"
#include "DilateImageCLP.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( inputMaskVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputMaskVolume Required! "  << std::endl;
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

  typedef itk::Image<PixelType,  Dimension>   ImageType;
  typedef itk::ImageFileReader<ImageType>     ReaderType;
  typedef itk::Image<char,  Dimension>        MaskImageType;
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;

  ReaderType::Pointer     imageReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  ImageType::Pointer image = ImageType::New();

  imageReader->SetFileName( inputVolume.c_str() );
  maskReader->SetFileName( inputMaskVolume.c_str() );

  typedef itk::MaskImageFilter<ImageType, MaskImageType, ImageType> MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  typedef itk::BinaryBallStructuringElement<PixelType, Dimension> StructuringElementType;

  typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType, StructuringElementType> DilateFilterType;
  DilateFilterType::Pointer grayscaleDilateFilter = DilateFilterType::New();

  StructuringElementType structuringElement;
  structuringElement.SetRadius(inputRadius);
  structuringElement.CreateStructuringElement();
  grayscaleDilateFilter->SetKernel(structuringElement);

  try
    {
    maskFilter->SetInput1( imageReader->GetOutput() );
    maskFilter->SetInput2( maskReader->GetOutput() );
    maskFilter->SetOutsideValue(0);

    grayscaleDilateFilter->SetInput( maskFilter->GetOutput() );
    grayscaleDilateFilter->Update();
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    return EXIT_FAILURE;
    }

  // Output Dilated Image
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName(outputVolume);
  imageWriter->SetInput( grayscaleDilateFilter->GetOutput() );
  imageWriter->Update();
  return EXIT_SUCCESS;
}
