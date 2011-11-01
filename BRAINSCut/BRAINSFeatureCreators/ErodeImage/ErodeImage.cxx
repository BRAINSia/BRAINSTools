#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMeanImageFilter.h"
#include "ErodeImageCLP.h"

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
  if( inputMaskVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputMaskVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef float PixelType;
  // typedef unsigned long       PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType,  Dimension>   ImageType;
  typedef itk::Image<char,  Dimension>        MaskImageType;
  typedef itk::ImageFileReader<ImageType>     ReaderType;
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;

  ReaderType::Pointer     imageReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  imageReader->SetFileName( inputVolume.c_str() );
  maskReader->SetFileName( inputMaskVolume.c_str() );

  typedef itk::MaskImageFilter<ImageType, MaskImageType, ImageType> MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  typedef itk::BinaryBallStructuringElement<PixelType, Dimension> StructuringElementType;

  typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType, StructuringElementType> ErodeFilterType;
  ErodeFilterType::Pointer grayscaleErodeFilter = ErodeFilterType::New();

  StructuringElementType structuringElement;
  structuringElement.SetRadius(inputRadius);
  structuringElement.CreateStructuringElement();
  grayscaleErodeFilter->SetKernel(structuringElement);

  try
    {
    maskFilter->SetInput1( imageReader->GetOutput() );
    maskFilter->SetInput2( maskReader->GetOutput() );
    maskFilter->SetOutsideValue(0);

    grayscaleErodeFilter->SetInput( maskFilter->GetOutput() );
    grayscaleErodeFilter->Update();
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
  imageWriter->SetInput( grayscaleErodeFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
