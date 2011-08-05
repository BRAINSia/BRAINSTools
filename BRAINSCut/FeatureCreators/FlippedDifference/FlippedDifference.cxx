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
#include "itkGrayscaleDilateImageFilter.h"
// #include "itkPluginFilterWatcher.h"
#include "itkMaskImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkConstrainedValueDifferenceImageFilter.h"
#include "itkOrientImageFilter.h"
#include "itkImageDuplicator.h"
#include "FlippedDifferenceCLP.h"

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

  //  typedef signed short       PixelType;
  typedef float PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType,  Dimension>   ImageType;
  typedef itk::ImageFileReader<ImageType>     ReaderType;
  typedef itk::Image<char,  Dimension>        MaskImageType;
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;

  ReaderType::Pointer     imageReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  imageReader->SetFileName( inputVolume.c_str() );
  maskReader->SetFileName( inputMaskVolume.c_str() );

  typedef itk::MaskImageFilter<ImageType, MaskImageType, ImageType> MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  typedef itk::ConstrainedValueDifferenceImageFilter<ImageType, ImageType, ImageType> DifferenceFilterType;
  DifferenceFilterType::Pointer differenceFilter = DifferenceFilterType::New();

  typedef itk::ImageDuplicator<ImageType> ImageDuplicatorType;
  ImageDuplicatorType::Pointer duplicateImageFilter = ImageDuplicatorType::New();

  try
    {
    maskReader->Update();

    maskFilter->SetInput1( imageReader->GetOutput() );
    maskFilter->SetInput2( maskReader->GetOutput() );
    maskFilter->SetOutsideValue(0);
    maskFilter->Update();

    duplicateImageFilter->SetInputImage( maskFilter->GetOutput() );
    duplicateImageFilter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << excep << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageRegionIterator<ImageType> ImageRegionIteratorType;
  ImageRegionIteratorType imgItr( maskFilter->GetOutput(), maskFilter->GetOutput()->GetRequestedRegion() );
  unsigned long           xMax = maskFilter->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  for( imgItr.GoToBegin(); !imgItr.IsAtEnd(); ++imgItr )
    {
    ImageType::IndexType idx = imgItr.GetIndex();
    ImageType::IndexType idxOpposite = imgItr.GetIndex();
    idxOpposite[0] = xMax - idx[0];
    if( maskReader->GetOutput()->GetPixel(idx) != 0 )
      {
      duplicateImageFilter->GetOutput()->SetPixel( idx, maskFilter->GetOutput()->GetPixel(idxOpposite) );
      }
    }

  try
    {
    differenceFilter->SetInput1( duplicateImageFilter->GetOutput() );
    differenceFilter->SetInput2( maskFilter->GetOutput() );
    differenceFilter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << excep << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter<ImageType> FileWriterType;
  FileWriterType::Pointer imageWriter = FileWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetInput( differenceFilter->GetOutput() );
  imageWriter->SetFileName( outputVolume.c_str() );
  imageWriter->Modified();
  imageWriter->Update();

  return EXIT_SUCCESS;
}
