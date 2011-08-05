#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include "itkMaskImageFilter.h"
#include "DistanceMapsCLP.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  bool violated = false;
  if( inputLabelVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputLabelVolume Required! "  << std::endl;
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

  typedef itk::Image<PixelType,  Dimension>    LabelImageType;
  typedef itk::Image<char,  Dimension>         MaskImageType;
  typedef itk::ImageFileReader<LabelImageType> ReaderType;
  typedef itk::ImageFileReader<MaskImageType>  MaskReaderType;

  ReaderType::Pointer     labelReader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  labelReader->SetFileName( inputLabelVolume.c_str() );
  maskReader->SetFileName( inputMaskVolume.c_str() );

  typedef itk::MaskImageFilter<LabelImageType, MaskImageType, LabelImageType> MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> BinaryThresholdFilterType;
  BinaryThresholdFilterType::Pointer binaryFilter = BinaryThresholdFilterType::New();

  typedef itk::SignedMaurerDistanceMapImageFilter<LabelImageType, LabelImageType> DistanceMapFilterType;
  DistanceMapFilterType::Pointer distanceMapFilter = DistanceMapFilterType::New();
  distanceMapFilter->SetInsideIsPositive(true); // Makes all distances positive
  typedef DistanceMapFilterType::OutputImageType DistanceMapImageType;

  try
    {
    maskFilter->SetInput1( labelReader->GetOutput() );
    maskFilter->SetInput2( maskReader->GetOutput() );
    maskFilter->SetOutsideValue(0);
    // binaryFilter->SetInput(labelReader->GetOutput());
    binaryFilter->SetInput( maskFilter->GetOutput() );
    binaryFilter->SetLowerThreshold(inputTissueLabel);
    binaryFilter->SetUpperThreshold(inputTissueLabel);
    binaryFilter->Update();
    distanceMapFilter->SetInput( binaryFilter->GetOutput() );
    distanceMapFilter->UseImageSpacingOff();
    distanceMapFilter->SquaredDistanceOff();
    distanceMapFilter->Update();
    maskReader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter<DistanceMapFilterType::OutputImageType> ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetFileName(outputVolume);
  imageWriter->SetInput( distanceMapFilter->GetOutput() );
  imageWriter->Update();

  return EXIT_SUCCESS;
}
