#include "itkImage.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"

#include "GenerateTestImageCLP.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const unsigned int dimension = 3;
  typedef itk::Image<double, dimension>        InputImageType;
  typedef itk::Image<unsigned char, dimension> OutputImageType;
  // Create input image
  typedef itk::ImageFileReader<InputImageType> ImageReaderType;

  ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
  inputImageReader->SetFileName(inputVolume);
  inputImageReader->Update();

  InputImageType::Pointer inputImage = inputImageReader->GetOutput();

  InputImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();

  std::cout << "Input size: " << inputSize << std::endl;

  // Resclaer
  typedef itk::RescaleIntensityImageFilter<InputImageType, InputImageType> RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();

  rescaler->SetInput( inputImage );
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

  // Resize
  InputImageType::SizeType outputSize;
  outputSize.Fill(outputVolumeSize);
  InputImageType::SpacingType outputSpacing;
  for( unsigned int i = 0; i < dimension; i++ )
    {
    outputSpacing[i] = inputImage->GetSpacing()[i]
      * (static_cast<double>(inputSize[i]) / static_cast<double>(outputSize[i]) );
    }

  typedef itk::IdentityTransform<double, dimension>                TransformType;
  typedef itk::ResampleImageFilter<InputImageType, InputImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();

  resample->SetInput( rescaler->GetOutput() );
  resample->SetSize( outputSize);
  resample->SetOutputSpacing(outputSpacing);
  resample->SetTransform(TransformType::New() );
  resample->UpdateLargestPossibleRegion();
  resample->SetOutputDirection( inputImage->GetDirection() );
  resample->SetOutputOrigin( inputImage->GetOrigin() );

  std::cout << "Output size: " << resample->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

  // casting
  typedef itk::CastImageFilter<InputImageType, OutputImageType> CasterType;

  CasterType::Pointer caster = CasterType::New();

  caster->SetInput( resample->GetOutput() );

  // writing
  typedef itk::ImageFileWriter<InputImageType> WriterType;
  std::cout << "Writing output... " << std::endl;
  WriterType::Pointer outputWriter = WriterType::New();
  outputWriter->SetFileName( outputVolume );
  outputWriter->SetInput( resample->GetOutput() );
  outputWriter->Update();

  return EXIT_SUCCESS;
}
