#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkSpeckleNoiseImageFilter.h"
#include "itkSaltAndPepperNoiseImageFilter.h"
#include "itkShotNoiseImageFilter.h"
#include "NoiseGeneratorCLP.h"
#include <BRAINSCommonLib.h>

int main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if( inputVolume == "" )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( outputVolume == "" )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  const int DIMENSION = 3;

  /** input iamges */
  typedef double                                 ReadInPixelType;
  typedef itk::Image<ReadInPixelType, DIMENSION> ReadInImageType;

  typedef itk::ImageFileReader<ReadInImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume );

  /** rescaler */
  typedef itk::RescaleIntensityImageFilter<ReadInImageType, ReadInImageType> RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();

  rescaler->SetInput( reader->GetOutput() );
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);

  /** casting */
  typedef unsigned char                    PixelType;
  typedef itk::Image<PixelType, DIMENSION> ImageType;

  typedef itk::CastImageFilter<ReadInImageType, ImageType> CastingType;
  CastingType::Pointer caster = CastingType::New();
  caster->SetInput( rescaler->GetOutput() );

  /** output image */
  ImageType::Pointer outputImage;

  if(  inputGaussianStandardDeviation < 0.0F &&
       inputShotNoiseScale < 0.0F  &&
       inputSpeckleNoiseStandardDeviation < 0.0F &&
       inputSaltAndPepperProbability < 0.0F )
    {
    std::cout << "ERROR:: No method is given " << std::endl;
    exit(EXIT_FAILURE);
    }

  /** gaussian */
  if( inputGaussianStandardDeviation > 0 )
    {
    typedef itk::AdditiveGaussianNoiseImageFilter<ImageType, ImageType> GaussianNoiseFilterType;
    GaussianNoiseFilterType::Pointer gaussianFilter = GaussianNoiseFilterType::New();
    gaussianFilter->SetInput( caster->GetOutput() );
    gaussianFilter->SetStandardDeviation( inputGaussianStandardDeviation );
    gaussianFilter->SetMean( inputGaussianMean );
    gaussianFilter->Update();
    outputImage = gaussianFilter->GetOutput();
    }

  /** shot noise */
  if( inputShotNoiseScale > 0 )
    {
    typedef itk::ShotNoiseImageFilter<ImageType, ImageType> ShotNoiseFilterType;
    ShotNoiseFilterType::Pointer shotNoiseFilter = ShotNoiseFilterType::New();
    shotNoiseFilter->SetInput( caster->GetOutput() );
    shotNoiseFilter->SetScale( inputShotNoiseScale );
    shotNoiseFilter->Update();
    outputImage = shotNoiseFilter->GetOutput();
    }

  /** speckle noise image filter */
  if( inputSpeckleNoiseStandardDeviation > 0 )
    {
    typedef itk::SpeckleNoiseImageFilter<ImageType, ImageType> SpeckleNoiseFilterType;
    SpeckleNoiseFilterType::Pointer speckleNoiseFilter = SpeckleNoiseFilterType::New();
    speckleNoiseFilter->SetInput( caster->GetOutput() );
    speckleNoiseFilter->SetStandardDeviation( inputSpeckleNoiseStandardDeviation );
    speckleNoiseFilter->Update();
    outputImage = speckleNoiseFilter->GetOutput();
    }

  /** salt and pepper noise image filter */
  if( inputSaltAndPepperProbability > 0 )
    {
    typedef itk::SaltAndPepperNoiseImageFilter<ImageType, ImageType> SaltAndPepperNoiseFilterType;
    SaltAndPepperNoiseFilterType::Pointer spFilter = SaltAndPepperNoiseFilterType::New();
    spFilter->SetInput( caster->GetOutput() );
    spFilter->SetProbability( inputSaltAndPepperProbability );
    spFilter->Update();
    outputImage = spFilter->GetOutput();
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(  outputImage );
  writer->SetFileName( outputVolume );
  writer->Update();

  return 0;
}
