#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "NoiseGeneratorCLP.h"

int main(int argc, char * argv[])
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  PARSE_ARGS;
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

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

  const int dim = 3;

  typedef float                      PixelType;
  typedef itk::Image<PixelType, dim> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume );

  typedef itk::AdditiveGaussianNoiseImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetStandardDeviation( inputGaussianStandardDeviation );

  filter->SetMean( inputGaussianMean );
  filter->Update();

  itk::SimpleFilterWatcher watcher(filter, "filter");

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  writer->SetFileName( outputVolume );
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  writer->Update();

  return 0;
}
