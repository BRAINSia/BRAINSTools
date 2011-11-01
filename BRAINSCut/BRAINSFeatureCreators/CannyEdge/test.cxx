#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include "itkMaskImageFilter.h"
#include "CannyEdgeCLP.h"

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
    return EXIT_FAILURE;
    }

  const unsigned int dimension = 3;
  typedef float                                PixelType;
  typedef itk::Image<float, dimension>         InputImage;
  typedef itk::Image<unsigned char, dimension> OutputImage;

  itk::ImageFileReader<InputImage>::Pointer input = itk::ImageFileReader<InputImage>::New();
  input->SetFileName(inputVolume);

  // Set up filter
  itk::CannyEdgeDetectionImageFilter<InputImage, InputImage>::Pointer
    filter =
    itk::CannyEdgeDetectionImageFilter<InputImage, InputImage>::New();
  itk::SimpleFilterWatcher watcher(filter);
  filter->SetInput( input->GetOutput() );
  filter->SetUpperThreshold(30);
  filter->SetLowerThreshold(10);
  filter->SetThreshold(30);
  filter->SetVariance(1.0f);
  filter->SetMaximumError(.01f);

  itk::RescaleIntensityImageFilter<InputImage, OutputImage>::Pointer
    rescale =
    itk::RescaleIntensityImageFilter<InputImage, OutputImage>::New();
  rescale->SetInput( filter->GetOutput() );
  rescale->SetOutputMinimum(0);
  rescale->SetOutputMaximum(255);

  try
    {
    // Generate test image
    itk::ImageFileWriter<OutputImage>::Pointer writer;
    writer = itk::ImageFileWriter<OutputImage>::New();
    writer->UseCompressionOn();
    writer->SetInput( rescale->GetOutput() );
    writer->SetFileName(outputVolume);
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    ( &err )->Print(std::cerr);
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
