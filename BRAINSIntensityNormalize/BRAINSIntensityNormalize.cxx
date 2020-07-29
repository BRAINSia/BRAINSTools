//
// Created by Alexander Powers on 27-07-2020.
// Contact: alexander-powers@uiowa.edu
//
#include "BRAINSIntensityNormalizeCLP.h"

// ITK Imports
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include <BRAINSIntensityTransform.h>

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  constexpr unsigned int Dimension = 3;

  const bool clip = !no_clip;
  const bool relative = !no_relative;
  auto       usage_message = [&]() -> void {
    std::cout << "inputVolume is: " << inputVolume << std::endl;
    std::cout << "outputVolume is: " << outputVolume << std::endl;
    std::cout << "lowerPercentile is: " << lowerPercentile << std::endl;
    std::cout << "upperPercentile is: " << upperPercentile << std::endl;
    std::cout << "lowerOutputIntensity is: " << lowerOutputIntensity << std::endl;
    std::cout << "upperOutputIntensity is: " << upperOutputIntensity << std::endl;
    std::cout << "clip is: " << (clip ? "true" : "false") << std::endl;
    std::cout << "relative is: " << (relative ? "true" : "false") << std::endl;
  };
  usage_message();

  using InputPixelType = double;
  using InputImageType = itk::Image<InputPixelType, Dimension>;

  using OutputPixelType = unsigned short;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  InputImageType::Pointer input_image = [](const std::string & input_volume_fn) -> InputImageType::Pointer {
    using ReaderType = itk::ImageFileReader<InputImageType>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(input_volume_fn);
    reader->Update();
    return reader->GetOutput();
  }(inputVolume);

  auto output_image = brains_intensity_normalize_quantiles<InputImageType, OutputImageType>(
    input_image, lowerPercentile, upperPercentile, lowerOutputIntensity, upperOutputIntensity, clip, relative);
  // write
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputVolume);
  writer->SetInput(output_image);
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Writing Error:" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
