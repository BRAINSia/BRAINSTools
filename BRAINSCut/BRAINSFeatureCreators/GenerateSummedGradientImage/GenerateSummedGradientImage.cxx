#include <iostream>
#include <fstream>
#include <sstream>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "BRAINSThreadControl.h"

#include "GenerateSummedGradientImageCLP.h"
#include <BRAINSCommonLib.h>
int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  typedef float PixelType;
  const unsigned int Dim = 3;
  typedef  itk::Image<PixelType, Dim> ImageType;

  typedef unsigned char                    OutputPixelType;
  typedef itk::Image<OutputPixelType, Dim> OutputImageType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  ImageReaderType::Pointer t1_imageReader = ImageReaderType::New();
  ImageReaderType::Pointer t2_imageReader = ImageReaderType::New();

  t1_imageReader->SetFileName(inputVolume1);
  t2_imageReader->SetFileName(inputVolume2);

  typedef itk::GradientMagnitudeImageFilter<ImageType,
                                            ImageType> GradientFilterType;

  GradientFilterType::Pointer t1_gradientFilter = GradientFilterType::New();
  GradientFilterType::Pointer t2_gradientFilter = GradientFilterType::New();

  t1_gradientFilter->SetInput( t1_imageReader->GetOutput() );
  t2_gradientFilter->SetInput( t2_imageReader->GetOutput() );
  t1_gradientFilter->Update();
  t2_gradientFilter->Update();

  typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxCalculatorType;
  MinMaxCalculatorType::Pointer t1_myMinMax = MinMaxCalculatorType::New();
  MinMaxCalculatorType::Pointer t2_myMinMax = MinMaxCalculatorType::New();
  t1_myMinMax->SetImage( t1_gradientFilter->GetOutput() );
  t2_myMinMax->SetImage( t2_gradientFilter->GetOutput() );
  t1_myMinMax->Compute();
  t2_myMinMax->Compute();

  typedef itk::Statistics::ScalarImageToHistogramGenerator<ImageType>
    HistogramGeneratorType;
  typedef HistogramGeneratorType::HistogramType
    HistogramType;
  HistogramGeneratorType::Pointer t1_HistogramGenerator =
    HistogramGeneratorType::New();
  t1_HistogramGenerator->SetInput( t1_gradientFilter->GetOutput() );
  t1_HistogramGenerator->SetNumberOfBins(1024);   // 4x oversampling to put into
                                                  // an unsigned char image.
  t1_HistogramGenerator->SetMarginalScale(10);

#if (ITK_VERSION_MAJOR < 4)
  // NOTHING TO DO IN ITKV3
#else
  t1_HistogramGenerator->SetHistogramMin( t1_myMinMax->GetMinimum() );
  t1_HistogramGenerator->SetHistogramMax( t1_myMinMax->GetMaximum() );
#endif
  t1_HistogramGenerator->Compute();
  HistogramType::ConstPointer t1_histogram = t1_HistogramGenerator->GetOutput();

  HistogramGeneratorType::Pointer t2_HistogramGenerator =
    HistogramGeneratorType::New();
  t2_HistogramGenerator->SetInput( t2_gradientFilter->GetOutput() );
  t2_HistogramGenerator->SetNumberOfBins(1024);   // 4x oversampling to put into
                                                  // an unsigned char image.
  t2_HistogramGenerator->SetMarginalScale(10);
#if (ITK_VERSION_MAJOR < 4)
  // NOTHING TO DO IN ITKV3
#else
  t2_HistogramGenerator->SetHistogramMin( t2_myMinMax->GetMinimum() );
  t2_HistogramGenerator->SetHistogramMax( t2_myMinMax->GetMaximum() );
#endif
  t2_HistogramGenerator->Compute();
  HistogramType::ConstPointer t2_histogram = t2_HistogramGenerator->GetOutput();

  const float UpperPercentileMatching = 0.95F; // Map 95th Quantile and above to
                                               // 127
  const float LowerPercentileMatching = 0.50F; // Map 25th Quantile and below to
                                               // 0
  // Now linear rescale gradient intensity values so that 10%=25/2 and 90%=230/2
  // so that when the two images are
  // added together they will fit into unsigned char
  typedef itk::IntensityWindowingImageFilter<ImageType,
                                             OutputImageType> RescalerType;
  RescalerType::Pointer t1_rescaler = RescalerType::New();
  t1_rescaler->SetInput( t1_gradientFilter->GetOutput() );
  t1_rescaler->SetOutputMinimum(0);
  t1_rescaler->SetOutputMaximum(127);
  const float T1Slope =
    ( t1_histogram->Quantile(0,
                             UpperPercentileMatching)
      - t1_histogram->Quantile(0, LowerPercentileMatching) ) / 80.0F;
  t1_rescaler->SetWindowMinimum( t1_histogram->Quantile(0,
                                                        LowerPercentileMatching) );
  t1_rescaler->SetWindowMaximum( t1_histogram->Quantile(0,
                                                        UpperPercentileMatching) + T1Slope * 100.0
                                 * ( 1.0 - UpperPercentileMatching ) );

  RescalerType::Pointer t2_rescaler = RescalerType::New();
  t2_rescaler->SetInput( t2_gradientFilter->GetOutput() );
  t2_rescaler->SetOutputMinimum(0);
  t2_rescaler->SetOutputMaximum(127);
  const float T2Slope =
    ( t2_histogram->Quantile(0,
                             UpperPercentileMatching)
      - t2_histogram->Quantile(0, LowerPercentileMatching) ) / 80.0F;
  t2_rescaler->SetWindowMinimum( t2_histogram->Quantile(0,
                                                        LowerPercentileMatching) );
  t2_rescaler->SetWindowMaximum( t2_histogram->Quantile(0,
                                                        UpperPercentileMatching) + T2Slope * 100.0
                                 * ( 1.0 - UpperPercentileMatching ) );

  // Type for Adder
  typedef itk::AddImageFilter<OutputImageType, OutputImageType,
                              OutputImageType> AddFilterType;
  AddFilterType::Pointer myAdder = AddFilterType::New();
  // Type for Maximum
  typedef itk::MaximumImageFilter<OutputImageType, OutputImageType,
                                  OutputImageType>
    MaximumFilterType;
  MaximumFilterType::Pointer myMax = MaximumFilterType::New();
  if( MaximumGradient == false )
    {
    // Set Output Image:: with Adder

    myAdder->SetInput1( t1_rescaler->GetOutput() );
    myAdder->SetInput2( t2_rescaler->GetOutput() );
    try
      {
      myAdder->Update();
      }
    catch( itk::ExceptionObject & exp )
      {
      std::cerr << "ExceptionObject with Iterator" << std::endl;
      std::cerr << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  else   // Set Output Image with Maximum of Gradient
    {
    myMax->SetInput1( t1_rescaler->GetOutput() );
    myMax->SetInput2( t2_rescaler->GetOutput() );
    try
      {
      myMax->Update();
      }
    catch( itk::ExceptionObject & exp )
      {
      std::cerr << "ExceptionObject with Iterator" << std::endl;
      std::cerr << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  // writer setting
  std::cout << "Writing output ... " << std::endl;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  typedef itk::IntensityWindowingImageFilter<OutputImageType,
                                             OutputImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);

  if( MaximumGradient == false )
    {
    rescaler->SetInput( myAdder->GetOutput() );
    }
  else
    {
    rescaler->SetInput( myMax->GetOutput() );
    }
  rescaler->Update();
  writer->SetFileName(outputFileName);
  writer->SetInput( rescaler->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "ExceptionObject with writer" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  return 0;
}
