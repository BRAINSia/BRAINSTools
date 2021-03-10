/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
//
// This is a simple program for [mattes mutual information/ means square error]
// metric evaluation over two input images and a BSpline transform.
//
//
#include "PerformMetricTestCLP.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkBSplineTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkTransformFileReader.h"


#define CHECK_PARAMETER_IS_SET(parameter, message)                                                                     \
  if (parameter == "")                                                                                                 \
  {                                                                                                                    \
    std::cerr << message << std::endl;                                                                                 \
    return EXIT_FAILURE;                                                                                               \
  }

int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  CHECK_PARAMETER_IS_SET(inputFixedImage, "Missing inputFixedImage parameter");
  CHECK_PARAMETER_IS_SET(inputMovingImage, "Missing inputMovingImage parameter");

  constexpr unsigned int Dimension = 3;
  using PixelType = float;

  using FixedImageType = itk::Image<PixelType, Dimension>;
  using MovingImageType = itk::Image<PixelType, Dimension>;

  using FixedImageReaderType = itk::ImageFileReader<FixedImageType>;
  using MovingImageReaderType = itk::ImageFileReader<MovingImageType>;

  using TransformType = itk::BSplineTransform<PixelType, Dimension, 3>;

  using MIMetricType =
    itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, PixelType>;
  using MSEMetricType =
    itk::MeanSquaresImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, PixelType>;
  using GenericMetricType = itk::ImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, PixelType>;
  using MetricSamplePointSetType = GenericMetricType::FixedSampledPointSetType;

  // Read input images
  FixedImageReaderType::Pointer fixedReader = FixedImageReaderType::New();
  fixedReader->SetFileName(inputFixedImage);
  FixedImageType::Pointer fixedImage = fixedReader->GetOutput();
  fixedReader->Update();

  MovingImageReaderType::Pointer movingReader = MovingImageReaderType::New();
  movingReader->SetFileName(inputMovingImage);
  MovingImageType::Pointer movingImage = movingReader->GetOutput();
  movingReader->Update();

  // Set input transform
  TransformType::Pointer transform = TransformType::New();
  if (!inputBSplineTransform.empty())
  {
    std::cout << "Read transform file from the disk ..." << std::endl;
    itk::TransformFileReader::Pointer transReader = itk::TransformFileReader::New();
    transReader->SetFileName(inputBSplineTransform);
    try
    {
      transReader->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception caught!" << std::endl;
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
    }
    transform = dynamic_cast<TransformType *>(transReader->GetTransformList()->begin()->GetPointer());
    if (transform.IsNull())
    {
      std::cerr << "ERROR: Needed a BSpline tranform as the input transform file." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cout << "No input transform file ... transform file is set to identity" << std::endl;
    transform->SetIdentity();
  }

  // Set metric
  GenericMetricType::Pointer metric;
  if (metricType == "MMI")
  {
    MIMetricType::Pointer mattesMetric = MIMetricType::New();
    mattesMetric->SetNumberOfHistogramBins(numberOfHistogramBins);
    mattesMetric->SetUseFixedImageGradientFilter(false);
    mattesMetric->SetUseMovingImageGradientFilter(false);
    metric = mattesMetric;
  }
  else if (metricType == "MSE")
  {
    MSEMetricType::Pointer msqMetric = MSEMetricType::New();
    metric = msqMetric;
  }
  else
  {
    std::cerr << "Error: Invalid parameter for metric type!" << std::endl;
    return EXIT_FAILURE;
  }

  metric->SetVirtualDomainFromImage(fixedImage);
  metric->SetFixedImage(fixedImage);
  metric->SetMovingImage(movingImage);
  metric->SetMovingTransform(transform);

  //  REGULAR sampling
  std::cout << "Use regular sampling ..." << std::endl;

  if (numberOfSamples == 0)
  {
    std::cerr << "Error: Valid number of samples needed!" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned long numberOfAllSamples = fixedImage->GetBufferedRegion().GetNumberOfPixels();
  const double        samplingPercentage = static_cast<double>(numberOfSamples) / numberOfAllSamples;
  std::cout << "Number of Requesete Samples: " << numberOfSamples << std::endl;
  std::cout << "Number of All Samples: " << numberOfAllSamples << std::endl;
  std::cout << "Sampling Percentage: " << samplingPercentage << std::endl;

  MetricSamplePointSetType::Pointer samplePointSet = MetricSamplePointSetType::New();
  samplePointSet->Initialize();
  using SamplePointType = MetricSamplePointSetType::PointType;
  unsigned long index = 0;

  const auto sampleCount = static_cast<unsigned long>(std::ceil(1.0 / samplingPercentage));
  // std::cout << "sample count: " << sampleCount << std::endl;
  unsigned long                                          count = sampleCount;
  itk::ImageRegionConstIteratorWithIndex<FixedImageType> It(fixedImage, fixedImage->GetBufferedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if (count == sampleCount)
    {
      count = 0; // Reset counter
      SamplePointType point;
      fixedImage->TransformIndexToPhysicalPoint(It.GetIndex(), point);
      samplePointSet->SetPoint(index, point);
      ++index;
    }
    ++count;
  }

  metric->SetFixedSampledPointSet(samplePointSet);
  metric->SetUseSampledPointSet(true);

  metric->Initialize();

  GenericMetricType::MeasureType    measure = NAN;
  GenericMetricType::DerivativeType derivative(metric->GetTransform()->GetNumberOfParameters());
  try
  {
    metric->GetValueAndDerivative(measure, derivative);
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Number of Valid Points: " << metric->GetNumberOfValidPoints() << " of "
            << metric->GetVirtualRegion().GetNumberOfPixels() << std::endl;
  std::cout << "measure: " << measure << std::endl;
  //  std::cout << "Derivatives: " << derivative << std::endl;

  return EXIT_SUCCESS;
}
