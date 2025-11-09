// \author Hans J. Johnson
#ifndef BRAINSIntensityTransform_h__
#define BRAINSIntensityTransform_h__

#include "itkHistogram.h"
#include "itkImageToHistogramFilter.h"
#include "itkStatisticsImageFilter.h"
#include <itkUnaryGeneratorImageFilter.h>

/**
 * Intensity normalize based on intensity intensity percentiles
 * @tparam InputImageType  Pixel type of input image
 * @tparam OutputImageType Pixel type of output image
 * @param input_image the input image
 * @param lowerPercentile  Reference point histogram percentile from input image
 * @param upperPercentile  Reference point histogram percentile from input image
 * @param lowerOutputIntensity Reference value for output image
 * @param upperOutputIntensity Reference value for output image
 * @param no_clip Do not clip values outside of (lowerOutputIntensity, upperOutputIntensity)
 * @param no_relative If true transfer function passes through points
 *  (lowerPercentile, lowerOutputIntensity) & (upperPercentile,upperOutputIntensity)
 *   else
 *   output_range=upperOutputIntensity-lowerOutputIntensity
 *   (lowerPercentile, lowerOutputIntensity+lowerPercentile*output_range  )
 *         & (upperPercentile,lowerOutputIntensity + uppperPercentile*output_range)
 * @return
 */
template <typename InputImageType, typename OutputImageType>
typename OutputImageType::Pointer
brains_intensity_normalize_quantiles(typename InputImageType::Pointer input_image,
                                     const double                     lowerPercentile,
                                     const double                     upperPercentile,
                                     const double                     lowerOutputIntensity,
                                     const double                     upperOutputIntensity,
                                     const bool                       clip,
                                     const bool                       relative,
                                     const bool                       print_diagnostics = false)
{
  // input validation
  if (lowerPercentile >= upperPercentile)
  {
    std::cerr << "ERROR: --upperPercentile must be greater than --lowerPercentile" << std::endl;
    return nullptr;
  }
  if (lowerOutputIntensity >= upperOutputIntensity)
  {
    std::cerr << "ERROR: --lowerOutputIntensity must be greater than --upperOutputIntensity" << std::endl;
    return nullptr;
  }

  // percentile based rescaling
  // find min/max pixels for image
  using StatisticsFilterType = itk::StatisticsImageFilter<InputImageType>;

  typename StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
  statsFilter->SetInput(input_image);
  try
  {
    statsFilter->Update();
  }
  catch (const itk::ExceptionObject & error)
  {
    std::cerr << "Stats Filter Error: " << error << std::endl;
    return nullptr;
  }
  const typename InputImageType::PixelType minPixel(statsFilter->GetMinimum());
  const typename InputImageType::PixelType maxPixel(statsFilter->GetMaximum());
  if (print_diagnostics)
  {
    const double variance = statsFilter->GetSigma();
    const double mean = statsFilter->GetMean();
    std::cout << "Input Image Min/Max/Mean/Variance = " << minPixel << "/" << maxPixel << "/" << mean << "/" << variance
              << std::endl;
  }
  // generate histogram to find quantiles
  constexpr unsigned int MeasurementVectorSize = 1; // Grayscale --> non-vector images
  const auto             binsPerDimension = static_cast<unsigned int>(1000);

  using ImageToHistogramFilterType = itk::Statistics::ImageToHistogramFilter<InputImageType>;
  typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType lowerBound(binsPerDimension);
  lowerBound.Fill(minPixel);

  typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType upperBound(binsPerDimension);
  upperBound.Fill(maxPixel);

  typename ImageToHistogramFilterType::HistogramType::SizeType size(MeasurementVectorSize);
  size.Fill(binsPerDimension);

  typename ImageToHistogramFilterType::Pointer imageToHistogramFilter = ImageToHistogramFilterType::New();
  imageToHistogramFilter->SetInput(input_image);
  imageToHistogramFilter->SetHistogramBinMinimum(lowerBound);
  imageToHistogramFilter->SetHistogramBinMaximum(upperBound);
  imageToHistogramFilter->SetHistogramSize(size);

  try
  {
    imageToHistogramFilter->Update();
  }
  catch (const itk::ExceptionObject & error)
  {
    std::cerr << "Histogram Error: " << error << std::endl;
    return nullptr;
  }

  typename ImageToHistogramFilterType::HistogramType::Pointer histogram = imageToHistogramFilter->GetOutput();
  const double lowerPercentileValue = histogram->Quantile(0, lowerPercentile);
  const double upperPercentileValue = histogram->Quantile(0, upperPercentile);

  // Compute relative if neccesary
  double relativeLowerOutputIntensity = lowerOutputIntensity;
  double relativeUpperOutputIntensity = upperOutputIntensity;
  if (relative)
  {
    double difference = upperOutputIntensity - lowerOutputIntensity;
    relativeLowerOutputIntensity = (difference * lowerPercentile) + lowerOutputIntensity;
    relativeUpperOutputIntensity = (difference * upperPercentile) + lowerOutputIntensity;
  }
  if (print_diagnostics)
  {
    std::cout << "Rescaling linearly from (" << lowerPercentileValue << ", " << upperPercentileValue << ") to ("
              << relativeLowerOutputIntensity << ", " << relativeUpperOutputIntensity << ")"
              << " clipping at (" << lowerOutputIntensity << ", " << upperOutputIntensity << ")" << std::endl;
  }

  const double slope =
    (relativeUpperOutputIntensity - relativeLowerOutputIntensity) / (upperPercentileValue - lowerPercentileValue);
  const double intercept = relativeUpperOutputIntensity - slope * upperPercentileValue;

  typename OutputImageType::Pointer output_image;
  using IntensityScalingFilter = typename itk::UnaryGeneratorImageFilter<InputImageType, OutputImageType>;
  typename IntensityScalingFilter::Pointer filter = IntensityScalingFilter::New();

  if (clip)
  {
    filter->SetFunctor([slope, intercept, upperOutputIntensity, lowerOutputIntensity](
                         const typename InputImageType::PixelType input_pixel_value_orig) ->
                       typename OutputImageType::PixelType {
                         const auto input_pixel_value = static_cast<double>(input_pixel_value_orig);
                         const auto temp = input_pixel_value * slope + intercept;
                         const auto return_value = (temp > upperOutputIntensity)
                                                     ? upperOutputIntensity
                                                     : ((temp < lowerOutputIntensity) ? lowerOutputIntensity : temp);
                         return static_cast<typename OutputImageType::PixelType>(return_value);
                       });
  }
  else
  {
    filter->SetFunctor([slope, intercept](const typename InputImageType::PixelType input_pixel_value_orig) ->
                       typename OutputImageType::PixelType {
                         const auto input_pixel_value = static_cast<double>(input_pixel_value_orig);
                         const auto temp = input_pixel_value * slope + intercept;
                         return static_cast<typename OutputImageType::PixelType>(temp);
                       });
  }
  filter->SetInput(input_image);
  filter->Update();
  output_image = filter->GetOutput();

  if (print_diagnostics)
  {
    using OutStatisticsFilterType = itk::StatisticsImageFilter<OutputImageType>;
    typename OutStatisticsFilterType::Pointer outStatsFilter = OutStatisticsFilterType::New();
    outStatsFilter->SetInput(output_image);
    try
    {
      outStatsFilter->Update();
    }
    catch (const itk::ExceptionObject & error)
    {
      std::cerr << "Stats Filter Error: " << error << std::endl;
      return nullptr;
    }
    const typename InputImageType::PixelType outMinPixel(outStatsFilter->GetMinimum());
    const typename InputImageType::PixelType outMaxPixel(outStatsFilter->GetMaximum());
    std::cout << "Output min/max: (" << outMinPixel << ", " << outMaxPixel << ")" << std::endl;
  }
  return output_image;
}

#endif // BRAINSIntensityTransform_h__
