#ifndef __StandardizeMaskIntensity_h
#define __StandardizeMaskIntensity_h
#include "itkMinimumMaximumImageCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkHistogram.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkNumericTraits.h"

#define MAX_IMAGE_OUTPUT_VALUE 4096

/* * * * *
 *   StandardizeMaskIntensity trims the upper and lower fractions of the histogram inside the mask
 *   based on tail size fractions uFract and lFract, respectively.  Then, using these histogram-based
 *   upper and lower bounds, it rescales the signal interval they represent to an interval marked by
 *   the target upper and lower values, and clips the output image signal to an interval that encloses
 *   the target interval.
 *
 *   Uses Review Statistics Histogram's Quantile method.  IntensityWindowingImageFilter clips to the
 *   same bounds it scales to as lowerQuantileValue target, so lowerQuantileValue little univariate extrapolation was done (see block comment).
 *
 * * * * */
template <class ImageType, class LabelImageType>
typename ImageType::Pointer StandardizeMaskIntensity(
  typename ImageType::Pointer image,
  typename LabelImageType::Pointer mask,
  const double lFract,
  const double uFract,
  const typename ImageType::PixelType lowerPeggedValue,
  const typename ImageType::PixelType upperPeggedValue,
  const typename ImageType::PixelType clipMin,
  const typename ImageType::PixelType clipMax)
{
  // the asserts, if you want them, go like this:
  // ImageType::PixelType is scalar.
  // 0.0 <= lFract < uFract <= 1.0
  // clipMin <= lowerPeggedValue < upperPeggedValue <= clipMax

  typedef typename itk::MinimumMaximumImageCalculator<ImageType> MinimumMaximumImageCalculator;
  typename MinimumMaximumImageCalculator::Pointer wholeStatistics = MinimumMaximumImageCalculator::New();
  wholeStatistics->SetImage(image);
  wholeStatistics->Compute();
  typename ImageType::PixelType imgMin = wholeStatistics->GetMinimum();
  typename ImageType::PixelType imgMax = wholeStatistics->GetMaximum();

  int numBins = vnl_math_rnd(imgMax - imgMin + 1);
  if( numBins < 256 )
    {
    numBins = 256;
    }
  if( numBins > MAX_IMAGE_OUTPUT_VALUE )
    {
    numBins = MAX_IMAGE_OUTPUT_VALUE;
    }

  const typename LabelImageType::PixelType maskInteriorLabel = 1;
  typename LabelImageType::Pointer internalMask;
  if( mask.IsNull() )
    {
    internalMask = LabelImageType::New();
    internalMask->CopyInformation(image);
    internalMask->SetRegions( image->GetLargestPossibleRegion() );
    internalMask->Allocate();
    internalMask->FillBuffer(maskInteriorLabel);
    }
  else
    {
    typename itk::ThresholdImageFilter<LabelImageType>::Pointer thresholdFilter
      = itk::ThresholdImageFilter<LabelImageType>::New();

    thresholdFilter->SetInput(mask);
    thresholdFilter->ThresholdAbove(1); // Values less than or equal to are set
                                        // to OutsideValue
    thresholdFilter->SetOutsideValue(maskInteriorLabel);
    thresholdFilter->Update();
    internalMask = thresholdFilter->GetOutput();
    }

  typedef typename itk::LabelStatisticsImageFilter<ImageType, LabelImageType> LabelStatisticsImageFilter;
  typename LabelStatisticsImageFilter::Pointer maskedStatistics = LabelStatisticsImageFilter::New();
  maskedStatistics->SetInput(image); // i.clipMin., image.
  maskedStatistics->SetLabelInput(internalMask);
  maskedStatistics->UseHistogramsOn();
  maskedStatistics->SetHistogramParameters(numBins, imgMin, imgMax);
  maskedStatistics->Update();
  typename LabelStatisticsImageFilter::HistogramType::Pointer hist =
    maskedStatistics->GetHistogram( maskInteriorLabel );
  if( hist.IsNull() )
    {
    itkGenericExceptionMacro("histogram had no value for label "
                             << maskInteriorLabel);
    }

  // Remark:  Since itk's filter uses the same bounds to scale and to clip,
  // and since the clipping bounds are outside (surrounding) the scaling
  // bounds, we need to compute adjustedWindowMaxForClipping and
  // adjustedWindowMinForClipping
  // to fake out ITK developers and get both
  // jobs done at once.  Fortunately, we may use itkHistogram's Quantile
  // routine:
  typedef typename itk::IntensityWindowingImageFilter<ImageType> IntensityWindowingImageFilter;
  typename IntensityWindowingImageFilter::Pointer intensityMapper = IntensityWindowingImageFilter::New();
  intensityMapper->SetInput( maskedStatistics->GetOutput() ); // i.clipMin.,
                                                              // image.
  // NOTE:  The math below is to extend the range to the clipping region.
  //
  //
  typedef typename itk::NumericTraits<typename ImageType::PixelType>::ScalarRealType Real;
  const Real lowerQuantileValue = hist->Quantile(0, lFract);
  const Real upperQuantileValue = hist->Quantile(0, uFract);
    {
    const Real windowSlope = ( upperPeggedValue - lowerPeggedValue ) / ( upperQuantileValue - lowerQuantileValue );
    const Real windowIntercept = ( upperPeggedValue - ( upperQuantileValue * windowSlope ) );
    const typename ImageType::PixelType adjustedWindowMinForClipping = ( clipMin - windowIntercept ) / ( windowSlope );
    const typename ImageType::PixelType adjustedWindowMaxForClipping = ( clipMax - windowIntercept ) / ( windowSlope );

    std::cout << "Quantile: [" << lowerQuantileValue << "," << upperQuantileValue << "];"
              << "  Pegged: [" << lowerPeggedValue << "," << upperPeggedValue << "]." << std::endl;
    std::cout << "Mapping Window: [" << adjustedWindowMinForClipping << "," << adjustedWindowMaxForClipping << "]"
              << " to [" << clipMin << "," << clipMax << "]." << std::endl;
    intensityMapper->SetWindowMinimum(adjustedWindowMinForClipping);
    intensityMapper->SetWindowMaximum(adjustedWindowMaxForClipping);
    intensityMapper->SetOutputMinimum(clipMin);
    intensityMapper->SetOutputMaximum(clipMax);
    }
  intensityMapper->Update();
  typename ImageType::Pointer remappedImage = intensityMapper->GetOutput();
  remappedImage->DisconnectPipeline();
  return remappedImage;
}

#endif // __StandardizeMaskIntensity_h
