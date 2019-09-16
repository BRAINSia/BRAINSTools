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
#ifndef __StandardizeMaskIntensity_h
#define __StandardizeMaskIntensity_h
#include "itkMinimumMaximumImageCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkHistogram.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkNumericTraits.h"

#include "itkResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkConstantBoundaryCondition.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

template <typename ImageType>
typename ImageType::Pointer
ResampleImageWithIdentityTransform(const std::string &                              resamplerInterpolatorType,
                                   const typename ImageType::PixelType              defaultPixelValue,
                                   const typename ImageType::ConstPointer &         inputImage,
                                   const typename itk::ImageBase<3>::ConstPointer & referenceImage)
{
  using ResampleType = itk::ResampleImageFilter<ImageType, ImageType>;
  using ResamplePointer = typename ResampleType::Pointer;
  ResamplePointer resampler = ResampleType::New();
  resampler->SetInput(inputImage);
  // resampler->SetTransform(); // default transform is identity

  if (resamplerInterpolatorType == "BSpline")
  {
    using SplineInterpolatorType = typename itk::BSplineInterpolateImageFunction<ImageType, double, double>;

    // Spline interpolation, only available for input images, not
    // atlas
    typename SplineInterpolatorType::Pointer splineInt = SplineInterpolatorType::New();
    splineInt->SetSplineOrder(5);
    resampler->SetInterpolator(splineInt);
  }
  else if (resamplerInterpolatorType == "WindowedSinc")
  {
    using BoundaryConditionType = typename itk::ConstantBoundaryCondition<ImageType>;
    static constexpr unsigned int WindowedSincHammingWindowRadius = 5;
    using WindowFunctionType = itk::Function::HammingWindowFunction<WindowedSincHammingWindowRadius, double, double>;
    typedef typename itk::WindowedSincInterpolateImageFunction<ImageType,
                                                               WindowedSincHammingWindowRadius,
                                                               WindowFunctionType,
                                                               BoundaryConditionType,
                                                               double>
                                                   WindowedSincInterpolatorType;
    typename WindowedSincInterpolatorType::Pointer windowInt = WindowedSincInterpolatorType::New();
    resampler->SetInterpolator(windowInt);
  }
  else if (resamplerInterpolatorType == "NearestNeighbor")
  {
    using NearestNeighborInterpolatorType = typename itk::NearestNeighborInterpolateImageFunction<ImageType, double>;
    typename NearestNeighborInterpolatorType::Pointer nearestNeighborInt = NearestNeighborInterpolatorType::New();
    resampler->SetInterpolator(nearestNeighborInt);
  }
  else // Default to m_UseNonLinearInterpolation == "Linear"
  {
    using LinearInterpolatorType = typename itk::LinearInterpolateImageFunction<ImageType, double>;
    typename LinearInterpolatorType::Pointer linearInt = LinearInterpolatorType::New();
    resampler->SetInterpolator(linearInt);
  }

  resampler->SetDefaultPixelValue(defaultPixelValue);
  resampler->SetOutputParametersFromImage(referenceImage);
  resampler->Update();

  typename ImageType::Pointer resimg = resampler->GetOutput();
  return resimg;
}


constexpr signed short MAX_IMAGE_OUTPUT_VALUE = 4096;

/* * * * *
 *   StandardizeMaskIntensity trims the upper and lower fractions of the histogram inside the mask
 *   based on tail size fractions uFract and lFract, respectively.  Then, using these histogram-based
 *   upper and lower bounds, it rescales the signal interval they represent to an interval marked by
 *   the target upper and lower values, and clips the output image signal to an interval that encloses
 *   the target interval.
 *
 *   Uses Review Statistics Histogram's Quantile method.  IntensityWindowingImageFilter clips to the
 *   same bounds it scales to as lowerQuantileValue target, so lowerQuantileValue little univariate extrapolation was
 * done (see block comment).
 *
 * * * * */
template <typename ImageType, typename LabelImageType>
typename ImageType::Pointer
StandardizeMaskIntensity(typename ImageType::Pointer         image,
                         typename LabelImageType::Pointer    mask,
                         const double                        lFract,
                         const double                        uFract,
                         const typename ImageType::PixelType lowerPeggedValue,
                         const typename ImageType::PixelType upperPeggedValue,
                         const typename ImageType::PixelType clipMin,
                         const typename ImageType::PixelType clipMax)
{
  // the asserts, if you want them, go like this:
  // ImageType::PixelType is scalar.
  // 0.0 <= lFract < uFract <= 1.0
  // clipMin <= lowerPeggedValue < upperPeggedValue <= clipMax

  using MinimumMaximumImageCalculator = typename itk::MinimumMaximumImageCalculator<ImageType>;
  typename MinimumMaximumImageCalculator::Pointer wholeStatistics = MinimumMaximumImageCalculator::New();
  wholeStatistics->SetImage(image);
  wholeStatistics->Compute();
  typename ImageType::PixelType imgMin = wholeStatistics->GetMinimum();
  typename ImageType::PixelType imgMax = wholeStatistics->GetMaximum();

  int numBins = itk::Math::rnd(imgMax - imgMin + 1);
  if (numBins < 256)
  {
    numBins = 256;
  }
  if (numBins > MAX_IMAGE_OUTPUT_VALUE)
  {
    numBins = MAX_IMAGE_OUTPUT_VALUE;
  }

  constexpr typename LabelImageType::PixelType maskInteriorLabel = 1;
  typename LabelImageType::Pointer             internalMask;
  if (mask.IsNull())
  {
    internalMask = LabelImageType::New();
    internalMask->CopyInformation(image);
    internalMask->SetRegions(image->GetLargestPossibleRegion());
    internalMask->Allocate();
    internalMask->FillBuffer(maskInteriorLabel);
  }
  else
  {
    // resample input mask to the voxel lattice of the input image
    // mask and input image are already at the same physical space
    typename LabelImageType::Pointer resampledMask =
      ResampleImageWithIdentityTransform<LabelImageType>("NearestNeighbor", 0, mask.GetPointer(), image.GetPointer());

    typename itk::ThresholdImageFilter<LabelImageType>::Pointer thresholdFilter =
      itk::ThresholdImageFilter<LabelImageType>::New();

    thresholdFilter->SetInput(resampledMask);
    thresholdFilter->ThresholdAbove(1); // Values less than or equal to are set
                                        // to OutsideValue
    thresholdFilter->SetOutsideValue(maskInteriorLabel);
    thresholdFilter->Update();
    internalMask = thresholdFilter->GetOutput();
  }

  using LabelStatisticsImageFilter = typename itk::LabelStatisticsImageFilter<ImageType, LabelImageType>;
  typename LabelStatisticsImageFilter::Pointer maskedStatistics = LabelStatisticsImageFilter::New();
  maskedStatistics->SetInput(image); // i.clipMin., image.
  maskedStatistics->SetLabelInput(internalMask);
  maskedStatistics->UseHistogramsOn();
  maskedStatistics->SetHistogramParameters(numBins, imgMin, imgMax);
  maskedStatistics->Update();
  typename LabelStatisticsImageFilter::HistogramType::Pointer hist = maskedStatistics->GetHistogram(maskInteriorLabel);
  if (hist.IsNull())
  {
    itkGenericExceptionMacro("histogram had no value for label " << maskInteriorLabel);
  }

  // Remark:  Since itk's filter uses the same bounds to scale and to clip,
  // and since the clipping bounds are outside (surrounding) the scaling
  // bounds, we need to compute adjustedWindowMaxForClipping and
  // adjustedWindowMinForClipping
  // to fake out ITK developers and get both
  // jobs done at once.  Fortunately, we may use itkHistogram's Quantile
  // routine:
  using IntensityWindowingImageFilter = typename itk::IntensityWindowingImageFilter<ImageType>;
  typename IntensityWindowingImageFilter::Pointer intensityMapper = IntensityWindowingImageFilter::New();
  intensityMapper->SetInput(image); // i.clipMin.,
                                    // image.
  // NOTE:  The math below is to extend the range to the clipping region.
  //
  //
  using Real = typename itk::NumericTraits<typename ImageType::PixelType>::ScalarRealType;
  const Real lowerQuantileValue = hist->Quantile(0, lFract);
  const Real upperQuantileValue = hist->Quantile(0, uFract);
  {
    const Real windowSlope = (upperPeggedValue - lowerPeggedValue) / (upperQuantileValue - lowerQuantileValue);
    const Real windowIntercept = (upperPeggedValue - (upperQuantileValue * windowSlope));
    const typename ImageType::PixelType adjustedWindowMinForClipping = (clipMin - windowIntercept) / (windowSlope);
    const typename ImageType::PixelType adjustedWindowMaxForClipping = (clipMax - windowIntercept) / (windowSlope);

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
