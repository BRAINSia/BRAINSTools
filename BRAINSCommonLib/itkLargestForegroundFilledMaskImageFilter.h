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
#ifndef __itkLargestForegroundFilledMaskImageFilter_h
#define __itkLargestForegroundFilledMaskImageFilter_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>

namespace itk
{
/**
 * \class LargestForegroundFilledMaskImageFilter
 * \author Hans J. Johnson
 *
 * This filter does a good job of finding a single largest connected
 * mask that separates the foreground object from the background.
 * It assumes that the corner voxels of the image belong to the backgound.
 * This filter was written for the purpose of finding the tissue
 * region of a brain image with no internal holes.
 *
 * The OtsuPercentile Thresholds are used to define the range of values
 * where the percentage of voxels falls beetween
 * (0+OtsuPercentileLowerThreshold) < "Intensities of Interest" <
 *(1-OtsuPercentileUpperThreshold).
 *
 * The ClosingSize specifies how many mm to dilate followed by
 * erode to fill holes that be present in the image.
 *
 * The DilateSize specifies how many mm to dilate
 * as a final step to include a small amount of surface background in addition
 *to the
 * tissue region present in the image.
 *
 * The image that is returned will be a binary image with foreground and
 *background
 * values specified by the user (defaults to 1 and 0 respectively).
 *
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class LargestForegroundFilledMaskImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Convenient type alias for simplifying declarations. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::ConstPointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputPixelType = typename InputImageType::PixelType;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using OutputPixelType = typename OutputImageType::PixelType;

  using Self = LargestForegroundFilledMaskImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using IntegerImageType = Image<unsigned short, OutputImageType::ImageDimension>;
  using IntegerPixelType = typename IntegerImageType::PixelType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LargestForegroundFilledMaskImageFilter, ImageToImageFilter);

  /** set Otsu Threshold */
  itkSetMacro(OtsuPercentileLowerThreshold, double);
  itkGetConstMacro(OtsuPercentileLowerThreshold, double);
  itkSetMacro(OtsuPercentileUpperThreshold, double);
  itkGetConstMacro(OtsuPercentileUpperThreshold, double);

  /** Short hand for setting both upper and lower
   * (0+OtsuPercentileThreshold) < "Intensities of Interest" <
   *(1-OtsuPercentileThreshold).
   */
  void
  SetOtsuPercentileThreshold(const double percentile)
  {
    this->SetOtsuPercentileLowerThreshold(percentile);
    this->SetOtsuPercentileUpperThreshold(1.0 - percentile);
  }

  double
  GetOtsuPercentileThreshold() const
  {
    return this->GetOtsuPercentileLowerThreshold();
  }

  /** The closing size in mm, this is rounded up to the next closest number of
   * voxel
   * by taking Spacing into account */
  itkSetMacro(ClosingSize, double);
  itkGetConstMacro(ClosingSize, double);
  /** The dilation size in mm, this is rounded up to the next closest number of
   * voxel
   * by taking Spacing into account */
  itkSetMacro(DilateSize, double);
  itkGetConstMacro(DilateSize, double);
  itkSetMacro(InsideValue, IntegerPixelType);
  itkGetMacro(InsideValue, IntegerPixelType);
  itkSetMacro(OutsideValue, IntegerPixelType);
  itkGetMacro(OutsideValue, IntegerPixelType);
  itkSetMacro(ThresholdCorrectionFactor, double);
  itkGetConstMacro(ThresholdCorrectionFactor, double);

protected:
  LargestForegroundFilledMaskImageFilter();
  ~LargestForegroundFilledMaskImageFilter() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  /** Returns true if more than two bins of informaiton are found,
   * returns false if only two bins of informaiton are found (i.e. found a
   *binary image).
   * Low and High are set to the ?????? */
  unsigned int
  SetLowHigh(InputPixelType & low, InputPixelType & high);

  void
  ImageMinMax(InputPixelType & imageMin, InputPixelType & imageMax) const;

  // No longer used  double m_OtsuPercentileThreshold;
  double           m_OtsuPercentileLowerThreshold{ 0.01 };
  double           m_OtsuPercentileUpperThreshold;
  double           m_ThresholdCorrectionFactor{ 1.0 };
  double           m_ClosingSize{ 9.0 };
  double           m_DilateSize{ 0.0 };
  IntegerPixelType m_InsideValue;
  IntegerPixelType m_OutsideValue;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkLargestForegroundFilledMaskImageFilter.hxx"
#endif

#endif
