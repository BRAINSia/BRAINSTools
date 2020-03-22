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
#ifndef __itkMultiModeHistogramThresholdBinaryImageFilter_h
#define __itkMultiModeHistogramThresholdBinaryImageFilter_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>
#include <itkArray.h>

namespace itk
{
/**
 * \author Hans J. Johnson
 *
 * This filter
 *
 */
template <typename TInputImage, typename TOutputImage = Image<unsigned short, TInputImage::ImageDimension>>
class MultiModeHistogramThresholdBinaryImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Extract dimension from input and output image. */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Convenient type alias for simplifying declarations. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::ConstPointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputPixelType = typename InputImageType::PixelType;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using OutputPixelType = typename OutputImageType::PixelType;

  using Self = MultiModeHistogramThresholdBinaryImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using IntegerImageType = TOutputImage;
  using IntegerPixelType = typename IntegerImageType::PixelType;

  using ThresholdArrayType = Array<double>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiModeHistogramThresholdBinaryImageFilter, ImageToImageFilter);

  itkSetMacro(LinearQuantileThreshold, double);
  itkGetConstMacro(LinearQuantileThreshold, double);

  /** set Quantile Threshold Arrays */
  itkSetMacro(QuantileLowerThreshold, ThresholdArrayType);
  itkGetConstMacro(QuantileLowerThreshold, ThresholdArrayType);
  itkSetMacro(QuantileUpperThreshold, ThresholdArrayType);
  itkGetConstMacro(QuantileUpperThreshold, ThresholdArrayType);

  itkGetConstObjectMacro(BinaryPortionImage, IntegerImageType);
  itkSetObjectMacro(BinaryPortionImage, IntegerImageType);

  itkSetMacro(InsideValue, IntegerPixelType);
  itkGetConstMacro(InsideValue, IntegerPixelType);
  itkSetMacro(OutsideValue, IntegerPixelType);
  itkGetConstMacro(OutsideValue, IntegerPixelType);

protected:
  MultiModeHistogramThresholdBinaryImageFilter();
  ~MultiModeHistogramThresholdBinaryImageFilter() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  ThresholdArrayType m_QuantileLowerThreshold;
  ThresholdArrayType m_QuantileUpperThreshold;
  double             m_LinearQuantileThreshold{ 0.01 };

  typename IntegerImageType::Pointer m_BinaryPortionImage;

  IntegerPixelType m_InsideValue;
  IntegerPixelType m_OutsideValue;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMultiModeHistogramThresholdBinaryImageFilter.hxx"
#endif

#endif
