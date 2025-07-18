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
#ifndef __itkComputeHistogramQuantileThresholds_h
#define __itkComputeHistogramQuantileThresholds_h

#include <itkImage.h>
#include <itkNumericTraits.h>

namespace itk
{
/**
 * \class ComputeHistogramQuantileThresholds
 * \author Hans J. Johnson
 *
 * This filter just computes Histogram Quantile Thresholds.  It does not apply
 *the thresholds.
 *
 */
template <typename TInputImage, typename TMaskImage>
class ComputeHistogramQuantileThresholds : public Object
{
public:
  /** Extract dimension from input and output image. */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

  /** Convenient type alias for simplifying declarations. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::ConstPointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputPixelType = typename InputImageType::PixelType;

  using Self = ComputeHistogramQuantileThresholds;
  using Superclass = Object;
  using Pointer = SmartPointer<Self>;
  using MaskPixelType = typename TMaskImage::PixelType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ComputeHistogramQuantileThresholds);

  /** set Quantile Threshold */
  itkSetMacro(QuantileLowerThreshold, double);
  itkGetConstMacro(QuantileLowerThreshold, double);
  itkSetMacro(QuantileUpperThreshold, double);
  itkGetConstMacro(QuantileUpperThreshold, double);

  itkGetConstMacro(LowerIntensityThresholdValue, typename InputImageType::PixelType);
  itkGetConstMacro(UpperIntensityThresholdValue, typename InputImageType::PixelType);
  itkGetConstMacro(NumberOfValidHistogramsEntries, unsigned int);

  itkGetConstObjectMacro(Image, InputImageType);
  itkSetConstObjectMacro(Image, InputImageType);

  itkSetMacro(ImageMin, typename TInputImage::PixelType);
  itkGetConstMacro(ImageMin, typename TInputImage::PixelType);
  itkSetMacro(ImageMax, typename TInputImage::PixelType);
  itkGetConstMacro(ImageMax, typename TInputImage::PixelType);

  itkGetConstObjectMacro(BinaryPortionImage, TMaskImage);
  itkSetObjectMacro(BinaryPortionImage, TMaskImage);

  void
  Calculate();

protected:
  ComputeHistogramQuantileThresholds();
  ~ComputeHistogramQuantileThresholds() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  void
  ImageMinMax(InputPixelType & ImageMin, InputPixelType & ImageMax);

  InputImagePointer            m_Image;
  typename TMaskImage::Pointer m_BinaryPortionImage;

  double       m_QuantileLowerThreshold{ 0.0 };
  double       m_QuantileUpperThreshold{ 1.0 };
  unsigned int m_NumberOfValidHistogramsEntries{ 0 };

  typename TInputImage::PixelType m_ImageMin;
  typename TInputImage::PixelType m_ImageMax;

  typename InputImageType::PixelType m_LowerIntensityThresholdValue;
  typename InputImageType::PixelType m_UpperIntensityThresholdValue;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkComputeHistogramQuantileThresholds.hxx"
#endif

#endif
