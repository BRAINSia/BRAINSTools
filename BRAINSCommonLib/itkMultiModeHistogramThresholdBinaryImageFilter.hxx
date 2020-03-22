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
#include "itkMultiModeHistogramThresholdBinaryImageFilter.h"
#include "itkComputeHistogramQuantileThresholds.h"

#include <itkBinaryThresholdImageFilter.h>
#include <itkMultiplyImageFilter.h>

#include <itkNumericTraits.h>
#include <itkMinimumMaximumImageFilter.h>
#include <itkCastImageFilter.h>

namespace itk
{
template <typename TInputImage, typename TOutputImage>
MultiModeHistogramThresholdBinaryImageFilter<TInputImage,
                                             TOutputImage>::MultiModeHistogramThresholdBinaryImageFilter()
  : m_QuantileLowerThreshold(1)
  , // temporarily estimate how many SetInput images
    // there are
  m_QuantileUpperThreshold(1)
  , m_InsideValue(NumericTraits<typename IntegerImageType::PixelType>::OneValue())
  , m_OutsideValue(NumericTraits<typename IntegerImageType::PixelType>::ZeroValue())
{
  m_QuantileLowerThreshold.Fill(0.0);
  m_QuantileUpperThreshold.Fill(1.0);
}

template <typename TInputImage, typename TOutputImage>
MultiModeHistogramThresholdBinaryImageFilter<TInputImage,
                                             TOutputImage>::~MultiModeHistogramThresholdBinaryImageFilter() = default;

template <typename TInputImage, typename TOutputImage>
void
MultiModeHistogramThresholdBinaryImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os,
                                                                                   Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "QuantileLowerThreshold " << m_QuantileLowerThreshold << " "
     << "QuantileUpperThreshold " << m_QuantileUpperThreshold << " "
     << "InsideValue " << m_InsideValue << " "
     << "OutsideValue " << m_OutsideValue << std::endl;
}

template <typename TInputImage, typename TOutputImage>
void
MultiModeHistogramThresholdBinaryImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  this->AllocateOutputs();

  typename IntegerImageType::Pointer accumulate = IntegerImageType::New();

  const unsigned int NumInputs = this->GetNumberOfInputs();
  for (unsigned int j = 0; j < NumInputs; ++j)
  {
    // Compute the quantile regions for linearizing the percentages.
    using ImageCalcType = ComputeHistogramQuantileThresholds<TInputImage, TOutputImage>;
    typename ImageCalcType::Pointer ImageCalc = ImageCalcType::New();
    ImageCalc->SetImage(this->GetInput(j));

    ImageCalc->SetQuantileLowerThreshold(m_LinearQuantileThreshold);
    ImageCalc->SetQuantileUpperThreshold(1.0 - m_LinearQuantileThreshold);

    std::cout << "Quantile Thresholds: [ " << m_QuantileLowerThreshold.GetElement(j) << ", "
              << m_QuantileUpperThreshold.GetElement(j) << " ]" << std::endl;

    ImageCalc->SetBinaryPortionImage(this->m_BinaryPortionImage);
    ImageCalc->Calculate();

    const typename InputImageType::PixelType thresholdLowerLinearRegion = ImageCalc->GetLowerIntensityThresholdValue();
    const typename InputImageType::PixelType thresholdUpperLinearRegion = ImageCalc->GetUpperIntensityThresholdValue();
    const typename InputImageType::PixelType imageMinValue = ImageCalc->GetImageMin();
    const typename InputImageType::PixelType imageMaxValue = ImageCalc->GetImageMax();
    const unsigned int                       numNonZeroHistogramBins = ImageCalc->GetNumberOfValidHistogramsEntries();

    typename InputImageType::PixelType thresholdLowerLinearRegion_foreground;
    if (numNonZeroHistogramBins <= 2)
    {
      thresholdLowerLinearRegion_foreground = thresholdUpperLinearRegion;
    }
    else
    {
      thresholdLowerLinearRegion_foreground = thresholdLowerLinearRegion;
    }

    std::cout << "LowHigh Thresholds: [ " << thresholdLowerLinearRegion << ", " << thresholdLowerLinearRegion_foreground
              << ", " << thresholdUpperLinearRegion << " ]" << std::endl;

    using ThresholdFilterType = BinaryThresholdImageFilter<InputImageType, IntegerImageType>;
    typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
    threshold->SetInput(this->GetInput(j));
    threshold->SetInsideValue(this->m_InsideValue);
    threshold->SetOutsideValue(this->m_OutsideValue);
    typename InputImageType::PixelType intensity_thresholdLowerLinearRegion;
    typename InputImageType::PixelType intensity_thresholdUpperLinearRegion;
    if (m_QuantileLowerThreshold.GetElement(j) < m_LinearQuantileThreshold)
    {
      const double range = (m_LinearQuantileThreshold - 0.0);
      const double percentValue = (m_QuantileLowerThreshold.GetElement(j) - 0.0) / range;
      intensity_thresholdLowerLinearRegion = static_cast<typename InputImageType::PixelType>(
        imageMinValue + (thresholdLowerLinearRegion_foreground - imageMinValue) * percentValue);
    }
    else
    {
      const double range = (1.0 - m_LinearQuantileThreshold) - m_LinearQuantileThreshold;
      const double percentValue = (m_QuantileLowerThreshold.GetElement(j) - m_LinearQuantileThreshold) / range;
      intensity_thresholdLowerLinearRegion = static_cast<typename InputImageType::PixelType>(
        thresholdLowerLinearRegion_foreground +
        (thresholdUpperLinearRegion - thresholdLowerLinearRegion_foreground) * percentValue);
    }
    if (m_QuantileUpperThreshold.GetElement(j) > (1.0 - m_LinearQuantileThreshold))
    {
      const double range = 1.0 - m_LinearQuantileThreshold;
      const double percentValue = (m_QuantileUpperThreshold.GetElement(j) - m_LinearQuantileThreshold) / range;
      intensity_thresholdUpperLinearRegion = static_cast<typename InputImageType::PixelType>(
        thresholdUpperLinearRegion + (imageMaxValue - thresholdUpperLinearRegion) * percentValue);
    }
    else
    {
      const double range = (1.0 - m_LinearQuantileThreshold) - m_LinearQuantileThreshold;
      const double percentValue = (m_QuantileUpperThreshold.GetElement(j) - m_LinearQuantileThreshold) / range;
      intensity_thresholdUpperLinearRegion = static_cast<typename InputImageType::PixelType>(
        thresholdLowerLinearRegion_foreground +
        (thresholdUpperLinearRegion - thresholdLowerLinearRegion_foreground) * percentValue);
    }
    std::cout << "DEBUG:MINMAX:DEBUG: [" << imageMinValue << "," << imageMaxValue << "]" << std::endl;
    std::cout << "DEBUG:LINLOWHIGH:DEBUG: [" << thresholdLowerLinearRegion << "," << thresholdUpperLinearRegion << "]"
              << std::endl;
    std::cout << "DEBUG:RANGE:DEBUG:  [" << intensity_thresholdLowerLinearRegion << ","
              << intensity_thresholdUpperLinearRegion << "]" << std::endl;
    threshold->SetLowerThreshold(intensity_thresholdLowerLinearRegion);
    threshold->SetUpperThreshold(intensity_thresholdUpperLinearRegion);
    // threshold->SetUpperThreshold( NumericTraits<typename
    // InputImageType::PixelType>::max() );
    threshold->Update();
    typename IntegerImageType::Pointer thresholdImage = threshold->GetOutput();

    if (j == 0)
    {
      accumulate = thresholdImage;
    }
    else
    {
      using IntersectMasksFilterType = MultiplyImageFilter<IntegerImageType, IntegerImageType>;
      if (accumulate->GetLargestPossibleRegion().GetSize() != thresholdImage->GetLargestPossibleRegion().GetSize())
      {
        itkExceptionMacro(<< "Image data size mismatch " << accumulate->GetLargestPossibleRegion().GetSize()
                          << " != " << thresholdImage->GetLargestPossibleRegion().GetSize() << "." << std::endl);
      }
      if (accumulate->GetSpacing() != thresholdImage->GetSpacing())
      {
        itkExceptionMacro(<< "Image data spacing mismatch " << accumulate->GetSpacing()
                          << " != " << thresholdImage->GetSpacing() << "." << std::endl);
      }
      if (accumulate->GetDirection() != thresholdImage->GetDirection())
      {
        itkExceptionMacro(<< "Image data spacing mismatch " << accumulate->GetDirection()
                          << " != " << thresholdImage->GetDirection() << "." << std::endl);
      }
      if (accumulate->GetOrigin() != thresholdImage->GetOrigin())
      {
        itkExceptionMacro(<< "Image data spacing mismatch " << accumulate->GetOrigin()
                          << " != " << thresholdImage->GetOrigin() << "." << std::endl);
      }
      typename IntersectMasksFilterType::Pointer intersect = IntersectMasksFilterType::New();
      intersect->SetInput1(accumulate);
      intersect->SetInput2(thresholdImage);
      intersect->Update();
      accumulate = intersect->GetOutput();
    }
  }

  using outputCasterType = CastImageFilter<IntegerImageType, OutputImageType>;
  typename outputCasterType::Pointer outputCaster = outputCasterType::New();
  outputCaster->SetInput(accumulate);

  outputCaster->GraftOutput(this->GetOutput());
  outputCaster->Update();
  this->GraftOutput(outputCaster->GetOutput());
  //  typename OutputImageType::Pointer outputMaskImage =
  // outputCaster->GetOutput();
  //  return outputMaskImage;
}
} // namespace itk
