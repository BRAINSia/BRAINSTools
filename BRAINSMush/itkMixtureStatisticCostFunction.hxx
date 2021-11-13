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
#ifndef itkMixtureStatisticsCostFunction_hxx_
#define itkMixtureStatisticsCostFunction_hxx_

#include <cassert>

namespace itk
{
template <typename TFirstImage, typename TSecondImage>
MixtureStatisticCostFunction<TFirstImage, TSecondImage>::MixtureStatisticCostFunction()
{}

template <typename TFirstImage, typename TSecondImage>
MixtureStatisticCostFunction<TFirstImage, TSecondImage>::~MixtureStatisticCostFunction()
{
}

template <typename TFirstImage, typename TSecondImage>
typename MixtureStatisticCostFunction<TFirstImage, TSecondImage>::MeasureType
MixtureStatisticCostFunction<TFirstImage, TSecondImage>::GetValue(const ParametersType & parameters) const
{
  const double & secondImageWeighting = parameters[0];
  // weighted mush_sum of additive statistics

  // measure each image and each squared image within the mask

  using FirstConstIteratorType = typename itk::ImageRegionConstIterator<typename Self::FirstImageType>;
  FirstConstIteratorType firstIt(m_FirstImage, m_FirstImage->GetLargestPossibleRegion());

  using SecondConstIteratorType = typename itk::ImageRegionConstIterator<typename Self::SecondImageType>;
  SecondConstIteratorType secondIt(m_SecondImage, m_SecondImage->GetLargestPossibleRegion());

  using MaskConstIteratorType = typename itk::ImageRegionConstIterator<typename Self::ImageMaskType>;
  MaskConstIteratorType maskIt(m_ImageMask, m_ImageMask->GetLargestPossibleRegion());

  if (m_FirstImage->GetLargestPossibleRegion().GetSize() != m_SecondImage->GetLargestPossibleRegion().GetSize() ||
      m_FirstImage->GetLargestPossibleRegion().GetSize() != m_ImageMask->GetLargestPossibleRegion().GetSize())
  {
    itkGenericExceptionMacro(<< "ERROR: Image or Mask size do not match");
  }

  double l_NumberOfMaskVoxels = 0.0;
  double mush_sum = 0.0;
  double mush_sumsq = 0.0;

  for (maskIt.GoToBegin(), firstIt.GoToBegin(), secondIt.GoToBegin(); !maskIt.IsAtEnd();
       ++maskIt, ++firstIt, ++secondIt)
  {
    if (maskIt.Get() != 0)
    {
      const typename TFirstImage::PixelType &  firstValue = firstIt.Get();
      const typename TSecondImage::PixelType & secondValue = secondIt.Get();

      l_NumberOfMaskVoxels += 1.0;
      const auto mush_value = std::fabs(firstValue - secondImageWeighting * secondValue);
      mush_sum += mush_value;
      mush_sumsq += mush_value * mush_value;
    }
  }

  // convert 3 statistics into mixture mean and variance
  const double mean = mush_sum / l_NumberOfMaskVoxels;
  const double variance = (mush_sumsq / l_NumberOfMaskVoxels) - mean * mean;

  // the measure is the squared distance to each of the goals.
  const double error_mean = mean;
  const double error_variance = variance;
  m_Measure = error_variance;

  std::cout << mean << " " << error_mean << " ** " << m_Measure << "=======" << parameters << std::endl;
  return m_Measure;
}

template <typename TFirstImage, typename TSecondImage>
unsigned int
MixtureStatisticCostFunction<TFirstImage, TSecondImage>::GetNumberOfParameters() const
{
  // Return the number of parameters.
  return number_of_unknowns;
}


template <typename TFirstImage, typename TSecondImage>
void
MixtureStatisticCostFunction<TFirstImage, TSecondImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk
#endif // itkMixtureStatisticsCostFunction_hxx_
