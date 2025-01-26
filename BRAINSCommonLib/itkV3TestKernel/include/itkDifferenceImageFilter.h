/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef __itkDifferenceImageFilter_h
#define __itkDifferenceImageFilter_h

#include "itkConfigure.h" // Needed to determine value of ITKV3_COMPATIBILITY
#ifdef ITKV3_COMPATIBILITY

#  include "itkTestingComparisonImageFilter.h"

namespace itk
{
/** \class DifferenceImageFilter
 * This filter is an alias to the TestingComparisonImageFilter
 * and is only here for backwards compatibility.
 *
 * This class has no implementation, thus no .hxx file is needed.
 * \ingroup ITKTestKernel
 */
template <typename TInputImage, typename TOutputImage>
class DifferenceImageFilter : public Testing::ComparisonImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(DifferenceImageFilter);

  /** Standard class type alias. */
  using Self = DifferenceImageFilter;
  using Superclass = Testing::ComparisonImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(DifferenceImageFilter);

  /** Some convenient type alias. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using OutputPixelType = typename OutputImageType::PixelType;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using RealType = typename NumericTraits<OutputPixelType>::RealType;
  using AccumulateType = typename NumericTraits<RealType>::AccumulateType;

  DifferenceImageFilter() {}

  virtual ~DifferenceImageFilter() {}

protected:
};
} // end namespace itk
#else
#  error For ITKv4 compatibility, use itk::Testing::ComparisonImageFilter instead of itk::DifferenceImageFilter
#endif

#endif
