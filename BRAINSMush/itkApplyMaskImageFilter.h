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
#ifndef __itkApplyMaskImageFilter_h
#define __itkApplyMaskImageFilter_h

#include <itkImageToImageFilter.h>

namespace itk
{
template <typename TInputImage, typename TOutputImage>
class ApplyMaskImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public: /* define methods available to everyone */
  /** Standard class type alias */
  using Self = ApplyMaskImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** method for creation through the object factory */
  itkNewMacro(Self);

  /** run-time type information (and related methods) */
  itkTypeMacro(ApplyMaskImageFilter, ImageToImageFilter);

  /** ImageDimension enumeration */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** type alias to describe the output image region type */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** inherited type alias */
  using InputImageType = typename Superclass::InputImageType;
  using InputImagePointer = typename Superclass::InputImagePointer;
  using InputImageConstPointer = typename Superclass::InputImageConstPointer;
  using OutputImageType = typename Superclass::OutputImageType;
  using OutputImagePointer = typename Superclass::OutputImagePointer;

  /** pixel related type alias */
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;

  itkSetMacro(InvertMask, bool);
  itkGetMacro(InvertMask, bool);

  /** Set/get the mask to be applied to the image */
  void
  SetMaskImage(const InputImageType * reference);

  const InputImageType *
  GetMaskImage(void);

protected: /* define methods available only to related classes */
  ApplyMaskImageFilter();
  ~ApplyMaskImageFilter() {}

  void
  PrintSelf(std::ostream & os, Indent indent) const;

  void
  GenerateData();

private: /* define methods available only to this class */
  bool m_InvertMask;

  ITK_DISALLOW_COPY_AND_ASSIGN(ApplyMaskImageFilter);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkApplyMaskImageFilter.hxx"
#endif

#endif
