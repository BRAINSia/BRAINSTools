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
#ifndef __itkBlendImageFilter_h
#define __itkBlendImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class BlendImageFilter
 *  \brief Blend 2 images based using weights for each images
 */
template <typename TInputImage, typename TOutputImage>
class BlendImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(BlendImageFilter);

  /** Standard class type alias. */
  using Self = BlendImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BlendImageFilter, ImageToImageFilter);

  /** Some convenient type alias. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputImagePixelType = typename InputImageType::PixelType;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using OutputImagePixelType = typename OutputImageType::PixelType;

  /** ImageDimension enumeration */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Input and output images must be the same dimension, or the output's
      dimension must be one less than that of the input. */
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(ImageDimensionCheck, (Concept::SameDimension<Self::InputImageDimension, Self::OutputImageDimension>));
  /** End concept checking */
#endif

  /** Set the blend amounts for each input image.
   * set before the update of the filter.
   */
  itkGetConstMacro(Blend1, double);
  itkSetMacro(Blend1, double);
  itkGetConstMacro(Blend2, double);
  itkSetMacro(Blend2, double);

  void
  SetInput1(const TInputImage * image1)
  {
    this->SetNthInput(0, const_cast<TInputImage *>(image1));
  }

  void
  SetInput2(const TInputImage * image2)
  {
    this->SetNthInput(1, const_cast<TInputImage *>(image2));
  }

protected:
  BlendImageFilter();
  ~BlendImageFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

private:
  double m_Blend1{ 1.0 }, m_Blend2{ 1.0 };
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlendImageFilter.hxx"
#endif

#endif
