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
#ifndef __itkBOBFFilter_h
#define __itkBOBFFilter_h

#include "itkImageToImageFilter.h"
#include "itkFixedArray.h"

namespace itk
{
/** \class BOBFilter
 */
template <typename TInputImage, typename TOutputImage>
class BOBFFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(BOBFFilter);

  /** Standard class type alias. */
  using Self = BOBFFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BOBFFilter, ImageToImageFilter);

  /** ImageDimension enumeration. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Typedef to describe the output image region type. */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Inherited type alias. */
  using InputImageType = typename Superclass::InputImageType;
  using InputImagePointer = typename Superclass::InputImagePointer;
  using InputImageConstPointer = typename Superclass::InputImageConstPointer;
  using OutputImageType = typename Superclass::OutputImageType;
  using OutputImagePointer = typename Superclass::OutputImagePointer;

  /** Pixel related type alias. */
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using IndexType = typename InputImageType::IndexType;
  using InputSizeType = typename InputImageType::SizeType;

  /** Set/Get the Input image. */
  void
  SetInputImage(const InputImageType * source)
  {
    this->SetInput(source);
  }

  const InputImageType *
  GetInputImage(void)
  {
    return this->GetInput();
  }

  /** Set the input mask */
  void
  SetInputMask(const InputImageType * image);

  /** Get the input mask */
  const InputImageType *
  GetInputMask(void);

  /** Set seed point. */

  /** Set/Get the lower threshold. The default is 0. */
  itkSetMacro(Lower, InputPixelType);
  itkGetMacro(Lower, InputPixelType);

  /** Set/Get the upper threshold. The default is the largest possible
   *  value for the InputPixelType. */
  itkSetMacro(Upper, InputPixelType);
  itkGetMacro(Upper, InputPixelType);

  /** Set/Get value to replace thresholded pixels. Pixels that lie *
   *  within Lower and Upper (inclusive) will be replaced with this
   *  value. The default is 1. */
  itkSetMacro(ReplaceValue, OutputPixelType);
  itkGetMacro(ReplaceValue, OutputPixelType);

  /** Set the radius of the neighborhood used for a mask. */
  itkSetMacro(Radius, InputSizeType);

  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Radius, InputSizeType);

  /** Set the Seed of the neighborhood used for a mask. */
  itkSetMacro(Seed, IndexType);

  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Seed, IndexType);

  /** Method to execute the Filter */
  void
  GenerateData() override;

protected:
  BOBFFilter();
  ~BOBFFilter() override {}

  /** Override VeriyInputInformation() since this filter does not expect
   * the input images to occupy the same physical space.
   *
   * \sa ProcessObject::VerifyInputInformation
   */
  void
  VerifyInputInformation() const override
  {}

private:
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  // std::vector<IndexType> m_Seeds;
  IndexType       m_Seed;
  InputPixelType  m_Lower;
  InputPixelType  m_Upper;
  OutputPixelType m_ReplaceValue;
  InputSizeType   m_Radius;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBOBFFilter.hxx"
#endif

#endif
