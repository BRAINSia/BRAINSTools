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
#ifndef __itkAverageImageFilter_h
#define __itkAverageImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{
/** \class AverageImageFilter
 *
 * \brief This filter performs pixelwise averaging among an arbitrary number
 * of input images.
 *
 * \ingroup ITKBRAINSFilterPack
 * \par INPUTS
 * Input volumes must all contain the same size RequestedRegions. All input
 * images must have the same pixel type. All pixel types are supported that
 * supply operator=(), operator+=(), operator*=(double), and a cast to the
 * chosen output type (if not identical to input pixel type).
 *
 * \par OUTPUTS
 * The averaging filter produces a single output volume. Each output pixel
 * is assigned the average of the values from the corresponding input image
 * pixels.
 *
 * \par LIMITATIONS
 * For integer output pixel types, the result of the averaging is converted
 * to that type by a cast operation. There is currently no rounding
 * implemented.
 *
 * \author Torsten Rohlfing, SRI International, Neuroscience Program
 *
 * Funding for the implementation of this class was provided by NIAAA under
 * Grant No. AA05965, "CNS Deficits: Interaction of Age and Alcoholism",
 * PI: A. Pfefferbaum; and Grant No. AA13521, "INIA: Imaging Core",
 * PI: A. Pfefferbaum.
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_EXPORT AverageImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(AverageImageFilter);

  /** Standard class type alias. */
  using Self = AverageImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(AverageImageFilter);

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  using OutputPixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** Image type alias support */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::ConstPointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass type alias. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;


protected:
  AverageImageFilter()
  {
    this->DynamicMultiThreadingOff(); // NEEDED FOR ITKv5 backwards compatibility
  }
  ~AverageImageFilter() override = default;

  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

  void
  PrintSelf(std::ostream &, Indent) const override;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkAverageImageFilter.hxx"
#endif

#endif
