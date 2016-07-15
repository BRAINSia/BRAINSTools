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
#ifndef __itkNaryRelabelImageFilter_h
#define __itkNaryRelabelImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkImageIterator.h"
#include "itkArray.h"

namespace itk
{
/** \class NaryRelabelImageFilter
 *
 * \brief relabel an combine the labels from several inputs
 *
 * This filter search all the label in the image, and give them a new value so they are all
 * consecutive. It then do the same with the second image, and give them the labels immediately
 * after the ones of the first image. It then do the same with the third image, etc.
 * The new labels are then copied to the output image.
 * Contrary to NaryLabelImageFilter, the user can't easily identify a label in the input image.
 * However, he is sure to get only consecutive labels in the output image (if the is no labeled
 * region collision), and no label collision.
 * This filter take one or more images as input, and produce a
 * single output image.
 * The SetIgnoreCollision(bool) method let the user choose to ignore the
 * labeled region collision or not. By default, they are ignored.
 * The SetBackgroundValue(OutputPixelType) let the user set the
 * value of the background label.
 *
 */

template <class TInputImage, class TOutputImage>
class NaryRelabelImageFilter :
  public         InPlaceImageFilter<TInputImage, TOutputImage>

{
public:
  /** Standard class typedefs. */
  typedef NaryRelabelImageFilter                        Self;
  typedef InPlaceImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NaryRelabelImageFilter, InPlaceImageFilter);

  /** Some typedefs. */
  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename InputImageType::PixelType   InputImagePixelType;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;
  typedef std::vector<InputImagePixelType>     NaryArrayType;

  /** ImageDimension constants */
  itkStaticConstMacro(
    InputImageDimension, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(
    OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck,
                  (Concept::SameDimension<InputImageDimension, OutputImageDimension> ) );
  itkConceptMacro(OutputHasZeroCheck,
                  (Concept::HasZero<OutputImagePixelType> ) );
  /** End concept checking */
#endif

  /**
   * Set/Get the value used as "background" in the images.
   * Defaults to NumericTraits<PixelType>::ZeroValue().
   */
  itkSetMacro(BackgroundValue, InputImagePixelType);
  itkGetConstMacro(BackgroundValue, InputImagePixelType);

  itkSetMacro(IgnoreCollision, bool);
  itkGetConstMacro(IgnoreCollision, bool);
  itkBooleanMacro(IgnoreCollision);
protected:
  NaryRelabelImageFilter();
  virtual ~NaryRelabelImageFilter()
  {
  };

  void GenerateData() ITK_OVERRIDE;

  void PrintSelf( std::ostream& os, Indent indent) const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(NaryRelabelImageFilter);

  InputImagePixelType m_BackgroundValue;
  bool                m_IgnoreCollision;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNaryRelabelImageFilter.hxx"
#endif

#endif
