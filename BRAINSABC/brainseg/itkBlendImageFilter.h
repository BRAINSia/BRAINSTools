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
template <class TInputImage, class TOutputImage>
class ITK_EXPORT BlendImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef BlendImageFilter                              Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BlendImageFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename    InputImageType::Pointer    InputImagePointer;
  typedef typename    InputImageType::RegionType InputImageRegionType;
  typedef typename    InputImageType::PixelType  InputImagePixelType;

  typedef TOutputImage                             OutputImageType;
  typedef typename     OutputImageType::Pointer    OutputImagePointer;
  typedef typename     OutputImageType::RegionType OutputImageRegionType;
  typedef typename     OutputImageType::PixelType  OutputImagePixelType;

  /** ImageDimension enumeration */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Input and output images must be the same dimension, or the output's
      dimension must be one less than that of the input. */
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( ImageDimensionCheck,
                   ( Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),
                                            itkGetStaticConstMacro(OutputImageDimension)> ) );
  /** End concept checking */
#endif

  /** Set the blend amounts for each input image.
   * set before the update of the filter.
   */
  itkGetConstMacro(Blend1, double);
  itkSetMacro(Blend1, double);
  itkGetConstMacro(Blend2, double);
  itkSetMacro(Blend2, double);

  void SetInput1(const TInputImage *image1)
  {
    this->SetNthInput(0, const_cast<TInputImage *>(image1) );
  }

  void SetInput2(const TInputImage *image2)
  {
    this->SetNthInput(1, const_cast<TInputImage *>(image2) );
  }

protected:
  BlendImageFilter();
  virtual ~BlendImageFilter()
  {
  }

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

  virtual void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

private:
  BlendImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);   // purposely not implemented

  double m_Blend1, m_Blend2;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlendImageFilter.hxx"
#endif

#endif
