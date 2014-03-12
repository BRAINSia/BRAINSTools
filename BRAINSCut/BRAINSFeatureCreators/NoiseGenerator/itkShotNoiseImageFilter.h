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
#ifndef __itkShotNoiseImageFilter_h
#define __itkShotNoiseImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
/** \class ShotNoiseImageFilter
 *
 * \brief Alter an image with shot noise.
 *
 * The shot noise follows a Poisson distribution.
 *
 * \author Gaetan Lehmann
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa InPlaceImageFilter
 */
template <class TInputImage, class TOutputImage = TInputImage>
class ShotNoiseImageFilter :
  public
  InPlaceImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ShotNoiseImageFilter Self;
  typedef InPlaceImageFilter<
      TInputImage, TOutputImage>                             Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ShotNoiseImageFilter, InPlaceImageFilter);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageType       OutputImageType;
  typedef typename Superclass::OutputImagePointer    OutputImagePointer;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;

  /** Some convenient typedefs. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;

  itkGetConstMacro(Scale, double);
  itkSetMacro(Scale, double);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
                  (Concept::Convertible<typename TInputImage::PixelType,
                                        typename TOutputImage::PixelType> ) );
  // The following concept check doesn't seem to work with vector immages
  // itkConceptMacro(Input1Input2OutputDivisionOperatorsCheck,
  //                (Concept::DivisionOperators<typename TInputImage::PixelType,
  //                 double,
  //                 typename TOutputImage::PixelType>));
  /** End concept checking */
#endif
protected:
  ShotNoiseImageFilter();
  virtual ~ShotNoiseImageFilter()
  {
  };

  void PrintSelf(std::ostream & os, Indent indent) const;

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId );

private:
  ShotNoiseImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);       // purposely not implemented

  double m_Scale;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShotNoiseImageFilter.txx"
#endif

#endif
