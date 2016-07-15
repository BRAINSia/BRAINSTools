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
/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: itkOppositeImageFilter.h,v $
 *  Language:  C++
 *  Date:      $Date: 2009-02-19 21:18:10 $
 *  Version:   $Revision: 0.0 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/

#ifndef __itkOppositeImageFilter_h
#define __itkOppositeImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
/** \class OppositeImageFilter
  *
  * \brief Take the opposite of the input pixels.
  *
  * This filter is templated over the input image type
  * and the output image type.
  *
  * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
  *
  * \ingroup IntensityImageFilters  Multithreaded
  * \sa UnaryFunctorImageFilter
  */
namespace Functor
{
template <class TInput, class TOutput>
class Opposite
{
public:
  Opposite()
  {
  }

  ~Opposite()
  {
  }

  bool operator!=(const Opposite & other) const
  {
    return false;
  }

  bool operator==(const Opposite & other) const
  {
    return true;
  }

  inline TOutput operator()(const TInput & A) const
  {
    // We don't check if the TOutput can be signed.
    // It's up to the user to decide whether this makes sense.
    return static_cast<TOutput>( -A );
  }
};
}

template <class TInputImage, class TOutputImage>
class OppositeImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::Opposite<
                            typename TInputImage::PixelType,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef OppositeImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputImage, TOutputImage,
                                  Functor::Opposite<typename TInputImage::PixelType,
                                                    typename TOutputImage::PixelType> >
    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OppositeImageFilter, UnaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputConvertibleToOutputCheck,
                   ( Concept::Convertible<typename TInputImage::PixelType,
                                          typename TOutputImage::PixelType> ) );
  /** End concept checking */
#endif
protected:
  OppositeImageFilter()
  {
  }

  virtual ~OppositeImageFilter()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(OppositeImageFilter);
};
} // end namespace itk

#endif
