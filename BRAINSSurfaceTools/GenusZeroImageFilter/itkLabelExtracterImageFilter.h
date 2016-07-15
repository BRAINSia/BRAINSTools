/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelExtracterImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/07/17 14:37:03 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabelExtracterImageFilter_h
#define __itkLabelExtracterImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkConceptChecking.h"
#include "itkSimpleDataObjectDecorator.h"

namespace itk
{
/** \class LabelExtracterImageFilter
 *
 * \brief Change Sets of Labels
 *
 * This filter produces an output image whose pixels
 * are either copied from the input if they are not being changed
 * or are rewritten based on the change parameters
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * The filter expect both images to have the same number of dimensions.
 *
 * This is a trivial change of the original filter by
 *
 * \author Tim Kelliher. GE Research, Niskayuna, NY.
 * \note This work was supported by a grant from DARPA, executed by the
 *  U.S. Army Medical Research and Materiel Command/TATRC Assistance
 *  Agreement, Contract#W81XWH-05-2-0059.
 *
 * The difference is that values not in the map are automatically set to zero
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */

#include <map>

namespace Functor
{
template <class TInput, class TOutput>
class LabelExtracter
{
public:
  LabelExtracter()
  {
  }

  ~LabelExtracter()
  {
  }

  typedef std::map<TInput, TOutput> ChangeMapType;

  bool operator!=( const LabelExtracter & other ) const
  {
    if( m_ChangeMap != other.m_ChangeMap )
      {
      return true;
      }
    return false;
  }

  bool operator==( const LabelExtracter & other ) const
  {
    return !( *this != other );
  }

  TInput GetChange( const TInput & original )
  {
    return m_ChangeMap[original];
  }

  void SetChange( const TInput & original, const TOutput & result )
  {
    m_ChangeMap[original] = result;
  }

  void SetChangeMap( const ChangeMapType & changeMap )
  {
    m_ChangeMap = changeMap;
  }

  void ClearChangeMap()
  {
    m_ChangeMap.clear();
  }

  inline TOutput operator()( const TInput & A )
  {
    if( m_ChangeMap.find(A) != m_ChangeMap.end() )
      {
      return m_ChangeMap[A];
      }
    // here is the change
    return 0;
    // return A;
  }

private:

  ChangeMapType m_ChangeMap;
};
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT LabelExtracterImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::LabelExtracter<
                            typename TInputImage::PixelType,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef LabelExtracterImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputImage, TOutputImage,
                                  Functor::LabelExtracter<
                                    typename TInputImage::PixelType,
                                    typename TOutputImage::PixelType>
                                  >  Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

/** Method for creation through the object factory. */
  itkNewMacro(Self);

/** Run-time type information (and related methods). */
  itkTypeMacro(LabelExtracterImageFilter, UnaryFunctorImageFilter);

/** Pixel types. */
  typedef typename TInputImage::PixelType  InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;

/** Type of the change map to use for change requests */
  typedef std::map<InputPixelType, OutputPixelType> ChangeMapType;

/** Set up a change of a single label */
  void SetChange( const InputPixelType & original, const OutputPixelType & result );

/** Set the entire change map */
  void SetChangeMap( const ChangeMapType & changeMap );

/** Clears the entire change map */
  void ClearChangeMap();

#ifdef ITK_USE_CONCEPT_CHECKING
/** Begin concept checking */
  itkConceptMacro( InputConvertibleToOutputCheck,
                   ( Concept::Convertible<InputPixelType, OutputPixelType> ) );
  itkConceptMacro( PixelTypeComparable,
                   ( Concept::Comparable<InputPixelType> ) );
/** End concept checking */
#endif
protected:
  LabelExtracterImageFilter();
  virtual ~LabelExtracterImageFilter()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(LabelExtracterImageFilter);

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelExtracterImageFilter.hxx"
#endif

#endif
