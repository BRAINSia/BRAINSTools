/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNaryRelabelImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/10/25 12:12:57 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
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
class ITK_EXPORT NaryRelabelImageFilter :
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
   * Defaults to NumericTraits<PixelType>::Zero.
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

  void GenerateData();

  void PrintSelf( std::ostream& os, Indent indent) const;

private:
  NaryRelabelImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);         // purposely not implemented

  InputImagePixelType m_BackgroundValue;
  bool                m_IgnoreCollision;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNaryRelabelImageFilter.hxx"
#endif

#endif
