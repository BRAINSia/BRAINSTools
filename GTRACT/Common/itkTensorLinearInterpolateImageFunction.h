/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkTensorLinearInterpolateImageFunction_h
#define _itkTensorLinearInterpolateImageFunction_h

#include "itkTensorInterpolateImageFunction.h"

namespace itk
{
template <class TInputImage, class TCoordRep = float>
class ITK_EXPORT TensorLinearInterpolateImageFunction :
  public         TensorInterpolateImageFunction<TInputImage, TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef TensorLinearInterpolateImageFunction                   Self;
  typedef TensorInterpolateImageFunction<TInputImage, TCoordRep> Superclass;
  typedef SmartPointer<Self>                                     Pointer;
  typedef SmartPointer<const Self>                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorLinearInterpolateImageFunction,
               TensorInterpolateImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::PixelType      PixelType;
  typedef typename Superclass::ValueType      ValueType;
  typedef typename Superclass::RealType       RealType;

  /** Grab the vector dimension from the superclass. */
  // itkStaticConstMacro(Dimension, unsigned int,
  //    Superclass::Dimension);

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is Vector<double,Dimension> */
  typedef typename Superclass::OutputType OutputType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(const ContinuousIndexType & index ) const;

protected:
  TensorLinearInterpolateImageFunction();
  ~TensorLinearInterpolateImageFunction()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  TensorLinearInterpolateImageFunction(const Self &); // purposely not
                                                      // implemented
  void operator=(const Self &);                       // purposely not

  // implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long m_Neighbors;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorLinearInterpolateImageFunction.txx"
#endif

#endif
