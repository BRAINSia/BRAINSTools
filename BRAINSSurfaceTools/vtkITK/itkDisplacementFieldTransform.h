/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDisplacementFieldTransform.h,v $
  Language:  C++
  Date:      $Date: 2008-07-27 16:28:11 $
  Version:   $Revision: 1.38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDisplacementFieldTransform_h
#define __itkDisplacementFieldTransform_h

#include <iostream>
#include "itkTransform.h"
#include "itkVectorImage.h"
#include "itkVector.h"
#include "itkImage.h"

namespace itk
{
template <
  typename TScalarType = double,             // Data type for scalars
  unsigned int NDimensions = 3>
// Number of dimensions
class DisplacementFieldTransform :
    public Transform<TScalarType, NDimensions, NDimensions>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(DisplacementFieldTransform);

  /** Standard class type alias. */
  using Self = DisplacementFieldTransform;
  using Superclass = Transform<TScalarType, NDimensions, NDimensions>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DisplacementFieldTransform, Transform );

  /** Dimension of the domain space. */
  static constexpr unsigned int SpaceDimension = NDimensions;

  /** Standard scalar type for this class. */
  using ScalarType = typename Superclass::ScalarType;

  /** Standard parameters container. */
  using ParametersType = typename Superclass::ParametersType;

  /** Standard Jacobian container. */
  using JacobianType = typename Superclass::JacobianType;

  /** Standard vector type for this class. */
  using InputVectorType = Vector<TScalarType,
                 Self::SpaceDimension>;
  using OutputVectorType = Vector<TScalarType,
                 Self::SpaceDimension>;

  /** Standard covariant vector type for this class. */
  using InputCovariantVectorType = CovariantVector<TScalarType,
                          Self::SpaceDimension>;
  using OutputCovariantVectorType = CovariantVector<TScalarType,
                          Self::SpaceDimension>;

  /** Standard vnl_vector type for this class. */
  using InputVnlVectorType = vnl_vector_fixed<TScalarType,
                           Self::SpaceDimension>;
  using OutputVnlVectorType = vnl_vector_fixed<TScalarType,
                           Self::SpaceDimension>;

  /** Standard coordinate point type for this class. */
  using InputPointType = Point<TScalarType,
                Self::SpaceDimension>;
  using OutputPointType = Point<TScalarType,
                Self::SpaceDimension>;

  /** This method sets the parameters of the transform.
   *
   * For efficiency, this transform does not make a copy of the parameters.
   * It only keeps a pointer to the input parameters. It assumes that the memory
   * is managed by the caller. Use SetParametersByValue to force the transform
   * to call copy the parameters.
   *
   * This method wraps each grid as itk::VectorImage's using the user-specified
   * grid region, spacing and origin.
   * NOTE: The grid region, spacing and origin must be set first by
   * SetFixedParameters().
   *
   */
  void SetParameters(const ParametersType & parameters);

  /** This method sets the fixed parameters of the transform.
   *
   * Fixed Parameters store the following information:
   *    Grid Size
   *    Grid Origin
   *    Grid Spacing
   *    Grid Direction -- optional
   *
   * The size of fixed parameters is either 3 * NDimensions or 3 *
   * NDimension + NDimension * NDimension.
   *
   * This function should be called when parameters are set by
   * SetParameters() or SetParametersByValue(). If the parameters are
   * set by SetImage(), then the size, origin, spacing, and direction
   * are extracted from the image.
   *
   */
  void SetFixedParameters(const ParametersType & parameters);

  /** This method sets the parameters of the transform.
   *
   * The parameters are (N * NDimension) number of N-D grid of
   * coefficients. Each N-D grid is represented as a flat array of
   * doubles (in the same configuration as an itk::VectorImage).  The
   * N arrays are then concatenated to form one parameter array.
   *
   * This methods makes a copy of the parameters while for
   * efficiency the SetParameters method does not.
   *
   * This method wraps each grid as itk::VectorImage's using the user
   * specified grid region, spacing and origin.
   * NOTE: The grid region, spacing and origin must be set first by
   * SetFixedParameters().
   *
   */
  void SetParametersByValue(const ParametersType & parameters);

  void SetIdentity();

  /** Get the Transformation Parameters. */
  virtual const ParametersType & GetParameters(void) const;

  /** Get the Transformation Fixed Parameters. */
  virtual const ParametersType & GetFixedParameters(void) const;

  /** Parameters as SpaceDimension number of images. */
  using InternalPixelType = typename ParametersType::ValueType;
  using ImageType = VectorImage<InternalPixelType, Self::SpaceDimension>;
  using ImagePointer = typename ImageType::Pointer;
  using PixelType = typename ImageType::PixelType;

  /** Get the array of coefficient images. */
  virtual ImagePointer GetImage()
  { return m_CoefficientImage; }
  virtual const ImagePointer GetImage() const
  { return m_CoefficientImage; }

  /** Set the array of coefficient images.
   *
   * This is for setting the Displacement Fields as an array of
   * SpaceDimension images. The grid region spacing and origin is
   * taken from the image. It is assume that the buffered region of
   * all the subsequent images are the same as the first image. Note
   * that no error checking is done.
   *
   * SetImage() also copies the parameters to m_InternalParametersBuffer.
   */
  virtual void SetImage( ImagePointer images );

  OutputPointType  TransformPoint(const InputPointType  & point ) const;

  using ContinuousIndexType = typename ContinuousIndex<ScalarType,
                                   Self::SpaceDimension>;

  /** Method to transform a vector -
   *  not applicable for this type of transform. */
  virtual OutputVectorType TransformVector(const InputVectorType &) const
  {
    itkExceptionMacro(<< "Method not applicable for deformable transform." );
    return OutputVectorType();
  }

  /** Method to transform a vnl_vector -
   *  not applicable for this type of transform */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const
  {
    itkExceptionMacro(<< "Method not applicable for deformable transform. ");
    return OutputVnlVectorType();
  }

  /** Method to transform a CovariantVector -
   *  not applicable for this type of transform */
  virtual OutputCovariantVectorType TransformCovariantVector(
    const InputCovariantVectorType &) const
  {
    itkExceptionMacro(<< "Method not applicable for deformable transfrom. ");
    return OutputCovariantVectorType();
  }

  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual const JacobianType & GetJacobian(const InputPointType  & point );

  /** Return the number of parameters that completely define the Transfom */
  virtual unsigned int GetNumberOfParameters(void) const;

  /** Indicates that this transform is linear. That is, given two
   * points P and Q, and scalar coefficients a and b, then
   *
   *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
   */
  virtual bool IsLinear() const { return false; }

  virtual bool GetEncodeParametersOnWrite() const
  {
    return true;
  }

protected:
  /** Print contents of an DisplacementFieldTransform. */
  void PrintSelf(std::ostream & os, Indent indent) const;

  DisplacementFieldTransform();
  virtual ~DisplacementFieldTransform();

  /** Wrap flat array into images of coefficients. */
  void WrapAsImages();

  /** Convert an input point to a continuous index inside the grid */
  void TransformPointToContinuousIndex(
    const InputPointType & point, ContinuousIndexType & index ) const;

private:

  /** Typedefs for specifying the extend to the grid. */
  using RegionType = ImageRegion<Self::SpaceDimension>;

  using IndexType = typename RegionType::IndexType;
  using SizeType = typename RegionType::SizeType;
  using SpacingType = typename ImageType::SpacingType;
  using DirectionType = typename ImageType::DirectionType;
  using OriginType = typename ImageType::PointType;

  DirectionType m_PointToIndex;
  DirectionType m_IndexToPoint;
  RegionType    m_LastSupportRegion;
  /** Array of images representing the displacement fields */
  ImagePointer m_CoefficientImage;

  /** Jacobian as SpaceDimension number of images. */
  using JacobianInternalPixelType = typename JacobianType::ValueType;
  using JacobianImageType = VectorImage<JacobianInternalPixelType,
                      Self::SpaceDimension>;
  using JacobianPixelType = typename JacobianImageType::PixelType;

  typename JacobianImageType::Pointer m_JacobianImage[NDimensions];

  /** Keep a pointer to the input parameters. */
  const ParametersType *m_InputParametersPointer;

  /** Internal parameters buffer. */
  ParametersType m_InternalParametersBuffer;

  OutputVectorType LinearInterpolateAtIndex( const ContinuousIndexType & index ) const;

  /** Check if a continuous index is inside the valid region. */
  bool InsideValidRegion( const ContinuousIndexType & index ) const;
}; // class DisplacementFieldTransform
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_DisplacementFieldTransform(_, EXPORT, x, y)        \
  namespace itk {                                                       \
  _( 3 ( class EXPORT DisplacementFieldTransform<ITK_TEMPLATE_3 x> ) )  \
  namespace Templates { typedef DisplacementFieldTransform<ITK_TEMPLATE_3 x> \
    DisplacementFieldTransform ## y; }                                  \
  }
#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDisplacementFieldTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkDisplacementFieldTransform.hxx"
#endif
#endif /* __itkDisplacementFieldTransform_h */
