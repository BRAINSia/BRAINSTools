/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLinearInterpolateDeformationFieldMeshFunction.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLinearInterpolateDeformationFieldMeshFunction_h
#define __itkLinearInterpolateDeformationFieldMeshFunction_h

#include "itkLinearInterpolateMeshFunction.h"

namespace itk
{
/** \class LinearInterpolateDeformationFieldMeshFunction
 * \brief Performs linear interpolation of the deformation field represented as
 * destination points.
 *
 * This class will first locate the cell that is closest to the evaluated
 * point, and then will compute on it the output value using linear
 * interpolation among the values at the deformation field of the cell
 * vertexes.
 *
 * \sa VectorLinearInterpolateDeformationFieldMeshFunction
 * \ingroup MeshFunctions MeshInterpolators
 *
 * */
template <class TInputMesh, class TDestinationPointsContainer = typename TInputMesh::PointsContainer>
class ITK_EXPORT LinearInterpolateDeformationFieldMeshFunction :
  public         LinearInterpolateMeshFunction<TInputMesh>
{
public:
  /** Standard class typedefs. */
  typedef LinearInterpolateDeformationFieldMeshFunction Self;
  typedef LinearInterpolateMeshFunction<TInputMesh>     Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LinearInterpolateDeformationFieldMeshFunction, LinearInterpolateMeshFunction);

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputMeshType typedef support. */
  typedef typename Superclass::InputMeshType InputMeshType;

  /** Dimension underlying input mesh. */
  itkStaticConstMacro(MeshDimension, unsigned int, Superclass::MeshDimension);

  /** Point typedef support. */
  typedef typename Superclass::PointType       PointType;
  typedef typename Superclass::PointIdentifier PointIdentifier;
  typedef typename Superclass::CellIdentifier  CellIdentifier;

  /** RealType typedef support. */
  typedef typename TInputMesh::PixelType PixelType;
  typedef typename Superclass::RealType  RealType;
  typedef typename PointType::VectorType VectorType;

  /** Type for the container of destination points of the deformation field. */
  typedef TDestinationPointsContainer DestinationPointsContainerType;

  /**
   * Interpolate the mesh at a point position.
   * Returns the interpolated mesh intensity at a specified point position. The
   * mesh cell is located based on proximity to the point to be evaluated.
   * Returns false when the point is far from the mesh or when not triangle
   * can be found to contain the point.
   *
   */
  virtual bool Evaluate( const DestinationPointsContainerType * field, const PointType & point,
                         PointType & outputPoint ) const;

  /** Provide empty implementation of virtual method from the base class.
      This method is not expected to be used by this current class. */
  virtual OutputType Evaluate( const PointType& point ) const;

protected:
  LinearInterpolateDeformationFieldMeshFunction();
  ~LinearInterpolateDeformationFieldMeshFunction();

  void PrintSelf(std::ostream& os, Indent indent) const;

  typedef typename Superclass::InstanceIdentifierVectorType InstanceIdentifierVectorType;
private:
  LinearInterpolateDeformationFieldMeshFunction( const Self & ); // purposely not implemented
  void operator=( const Self & );                                // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLinearInterpolateDeformationFieldMeshFunction.txx"
#endif

#endif
