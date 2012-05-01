/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkInterpolateMeshFunction.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkInterpolateMeshFunction_h
#define __itkInterpolateMeshFunction_h

#include "itkMeshFunction.h"
#include "itkPointLocator2.h"

namespace itk
{
/** \class InterpolateMeshFunction
 * \brief Base class for all mesh interpolators.
 *
 * InterpolateMeshFunction is the base for all MeshFunctions that
 * interpolates mesh intensity at a given point position.
 * This class is templated over the input mesh type and the
 * coordinate representation type (e.g. float or double ).
 *
 * \warning This hierarchy of functions work only for meshes
 * with scalar pixel types. For meshes of vector pixel types
 * use VectorInterpolateMeshFunctions.
 *
 * \sa VectorInterpolateMeshFunction
 * \ingroup MeshFunctions MeshInterpolators
 *
 * */
template <class TInputMesh>
class ITK_EXPORT InterpolateMeshFunction :
  public         MeshFunction<TInputMesh,
                              ITK_TYPENAME NumericTraits<typename TInputMesh::PixelType>::RealType>
{
public:
  /** Standard class typedefs. */
  typedef InterpolateMeshFunction Self;
  typedef MeshFunction<TInputMesh,
                       typename NumericTraits<
                         typename TInputMesh::PixelType>::RealType>     Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(InterpolateMeshFunction, MeshFunction);

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputMeshType typedef support. */
  typedef typename Superclass::InputMeshType      InputMeshType;
  typedef typename InputMeshType::PointIdentifier PointIdentifier;
  typedef typename InputMeshType::CellIdentifier  CellIdentifier;

  /** Dimension underlying input mesh. */
  itkStaticConstMacro(MeshDimension, unsigned int,
                      Superclass::MeshDimension);

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** RealType typedef support. */
  typedef typename TInputMesh::PixelType                PixelType;
  typedef typename NumericTraits<PixelType>::RealType   RealType;
  typedef itk::CovariantVector<RealType, MeshDimension> DerivativeType;

  /**
   * Interpolate the mesh at a point position.
   * Returns the interpolated mesh intensity at a specified point position. The
   * mesh cell is located based on proximity to the point to be evaluated.
   *
   * FIXME: What to do if the point is far from the Mesh ?
   *
   */
  virtual OutputType Evaluate( const PointType& point ) const = 0;

  /** Evaluate the derivative of the scalar function at the
   *  specified point. */
  virtual void EvaluateDerivative( const PointType& point, DerivativeType & derivative ) const = 0;

  /** Prepare internal data structures of the PointLocator. This method must be
   * called before performing any call to Evaluate. */
  void Initialize();

protected:
  InterpolateMeshFunction();
  ~InterpolateMeshFunction();

  void PrintSelf(std::ostream& os, Indent indent) const;

  typedef PointLocator2<TInputMesh>          PointLocatorType;
  typedef typename PointLocatorType::Pointer PointLocatorPointer;

  typedef typename PointLocatorType::InstanceIdentifierVectorType InstanceIdentifierVectorType;

  /** Searches the k-nearest neighbors */
  void Search(const PointType & query, unsigned int numberOfNeighborsRequested,
              InstanceIdentifierVectorType& result) const;

  /** Searches the neighbors fallen into a hypersphere */
  void Search(const PointType & query, double radius, InstanceIdentifierVectorType& result) const;

  /** Return the value associated with the point identified by pointId. */
  void GetPointData( PointIdentifier pointId, PixelType * value ) const;

private:
  InterpolateMeshFunction( const Self & ); // purposely not implemented
  void operator=( const Self & );          // purposely not implemented

  PointLocatorPointer m_PointLocator;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInterpolateMeshFunction.txx"
#endif

#endif
