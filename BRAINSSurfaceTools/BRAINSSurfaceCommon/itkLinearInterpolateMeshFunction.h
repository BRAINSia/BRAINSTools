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
#ifndef __itkLinearInterpolateMeshFunction_h
#define __itkLinearInterpolateMeshFunction_h

#include "itkInterpolateMeshFunction.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleBasisSystemCalculator.h"

namespace itk
{
/** \class LinearInterpolateMeshFunction
 * \brief Performs linear interpolation in the cell closest to the evaluated point.
 *
 * This class will first locate the cell that is closest to the evaluated
 * point, and then will compute on it the output value using linear
 * interpolation among the values at the points of the cell.
 *
 * \sa VectorLinearInterpolateMeshFunction
 * \ingroup MeshFunctions MeshInterpolators
 *
 * */
template <class TInputMesh>
class LinearInterpolateMeshFunction :
  public         InterpolateMeshFunction<TInputMesh>
{
public:
  /** Standard class typedefs. */
  typedef LinearInterpolateMeshFunction       Self;
  typedef InterpolateMeshFunction<TInputMesh> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LinearInterpolateMeshFunction, InterpolateMeshFunction);

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
  typedef typename TInputMesh::PixelType      PixelType;
  typedef typename Superclass::RealType       RealType;
  typedef typename Superclass::DerivativeType DerivativeType;
  typedef typename PointType::VectorType      VectorType;

  /**
   * Interpolate the mesh at a point position.
   * Returns the interpolated mesh intensity at a specified point position. The
   * mesh cell is located based on proximity to the point to be evaluated.
   *
   * FIXME: What to do if the point is far from the Mesh ?
   *
   */
  virtual OutputType Evaluate( const PointType& point ) const ITK_OVERRIDE;

  virtual void EvaluateDerivative( const PointType& point, DerivativeType & derivative ) const ITK_OVERRIDE;

  static void GetDerivativeFromPixelsAndBasis(PixelType pixelValue1, PixelType pixelValue2, PixelType pixelValue3,
                                              const VectorType & u12, const VectorType & u32,
                                              DerivativeType & derivative);

  template <class TArray, class TMatrix>
  static void GetJacobianFromVectorAndBasis(const TArray & pixelArray1, const TArray & pixelArray2,
                                            const TArray & pixelArray3, const VectorType & u12, const VectorType & u32,
                                            TMatrix & jacobian);

  typedef typename Superclass::InstanceIdentifierVectorType InstanceIdentifierVectorType;

  /** Find the triangle that contains the input point. Return the point Ids of the triangle vertices. */
  virtual bool FindTriangle( const PointType& point, InstanceIdentifierVectorType & pointIds ) const;

  /** Find the first triangle touching the closest point to the input point,
   * and return the point Ids of the triangle vertices. Note that this triangle
   * will not necessarily contain the input point. Thie method is a good fallback
   * position to use when the FindTriangle() fails to find a triangle that
   * contains the input point. */
  virtual void FindTriangleOfClosestPoint( const PointType& point, InstanceIdentifierVectorType & pointIds ) const;

  /** Set Sphere Center.  The implementation of this interpolator assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the coordinates of the center of the sphere
   * represented by the Mesh.
   */
  itkSetMacro( SphereCenter, PointType );
  itkGetConstMacro( SphereCenter, PointType );

  /** When this boolean flag is ON, and the Evaluate() method fails to locate a
   * triangle that contains the queried point position, the interpolator will
   * fall-back to using nearest neighbor interpolation, and will make the
   * Evaluate() method to simply return the value of the mesh point that is
   * closest to the queried point. This option is useful when dealing with
   * meshses whose geometry may have folding or poorly shaped triangles. */
  itkSetMacro( UseNearestNeighborInterpolationAsBackup, bool );
  itkGetConstMacro( UseNearestNeighborInterpolationAsBackup, bool );
  itkBooleanMacro( UseNearestNeighborInterpolationAsBackup );
protected:
  LinearInterpolateMeshFunction();
  ~LinearInterpolateMeshFunction();

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  virtual bool ComputeWeights( const PointType & point, const InstanceIdentifierVectorType & pointIds ) const;

  const RealType & GetInterpolationWeight( unsigned int ) const;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(LinearInterpolateMeshFunction);

  mutable VectorType m_V12;
  mutable VectorType m_V32;

  mutable VectorType m_U12;
  mutable VectorType m_U32;

  mutable RealType m_InterpolationWeights[MeshDimension];

  itkStaticConstMacro( SurfaceDimension, unsigned int, 2 );

  typedef TriangleBasisSystem<VectorType, SurfaceDimension>                  TriangleBasisSystemType;
  typedef TriangleBasisSystemCalculator<TInputMesh, TriangleBasisSystemType> TriangleBasisSystemCalculatorType;

  typename TriangleBasisSystemCalculatorType::Pointer m_TriangleBasisSystemCalculator;

  bool      m_UseNearestNeighborInterpolationAsBackup;
  PointType m_SphereCenter;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLinearInterpolateMeshFunction.hxx"
#endif

#endif
