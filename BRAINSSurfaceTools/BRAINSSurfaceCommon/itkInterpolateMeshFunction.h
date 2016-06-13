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
class InterpolateMeshFunction :
  public         MeshFunction<TInputMesh,
                              typename NumericTraits<typename TInputMesh::PixelType>::RealType>
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
  virtual OutputType Evaluate( const PointType& point ) const ITK_OVERRIDE = 0;

  /** Evaluate the derivative of the scalar function at the
   *  specified point. */
  virtual void EvaluateDerivative( const PointType& point, DerivativeType & derivative ) const = 0;

  /** Prepare internal data structures of the PointLocator. This method must be
   * called before performing any call to Evaluate. */
  void Initialize();

protected:
  InterpolateMeshFunction();
  ~InterpolateMeshFunction();

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

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
  ITK_DISALLOW_COPY_AND_ASSIGN(InterpolateMeshFunction);

  PointLocatorPointer m_PointLocator;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInterpolateMeshFunction.hxx"
#endif

#endif
