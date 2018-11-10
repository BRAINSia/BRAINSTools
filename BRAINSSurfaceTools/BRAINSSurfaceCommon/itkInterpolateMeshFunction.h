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
template <typename TInputMesh>
class InterpolateMeshFunction :
  public         MeshFunction<TInputMesh,
                              typename NumericTraits<typename TInputMesh::PixelType>::RealType>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(InterpolateMeshFunction);

  /** Standard class type alias. */
  using Self = InterpolateMeshFunction;
  using Superclass = MeshFunction<TInputMesh,
                       typename NumericTraits<
                         typename TInputMesh::PixelType>::RealType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(InterpolateMeshFunction, MeshFunction);

  /** OutputType type alias support. */
  using OutputType = typename Superclass::OutputType;

  /** InputMeshType type alias support. */
  using InputMeshType = typename Superclass::InputMeshType;
  using PointIdentifier = typename InputMeshType::PointIdentifier;
  using CellIdentifier = typename InputMeshType::CellIdentifier;

  /** Dimension underlying input mesh. */
  static constexpr unsigned int MeshDimension = Superclass::MeshDimension;

  /** Point type alias support. */
  using PointType = typename Superclass::PointType;

  /** RealType type alias support. */
  using PixelType = typename TInputMesh::PixelType;
  using RealType = typename NumericTraits<PixelType>::RealType;
  using DerivativeType = itk::CovariantVector<RealType, MeshDimension>;

  /**
   * Interpolate the mesh at a point position.
   * Returns the interpolated mesh intensity at a specified point position. The
   * mesh cell is located based on proximity to the point to be evaluated.
   *
   * FIXME: What to do if the point is far from the Mesh ?
   *
   */
  virtual OutputType Evaluate( constexpr PointType& point ) const override   = 0;

  /** Evaluate the derivative of the scalar function at the
   *  specified point. */
  virtual void EvaluateDerivative( constexpr PointType& point, DerivativeType & derivative ) const = 0;

  /** Prepare internal data structures of the PointLocator. This method must be
   * called before performing any call to Evaluate. */
  void Initialize();

protected:
  InterpolateMeshFunction();
  ~InterpolateMeshFunction();

  void PrintSelf(std::ostream& os, Indent indent) const override;

  using PointLocatorType = PointLocator2<TInputMesh>;
  using PointLocatorPointer = typename PointLocatorType::Pointer;

  using InstanceIdentifierVectorType = typename PointLocatorType::InstanceIdentifierVectorType;

  /** Searches the k-nearest neighbors */
  void Search(const PointType & query, unsigned int numberOfNeighborsRequested,
              InstanceIdentifierVectorType& result) const;

  /** Searches the neighbors fallen into a hypersphere */
  void Search(const PointType & query, double radius, InstanceIdentifierVectorType& result) const;

  /** Return the value associated with the point identified by pointId. */
  void GetPointData( PointIdentifier pointId, PixelType * value ) const;

private:
  PointLocatorPointer m_PointLocator;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInterpolateMeshFunction.hxx"
#endif

#endif
