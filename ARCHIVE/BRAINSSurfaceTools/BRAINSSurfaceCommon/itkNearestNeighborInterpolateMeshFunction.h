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
#ifndef __itkNearestNeighborInterpolateMeshFunction_h
#define __itkNearestNeighborInterpolateMeshFunction_h

#include "itkInterpolateMeshFunction.h"

namespace itk
{
/** \class NearestNeighborInterpolateMeshFunction
 * \brief Performs linear interpolation in the cell closest to the evaluated point.
 *
 * This class will first locate the cell that is closest to the evaluated
 * point, and then will compute on it the output value using linear
 * interpolation among the values at the points of the cell.
 *
 * \sa VectorNearestNeighborInterpolateMeshFunction
 * \ingroup MeshFunctions MeshInterpolators
 *
 * */
template < typename TInputMesh >
class NearestNeighborInterpolateMeshFunction : public InterpolateMeshFunction< TInputMesh >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( NearestNeighborInterpolateMeshFunction );

  /** Standard class type alias. */
  using Self = NearestNeighborInterpolateMeshFunction;
  using Superclass = InterpolateMeshFunction< TInputMesh >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( NearestNeighborInterpolateMeshFunction, InterpolateMeshFunction );

  /** OutputType type alias support. */
  using OutputType = typename Superclass::OutputType;

  /** InputMeshType type alias support. */
  using InputMeshType = typename Superclass::InputMeshType;

  /** Dimension underlying input mesh. */
  static constexpr unsigned int MeshDimension = Superclass::MeshDimension;

  /** Point type alias support. */
  using PointType = typename Superclass::PointType;
  using PointIdentifier = typename Superclass::PointIdentifier;

  /** RealType type alias support. */
  using PixelType = typename TInputMesh::PixelType;
  using RealType = typename Superclass::RealType;
  using DerivativeType = typename Superclass::DerivativeType;

  /**
   * Interpolate the mesh at a point position
   * Returns the interpolated mesh intensity at a specified point position. The
   * mesh cell is located based on proximity to the point to be evaluated.
   *
   * FIXME: What to do if the point is far from the Mesh ?
   *
   */
  virtual OutputType
  Evaluate( const PointType & point ) const override;

  virtual void
  EvaluateDerivative( const PointType & point, DerivativeType & derivative ) const override;

  using InstanceIdentifierVectorType = typename Superclass::InstanceIdentifierVectorType;

protected:
  NearestNeighborInterpolateMeshFunction();
  ~NearestNeighborInterpolateMeshFunction();

  void
  PrintSelf( std::ostream & os, Indent indent ) const override;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkNearestNeighborInterpolateMeshFunction.hxx"
#endif

#endif
