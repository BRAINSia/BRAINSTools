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
#ifndef __itkMeanSquaresMeshToMeshMetric_h
#define __itkMeanSquaresMeshToMeshMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"
#include "itkMeshToMeshMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"

namespace itk
{
/** \class MeanSquaresMeshToMeshMetric
 * \brief Computes similarity between two meshes to be registered
 *
 * This Class is templated over the type of the fixed and moving
 * meshes to be compared.
 *
 * This metric computes the sum of squared differences between point values
 * in the moving mesh and point values in the fixed mesh. The spatial
 * correspondance between both images is established through a Transform.
 * Point values are taken from the Moving mesh. Their positions are mapped to
 * the Fixed mesh and result in general in non-vertex position on it. Values
 * at these non-vertex positions of the Fixed mesh are interpolated using a
 * user-selected Interpolator.
 *
 * \ingroup RegistrationMetrics
 */
template < typename TFixedMesh, typename TMovingMesh >
class MeanSquaresMeshToMeshMetric : public MeshToMeshMetric< TFixedMesh, TMovingMesh >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( MeanSquaresMeshToMeshMetric );

  /** Standard class type alias. */
  using Self = MeanSquaresMeshToMeshMetric;
  using Superclass = MeshToMeshMetric< TFixedMesh, TMovingMesh >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MeanSquaresMeshToMeshMetric, MeshToMeshMetric );

  /** Types transferred from the base class */
  using TransformType = typename Superclass::TransformType;
  using TransformPointer = typename Superclass::TransformPointer;
  using TransformParametersType = typename Superclass::TransformParametersType;
  using TransformJacobianType = typename Superclass::TransformJacobianType;
  using InputPointType = typename Superclass::InputPointType;
  using OutputPointType = typename Superclass::OutputPointType;

  using MeasureType = typename Superclass::MeasureType;
  using DerivativeType = typename Superclass::DerivativeType;
  using FixedMeshType = typename Superclass::FixedMeshType;
  using MovingMeshType = typename Superclass::MovingMeshType;
  using FixedMeshConstPointer = typename Superclass::FixedMeshConstPointer;
  using MovingMeshConstPointer = typename Superclass::MovingMeshConstPointer;
  using PointIterator = typename Superclass::PointIterator;
  using PointDataIterator = typename Superclass::PointDataIterator;

  using InterpolatorType = typename Superclass::InterpolatorType;

  using RealDataType = typename Superclass::RealDataType;
  using DerivativeDataType = typename Superclass::DerivativeDataType;

  /** Constants for the pointset dimensions */
  static constexpr unsigned int MovingMeshDimension = Superclass::MovingMeshDimension;
  static constexpr unsigned int FixedMeshDimension = Superclass::FixedMeshDimension;

  /** Get the derivatives of the match measure. */
  void
  GetDerivative( const TransformParametersType & parameters, DerivativeType & derivative ) const override;

  /**  Get the value for single valued optimizers. */
  MeasureType
  GetValue( const TransformParametersType & parameters ) const override;

  /**  Get value and derivatives for multiple valued optimizers. */
  void
  GetValueAndDerivative( const TransformParametersType & parameters, MeasureType & Value,
                         DerivativeType & Derivative ) const override;

protected:
  MeanSquaresMeshToMeshMetric();
  virtual ~MeanSquaresMeshToMeshMetric(){};

private:
  mutable unsigned int m_NumberOfPixelsCounted;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMeanSquaresMeshToMeshMetric.hxx"
#endif

#endif
