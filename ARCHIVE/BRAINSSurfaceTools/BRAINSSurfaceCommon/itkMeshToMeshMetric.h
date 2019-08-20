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
#ifndef __itkMeshToMeshMetric_h
#define __itkMeshToMeshMetric_h

#include "itkTransform.h"
#include "itkSingleValuedCostFunction.h"
#include "itkInterpolateMeshFunction.h"
#include "itkSpatialObject.h"

namespace itk
{
/** \class MeshToMeshMetric
 * \brief Computes similarity between two meshes.
 *
 * This Class is templated over the type of the two meshes.  It expects a
 * Transform to be plugged in.  This particular class is the base class for a
 * hierarchy of mesh to mesh metrics.
 *
 * This class computes a value that measures the similarity between the fixed
 * mesh and the transformed moving mesh.
 *
 * \ingroup RegistrationMetrics
 *
 */

template < typename TFixedMesh, typename TMovingMesh >
class MeshToMeshMetric : public SingleValuedCostFunction
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( MeshToMeshMetric );

  /** Standard class type alias. */
  using Self = MeshToMeshMetric;
  using Superclass = SingleValuedCostFunction;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Type used for representing point components  */
  using CoordinateRepresentationType = typename TFixedMesh::CoordRepType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MeshToMeshMetric, SingleValuedCostFunction );

  /**  Type of the moving Mesh. */
  using MovingMeshType = TMovingMesh;
  using MovingMeshPixelType = typename TMovingMesh::PixelType;
  using MovingMeshConstPointer = typename MovingMeshType::ConstPointer;

  /**  Type of the fixed Mesh. */
  using FixedMeshType = TFixedMesh;
  using FixedMeshConstPointer = typename FixedMeshType::ConstPointer;

  /** Constants for the pointset dimensions */
  static constexpr unsigned int MovingMeshDimension = TMovingMesh::PointDimension;
  static constexpr unsigned int FixedMeshDimension = TFixedMesh::PointDimension;

  using PointIterator = typename FixedMeshType::PointsContainer::ConstIterator;
  using PointDataIterator = typename FixedMeshType::PointDataContainer::ConstIterator;

  /**  Type of the Transform Base class */
  using TransformComputationType = typename NumericTraits< CoordinateRepresentationType >::RealType;
  using TransformType = Transform< TransformComputationType, Self::MovingMeshDimension, Self::FixedMeshDimension >;

  using TransformPointer = typename TransformType::Pointer;
  using InputPointType = typename TransformType::InputPointType;
  using OutputPointType = typename TransformType::OutputPointType;
  using TransformParametersType = typename TransformType::ParametersType;
  using TransformJacobianType = typename TransformType::JacobianType;

  /**  Type of the Interpolator Base class */
  using InterpolatorType = InterpolateMeshFunction< MovingMeshType >;
  using InterpolatorPointer = typename InterpolatorType::Pointer;
  using RealDataType = typename InterpolatorType::RealType;
  using DerivativeDataType = typename InterpolatorType::DerivativeType;

  /**  Type of the measure. */
  using MeasureType = Superclass::MeasureType;

  /**  Type of the derivative. */
  using DerivativeType = Superclass::DerivativeType;

  /**  Type of the parameters. */
  using ParametersType = Superclass::ParametersType;

  /**  Type for the mask of the fixed image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  using FixedMaskType = SpatialObject< Self::FixedMeshDimension >;
  using FixedMaskPointer = typename FixedMaskType::ConstPointer;

  /**  Type for the mask of the moving image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  using MovingMaskType = SpatialObject< Self::MovingMeshDimension >;
  using MovingMaskPointer = typename MovingMaskType::ConstPointer;

  /** Connect the Fixed Pointset.  */
  itkSetConstObjectMacro( FixedMesh, FixedMeshType );

  /** Get the Fixed Pointset. */
  itkGetConstObjectMacro( FixedMesh, FixedMeshType );

  /** Connect the Moving Pointset.  */
  itkSetConstObjectMacro( MovingMesh, MovingMeshType );

  /** Get the Moving Pointset. */
  itkGetConstObjectMacro( MovingMesh, MovingMeshType );

  /** Connect the Transform. */
  itkSetObjectMacro( Transform, TransformType );

  /** Get a pointer to the Transform.  */
  itkGetConstObjectMacro( Transform, TransformType );

  /** Set the parameters defining the Transform. */
  void
  SetTransformParameters( const ParametersType & parameters ) const;

  /** Connect the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the Interpolator.  */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Return the number of parameters required by the Transform */
  unsigned int
  GetNumberOfParameters( void ) const override
  {
    return m_Transform->GetNumberOfParameters();
  }

  /** Connect the FixedMask .  */
  itkSetConstObjectMacro( FixedMask, FixedMaskType );

  /** Get the Fixed Pointset. */
  itkGetConstObjectMacro( FixedMask, FixedMaskType );

  /** Connect the Moving Pointset.  */
  itkSetConstObjectMacro( MovingMask, MovingMaskType );

  /** Get the Moving Pointset. */
  itkGetConstObjectMacro( MovingMask, MovingMaskType );

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void
  Initialize( void ) throw( ExceptionObject );

protected:
  MeshToMeshMetric();
  virtual ~MeshToMeshMetric(){};
  void
  PrintSelf( std::ostream & os, Indent indent ) const override;

  FixedMeshConstPointer  m_FixedMesh;
  MovingMeshConstPointer m_MovingMesh;

  mutable TransformPointer m_Transform;
  InterpolatorPointer      m_Interpolator;

  mutable FixedMaskPointer  m_FixedMask;
  mutable MovingMaskPointer m_MovingMask;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMeshToMeshMetric.hxx"
#endif

#endif
