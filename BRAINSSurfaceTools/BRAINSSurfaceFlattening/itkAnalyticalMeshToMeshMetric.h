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
#ifndef __itkAnalyticalMeshToMeshMetric_h
#define __itkAnalyticalMeshToMeshMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"
#include "itkMeshToMeshMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"

namespace itk
{
/** \class AnalyticalMeshToMeshMetric
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
template <class TFixedMesh, class TMovingMesh>
class AnalyticalMeshToMeshMetric :
  public         MeshToMeshMetric<TFixedMesh, TMovingMesh>
{
public:

  /** Standard class typedefs. */
  typedef AnalyticalMeshToMeshMetric                Self;
  typedef MeshToMeshMetric<TFixedMesh, TMovingMesh> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AnalyticalMeshToMeshMetric, MeshToMeshMetric);

  /** Types transferred from the base class */
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;
  typedef typename Superclass::InputPointType          InputPointType;
  typedef typename Superclass::OutputPointType         OutputPointType;

  typedef typename Superclass::MeasureType            MeasureType;
  typedef typename Superclass::DerivativeType         DerivativeType;
  typedef typename Superclass::FixedMeshType          FixedMeshType;
  typedef typename Superclass::MovingMeshType         MovingMeshType;
  typedef typename Superclass::FixedMeshConstPointer  FixedMeshConstPointer;
  typedef typename Superclass::MovingMeshConstPointer MovingMeshConstPointer;
  typedef typename Superclass::PointIterator          PointIterator;
  typedef typename Superclass::PointDataIterator      PointDataIterator;

  typedef typename Superclass::InterpolatorType InterpolatorType;

  typedef typename Superclass::RealDataType       RealDataType;
  typedef typename Superclass::DerivativeDataType DerivativeDataType;

  /** Constants for the pointset dimensions */
  itkStaticConstMacro(MovingMeshDimension, unsigned int,
                      Superclass::MovingMeshDimension);
  itkStaticConstMacro(FixedMeshDimension, unsigned int,
                      Superclass::FixedMeshDimension);

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters, DerivativeType & derivative ) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters, MeasureType& Value,
                              DerivativeType& Derivative ) const;

protected:
  AnalyticalMeshToMeshMetric();
  virtual ~AnalyticalMeshToMeshMetric()
  {
  };
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(AnalyticalMeshToMeshMetric);

  mutable unsigned int m_NumberOfPixelsCounted;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnalyticalMeshToMeshMetric.hxx"
#endif

#endif
