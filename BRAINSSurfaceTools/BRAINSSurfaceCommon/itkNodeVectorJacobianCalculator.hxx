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
#ifndef __itkNodeVectorJacobianCalculator_hxx
#define __itkNodeVectorJacobianCalculator_hxx

#include "itkNodeScalarGradientCalculator.h"
#include "itkVersor.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputMesh, class TScalar>
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::NodeVectorJacobianCalculator()
{
  this->m_AreaList = AreaListType::New();
  this->m_PointAreaAccumulatorList = CoordRepListType::New();
  this->m_PointJacobianAccumulatorList = JacobianListType::New();

  this->m_SphereCenter.Fill( 0.0 );
  this->m_SphereRadius = 1.0;
}

/**
 * Destructor
 */
template <class TInputMesh, class TScalar>
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::~NodeVectorJacobianCalculator()
{
}

/**
 * Check inputs
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::Initialize()
{
  this->VerifyInputs();
  this->AllocateInternalContainers();
  this->ComputeAreaForAllCells();
}

/**
 * Check inputs
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::VerifyInputs() const
{
  if( this->m_InputMesh.IsNull() )
    {
    itkExceptionMacro(<< "NodeVectorJacobianCalculator Initialize  m_InputMesh is NULL.");
    }

  if( this->m_VectorContainer.IsNull() )
    {
    itkExceptionMacro(<< "NodeVectorJacobianCalculator Initialize  m_VectorContainer is NULL.");
    }

  if( this->m_BasisSystemList.IsNull() )
    {
    itkExceptionMacro(<< "NodeVectorJacobianCalculator Initialize  m_BasisSystemList is NULL.");
    }
}

/**
 * Allocate internal Containers
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::AllocateInternalContainers()
{
  this->m_AreaList->Reserve( this->m_InputMesh->GetNumberOfCells() );
  this->m_PointAreaAccumulatorList->Reserve( this->m_InputMesh->GetNumberOfPoints() );
  this->m_PointJacobianAccumulatorList->Reserve( this->m_InputMesh->GetNumberOfPoints() );

  this->m_AreaList->Squeeze();
  this->m_PointAreaAccumulatorList->Squeeze();
  this->m_PointJacobianAccumulatorList->Squeeze();
}

/**
 * Compute the function
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::Compute()
{
  this->SetContainersToNullValues();

  const CellsContainer * cells =  this->m_InputMesh->GetCells();

  CellsContainerConstIterator cellIterator = cells->Begin();
  CellsContainerConstIterator cellEnd = cells->End();

  AreaListConstIterator areaIterator = this->m_AreaList->Begin();

  BasisSystemListIterator basisSystemListIterator = m_BasisSystemList->Begin();

  const unsigned int numberOfVerticesInTriangle = 3;

  DerivativeType derivative;
  DerivativeType parallelTransportedDerivative[numberOfVerticesInTriangle];

  JacobianType   projectedJacobian;
  DerivativeType projectedDerivative[numberOfVerticesInTriangle];

  //
  // Look at all triangular cells, re-use the basis of each, and new scalar values.
  //
  ArrayType       pixelValue[numberOfVerticesInTriangle];
  PointIdentifier pointIds[numberOfVerticesInTriangle];
  PointType       point[numberOfVerticesInTriangle];

  while( cellIterator != cellEnd )
    {
    CellType * cellPointer = cellIterator.Value();

    // Consider current cell. Iterate through its points.

    PointIdIterator pointIdIterator = cellPointer->PointIdsBegin();
    PointIdIterator pointIdEnd = cellPointer->PointIdsEnd();

    unsigned int i = 0;
    while( pointIdIterator != pointIdEnd )
      {
      const PointIdentifier pointId = *pointIdIterator;
      pointIds[i] = pointId;
      pixelValue[i] = this->m_VectorContainer->GetElement( pointId );
      point[i] = this->m_InputMesh->GetPoint( pointId );
      i++;
      ++pointIdIterator;
      }

    const AreaType area = areaIterator.Value();
    //
    // contribute to accumulated area around each point of triangle
    //
    for( unsigned int ik = 0; ik < numberOfVerticesInTriangle; ik++ )
      {
      this->m_PointAreaAccumulatorList->ElementAt( pointIds[ik] ) += area;
      }

    const TriangleBasisSystemType & basisSystem = basisSystemListIterator.Value();

    const VectorType & u12 = basisSystem.GetVector(0);
    const VectorType & u32 = basisSystem.GetVector(1);

    //
    //  Compute the coordinates of the point at the center of the cell.
    //
    const double equalWeight = 1.0 / 3.0;
    PointType    cellCenter;
    cellCenter.SetToBarycentricCombination(
      point[0], point[1], point[2], equalWeight, equalWeight );

    VectorType vectorToCenter = cellCenter - this->m_SphereCenter;
    vectorToCenter *= this->m_SphereRadius / vectorToCenter.GetNorm();

    const PointType cellCenterProjectedInSphere = this->m_SphereCenter + vectorToCenter;

    //
    // Project derivative in the plane tangent to the sphere at the projected
    // point of the center cell
    //
    VectorType vectorToCenterNormalized = vectorToCenter;
    vectorToCenterNormalized.Normalize();
    for( unsigned int il = 0; il < numberOfVerticesInTriangle; il++ )
      {
      ArrayValueType pixelValue1 = pixelValue[0][il];
      ArrayValueType pixelValue2 = pixelValue[1][il];
      ArrayValueType pixelValue3 = pixelValue[2][il];

      InterpolatorType::GetDerivativeFromPixelsAndBasis(
        pixelValue1, pixelValue2, pixelValue3, u12, u32, derivative );

      const double radialComponent = derivative * vectorToCenterNormalized;
      for( unsigned int k = 0; k < MeshDimension; k++ )
        {
        projectedDerivative[il][k] = derivative[k] - vectorToCenterNormalized[k] * radialComponent;
        }

      projectedDerivative[il] *= area;
      }
    // Store at each vertex the value equal to triangle area x
    // Jacobian. Later, area-based weighting will use the total value
    // and divide by the sum of triangle areas about each vertex.
    // Begin by weighting contribution to point of this triangle by its area.
    for( unsigned int im = 0; im < numberOfVerticesInTriangle; im++ )
      {
      for( unsigned int k = 0; k < MeshDimension; k++ )
        {
        //
        // Parallel transport the derivative vector to each neighbor point
        //
        this->ParallelTransport( cellCenterProjectedInSphere, point[im],
                                 projectedDerivative[k], parallelTransportedDerivative[k] );
        for( unsigned int h = 0; h < MeshDimension; h++ )
          {
          projectedJacobian[k][h] = parallelTransportedDerivative[k][h];
          }
        }

      //
      // then accumulate them there.
      //
      this->m_PointJacobianAccumulatorList->ElementAt( pointIds[im] ) += projectedJacobian;
      }

    ++cellIterator;
    ++areaIterator;
    ++basisSystemListIterator;
    }

  this->NormalizeDerivativesByTotalArea();
}

/**
 * Initialize several containers with null values in all their elements.
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::SetContainersToNullValues()
{
  typedef typename DerivativeType::ValueType DerivativeValueType;

  JacobianType nullJacobian;
  nullJacobian.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

  JacobianListIterator jacobianItr = this->m_PointJacobianAccumulatorList->Begin();
  JacobianListIterator jacobianEnd = this->m_PointJacobianAccumulatorList->End();
  AreaListIterator     areaItr = this->m_PointAreaAccumulatorList->Begin();

  while( jacobianItr != jacobianEnd )
    {
    jacobianItr.Value() = nullJacobian;
    areaItr.Value() = NumericTraits<AreaType>::ZeroValue();
    ++jacobianItr;
    ++areaItr;
    }
}

/**
 * Compute the area of each cell and store them in a container
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::ComputeAreaForAllCells()
{
  const CellsContainer * cells =  this->m_InputMesh->GetCells();

  const unsigned int numberOfVerticesInTriangle = 3;
  PointType          point[numberOfVerticesInTriangle];

  AreaListIterator areaIterator = this->m_AreaList->Begin();

  CellsContainerConstIterator cellIterator = cells->Begin();
  CellsContainerConstIterator cellEnd = cells->End();

  while( cellIterator != cellEnd )
    {
    CellType * cellPointer = cellIterator.Value();

    // Consider current cell. Iterate through its points.
    PointIdIterator pointIdIterator = cellPointer->PointIdsBegin();
    PointIdIterator pointIdEnd = cellPointer->PointIdsEnd();

    unsigned int i = 0;
    while( pointIdIterator != pointIdEnd )
      {
      const PointIdentifier pointId = *pointIdIterator;
      point[i] = m_InputMesh->GetPoint( pointId );
      i++;
      ++pointIdIterator;
      }

    const AreaType area = TriangleHelperType::ComputeArea( point[0], point[1], point[2] );

    areaIterator.Value() = area;

    ++cellIterator;
    ++areaIterator;
    }
}

/**
 * Compute the function
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::NormalizeDerivativesByTotalArea()
{
  JacobianListIterator  jacobianItr = this->m_PointJacobianAccumulatorList->Begin();
  JacobianListIterator  jacobianEnd = this->m_PointJacobianAccumulatorList->End();
  AreaListConstIterator areaItr = this->m_PointAreaAccumulatorList->Begin();

  while( jacobianItr != jacobianEnd )
    {
    jacobianItr.Value() /= areaItr.Value();
    ++jacobianItr;
    ++areaItr;
    }
}

template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::ParallelTransport(
  const PointType src, const PointType dst,
  const DerivativeType & inputVector,
  DerivativeType & transportedVector ) const
{
  VectorType vsrc = src - this->m_SphereCenter;
  VectorType vdst = dst - this->m_SphereCenter;

  VectorType axis = CrossProduct( vsrc, vdst );

  const double scaledSinus   = axis.GetNorm();
  const double scaledCosinus = vsrc * vdst;

  double angle = vcl_atan2( scaledSinus, scaledCosinus );

  typedef Versor<double> VersorType;

  VersorType versor;
  versor.Set( axis, angle );

  transportedVector = versor.Transform( inputVector );
}

/**
 * Compute the function
 */
template <class TInputMesh, class TScalar>
typename NodeVectorJacobianCalculator<TInputMesh, TScalar>::OutputType
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::Evaluate( const InputType & input  ) const
{
  return this->m_PointJacobianAccumulatorList->ElementAt( input );
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputMesh, class TScalar>
void
NodeVectorJacobianCalculator<TInputMesh, TScalar>
::PrintSelf( std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Input mesh = " << this->m_InputMesh.GetPointer() << std::endl;
  os << indent << "Data Container = " << this->m_VectorContainer.GetPointer() << std::endl;
  os << indent << "Basis System List = " << this->m_BasisSystemList.GetPointer() << std::endl;
  os << indent << "Point Area Accumulator = " << this->m_PointAreaAccumulatorList.GetPointer() << std::endl;
  os << indent << "Point Derivative Accumulator = " << this->m_PointJacobianAccumulatorList.GetPointer() << std::endl;
}

#include <iostream>
} // end namespace itk

#endif
