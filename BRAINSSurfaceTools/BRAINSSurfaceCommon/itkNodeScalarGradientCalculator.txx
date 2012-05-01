/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculator.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNodeScalarGradientCalculator_txx
#define __itkNodeScalarGradientCalculator_txx

#include "itkNodeScalarGradientCalculator.h"
#include "itkVersor.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputMesh, class TScalar>
NodeScalarGradientCalculator<TInputMesh, TScalar>
::NodeScalarGradientCalculator()
{
  this->m_AreaList = AreaListType::New();
  this->m_PointAreaAccumulatorList = CoordRepListType::New();
  this->m_PointDerivativeAccumulatorList = DerivativeListType::New();

  this->m_SphereCenter.Fill( 0.0 );
  this->m_SphereRadius = 1.0;
}

/**
 * Destructor
 */
template <class TInputMesh, class TScalar>
NodeScalarGradientCalculator<TInputMesh, TScalar>
::~NodeScalarGradientCalculator()
{
}

/**
 * Check inputs
 */
template <class TInputMesh, class TScalar>
void
NodeScalarGradientCalculator<TInputMesh, TScalar>
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
NodeScalarGradientCalculator<TInputMesh, TScalar>
::VerifyInputs() const
{
  if( this->m_InputMesh.IsNull() )
    {
    itkExceptionMacro(<< "NodeScalarGradientCalculator Initialize  m_InputMesh is NULL.");
    }

  if( this->m_DataContainer.IsNull() )
    {
    itkExceptionMacro(<< "NodeScalarGradientCalculator Initialize  m_DataContainer is NULL.");
    }

  if( this->m_BasisSystemList.IsNull() )
    {
    itkExceptionMacro(<< "NodeScalarGradientCalculator Initialize  m_BasisSystemList is NULL.");
    }
}

/**
 * Allocate internal Containers
 */
template <class TInputMesh, class TScalar>
void
NodeScalarGradientCalculator<TInputMesh, TScalar>
::AllocateInternalContainers()
{
  this->m_AreaList->Reserve( this->m_InputMesh->GetNumberOfCells() );
  this->m_PointAreaAccumulatorList->Reserve( this->m_InputMesh->GetNumberOfPoints() );
  this->m_PointDerivativeAccumulatorList->Reserve( this->m_InputMesh->GetNumberOfPoints() );

  this->m_AreaList->Squeeze();
  this->m_PointAreaAccumulatorList->Squeeze();
  this->m_PointDerivativeAccumulatorList->Squeeze();
}

/**
 * Compute the function
 */
template <class TInputMesh, class TScalar>
void
NodeScalarGradientCalculator<TInputMesh, TScalar>
::Compute()
{
  this->SetContainersToNullValues();

  const CellsContainer * cells =  this->m_InputMesh->GetCells();

  CellsContainerConstIterator cellIterator = cells->Begin();
  CellsContainerConstIterator cellEnd = cells->End();

  AreaListConstIterator areaIterator = this->m_AreaList->Begin();

  BasisSystemListIterator basisSystemListIterator = m_BasisSystemList->Begin();

  DerivativeType derivative;
  DerivativeType parallelTransportedDerivative;

  //
  // Look at all triangular cells, re-use the basis of each, and new scalar values.
  //
  const unsigned int numberOfVerticesInTriangle = 3;
  PixelType          pixelValue[numberOfVerticesInTriangle];
  PointIdentifier    pointIds[numberOfVerticesInTriangle];
  PointType          point[numberOfVerticesInTriangle];

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
      pixelValue[i] = this->m_DataContainer->GetElement( pointId );
      point[i] = this->m_InputMesh->GetPoint( pointId );
      i++;
      ++pointIdIterator;
      }

    const AreaType area = areaIterator.Value();
    //
    // contribute to accumulated area around each point of triangle
    //
    for( unsigned int ii = 0; ii < numberOfVerticesInTriangle; ii++ )
      {
      this->m_PointAreaAccumulatorList->ElementAt( pointIds[ii] ) += area;
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

    InterpolatorType::GetDerivativeFromPixelsAndBasis(
      pixelValue[0], pixelValue[1], pixelValue[2], u12, u32, derivative);

    const double radialComponent = derivative * vectorToCenterNormalized;

    DerivativeType projectedDerivative;
    for( unsigned int k = 0; k < MeshDimension; k++ )
      {
      projectedDerivative[k] = derivative[k] - vectorToCenterNormalized[k] * radialComponent;
      }

    // Store at each vertex the value equal to triangle area x
    // derivative. Later, area-based weighting will use the total value
    // and divide by the sum of triangle areas about each vertex.
    // Begin by weighting contribution to point of this triangle by its area.

    projectedDerivative *= area;
    for( unsigned int ij = 0; ij < numberOfVerticesInTriangle; ij++ )
      {
      //
      // Parallel transport the derivative vector to each neighbor point
      //
      this->ParalelTransport( cellCenterProjectedInSphere, point[ij],
                              projectedDerivative, parallelTransportedDerivative );

      //
      // then accumulate them there.
      //
      this->m_PointDerivativeAccumulatorList->ElementAt( pointIds[ij] ) += parallelTransportedDerivative;
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
NodeScalarGradientCalculator<TInputMesh, TScalar>
::SetContainersToNullValues()
{
  typedef typename DerivativeType::ValueType DerivativeValueType;

  DerivativeType nullDerivative;
  nullDerivative.Fill( NumericTraits<DerivativeValueType>::Zero );

  DerivativeListIterator derivativeItr = this->m_PointDerivativeAccumulatorList->Begin();
  DerivativeListIterator derivativeEnd = this->m_PointDerivativeAccumulatorList->End();
  AreaListIterator       areaItr = this->m_PointAreaAccumulatorList->Begin();

  while( derivativeItr != derivativeEnd )
    {
    derivativeItr.Value() = nullDerivative;
    areaItr.Value() = NumericTraits<AreaType>::Zero;
    ++derivativeItr;
    ++areaItr;
    }
}

/**
 * Compute the area of each cell and store them in a container
 */
template <class TInputMesh, class TScalar>
void
NodeScalarGradientCalculator<TInputMesh, TScalar>
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

    const AreaType area = TriangleType::ComputeArea( point[0], point[1], point[2] );

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
NodeScalarGradientCalculator<TInputMesh, TScalar>
::NormalizeDerivativesByTotalArea()
{
  DerivativeListIterator derivativeItr = this->m_PointDerivativeAccumulatorList->Begin();
  DerivativeListIterator derivativeEnd = this->m_PointDerivativeAccumulatorList->End();
  AreaListConstIterator  areaItr = this->m_PointAreaAccumulatorList->Begin();

  while( derivativeItr != derivativeEnd )
    {
    derivativeItr.Value() /= areaItr.Value();
    ++derivativeItr;
    ++areaItr;
    }
}

template <class TInputMesh, class TScalar>
void
NodeScalarGradientCalculator<TInputMesh, TScalar>
::ParalelTransport(
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
typename NodeScalarGradientCalculator<TInputMesh, TScalar>::OutputType
NodeScalarGradientCalculator<TInputMesh, TScalar>
::Evaluate( const InputType & input  ) const
{
  return this->m_PointDerivativeAccumulatorList->ElementAt( input );
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputMesh, class TScalar>
void
NodeScalarGradientCalculator<TInputMesh, TScalar>
::PrintSelf( std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Input mesh = " << this->m_InputMesh.GetPointer() << std::endl;
  os << indent << "Data Container = " << this->m_DataContainer.GetPointer() << std::endl;
  os << indent << "Basis System List = " << this->m_BasisSystemList.GetPointer() << std::endl;
  os << indent << "Point Area Accumulator = " << this->m_PointAreaAccumulatorList.GetPointer() << std::endl;
  os << indent << "Point Derivative Accumulator = " << this->m_PointDerivativeAccumulatorList.GetPointer() << std::endl;
}

#include <iostream>
} // end namespace itk

#endif
