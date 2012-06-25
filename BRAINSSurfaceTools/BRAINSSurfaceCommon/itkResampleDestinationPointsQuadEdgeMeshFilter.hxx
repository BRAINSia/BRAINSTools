/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculator.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkResampleDestinationPointsQuadEdgeMeshFilter_hxx
#define __itkResampleDestinationPointsQuadEdgeMeshFilter_hxx

#include "itkResampleDestinationPointsQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkIdentityTransform.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::ResampleDestinationPointsQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 3 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputPointSetType::New() );

  this->m_Interpolator = InterpolatorType::New();
  // this->m_Interpolator->SetUseNearestNeighborInterpolationAsBackup(true);

  this->m_Transform = itk::IdentityTransform<double>::New();

  this->m_SphereCenter.Fill( 0.0 );
  this->m_SphereRadius = 1.0;
}

template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::~ResampleDestinationPointsQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
void
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::SetFixedMesh( const FixedMeshType * mesh )
{
  itkDebugMacro("setting input FixedMesh to " << mesh);
  if( mesh != static_cast<const FixedMeshType *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<FixedMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
const typename
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>::FixedMeshType
* ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::GetFixedMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const FixedMeshType * referenceMesh =
    static_cast<const FixedMeshType *>( surrogate->ProcessObject::GetInput(1) );
  return referenceMesh;
  }

template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
void
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::SetReferenceMesh( const ReferenceMeshType * mesh )
{
  itkDebugMacro("setting input ReferenceMesh to " << mesh);
  if( mesh != static_cast<const ReferenceMeshType *>(this->ProcessObject::GetInput( 2 ) ) )
    {
    this->ProcessObject::SetNthInput(2, const_cast<ReferenceMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
const typename
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>::ReferenceMeshType
* ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::GetReferenceMesh() const
  {
  Self *                    surrogate = const_cast<Self *>( this );
  const ReferenceMeshType * referenceMesh =
    static_cast<const ReferenceMeshType *>( surrogate->ProcessObject::GetInput(2) );
  return referenceMesh;
  }

template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
void
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::ProjectPointToSphereSurface( OutputPointType & point ) const
{
  typedef typename OutputPointType::VectorType VectorType;

  VectorType vectorToCenter( point - this->m_SphereCenter );

  const double radialDistance = vectorToCenter.GetNorm();
  vectorToCenter *= this->m_SphereRadius / radialDistance;
  point = this->m_SphereCenter + vectorToCenter;
}

template <class TInputMesh, class TFixedMesh, class TReferenceMesh, class TOutputMesh>
void
ResampleDestinationPointsQuadEdgeMeshFilter<TInputMesh, TFixedMesh, TReferenceMesh, TOutputMesh>
::GenerateData()
{
  const ReferenceMeshType * referenceMesh = this->GetReferenceMesh();

  if( !referenceMesh )
    {
    itkExceptionMacro("Reference mesh is missing");
    }

  const FixedMeshType * fixedMesh = this->GetFixedMesh();

  if( !fixedMesh )
    {
    itkExceptionMacro("Fixed mesh is missing");
    }

  const InputPointSetType * inputPointSet = this->GetInput();

  if( !inputPointSet )
    {
    itkExceptionMacro("Input PointSet is missing");
    }

  const InputPointsContainer * inputPoints = inputPointSet->GetPoints();

  if( !inputPoints )
    {
    itkExceptionMacro("Input PointSet has no points");
    }

  const unsigned int numberOfPoints = referenceMesh->GetNumberOfPoints();

  OutputPointSetPointer outputPointSet = this->GetOutput();

  OutputPointsContainerPointer points = outputPointSet->GetPoints();

  if( points.IsNull() || ( points->Size() != numberOfPoints ) )
    {
    // We need to reallocate the points container for the output.
    points = OutputPointsContainer::New();
    points->Reserve( numberOfPoints );
    outputPointSet->SetPoints( points );
    }

  ProgressReporter progress(this, 0, numberOfPoints);

  this->m_Interpolator->SetSphereCenter( this->m_SphereCenter );

  this->m_Interpolator->SetInputMesh( fixedMesh );
  this->m_Interpolator->Initialize();

  const ReferencePointsContainer * referencePoints = referenceMesh->GetPoints();

  ReferencePointsContainerConstIterator referenceItr = referencePoints->Begin();
  ReferencePointsContainerConstIterator referenceEnd = referencePoints->End();

  OutputPointsContainerIterator outputPointItr = points->Begin();

  typedef typename TransformType::InputPointType InputPointType;
  typedef typename InterpolatorType::PointType   TransformInputPointType;

  InputPointType inputPoint;

  TransformInputPointType pointToEvaluate;
  TransformInputPointType evaluatedPoint;

  OutputPointType resultingPoint;

  while( referenceItr != referenceEnd )
    {
    inputPoint.CastFrom( referenceItr.Value() );

    pointToEvaluate.CastFrom( this->m_Transform->TransformPoint( inputPoint ) );

    this->m_Interpolator->Evaluate( inputPoints, pointToEvaluate, evaluatedPoint );

    resultingPoint.CastFrom( evaluatedPoint );

    this->ProjectPointToSphereSurface( resultingPoint );

    outputPointItr.Value() = resultingPoint;

    progress.CompletedPixel();

    ++outputPointItr;
    ++referenceItr;
    }
}
} // end namespace itk

#endif
