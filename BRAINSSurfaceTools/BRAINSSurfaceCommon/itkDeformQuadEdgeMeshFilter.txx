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
#ifndef __itkDeformQuadEdgeMeshFilter_txx
#define __itkDeformQuadEdgeMeshFilter_txx

#include "itkDeformQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::DeformQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 3 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  this->m_Interpolator = InterpolatorType::New();
  this->m_Interpolator->SetUseNearestNeighborInterpolationAsBackup(true);

  this->m_SphereRadius = 1.0;
  this->m_SphereCenter.Fill(0.0);
}

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::~DeformQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro("setting input mesh to " << mesh);
  if( mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
const typename
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::InputMeshType
* DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::GetInputMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * inputMesh =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(0) );
  return inputMesh;
  }

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::SetReferenceMesh( const ReferenceMeshType * mesh )
{
  itkDebugMacro("setting input deformation mesh to " << mesh);
  if( mesh != static_cast<const ReferenceMeshType *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<ReferenceMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
const typename
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::ReferenceMeshType
* DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::GetReferenceMesh() const
  {
  Self *                    surrogate = const_cast<Self *>( this );
  const ReferenceMeshType * deformationMesh =
    static_cast<const ReferenceMeshType *>( surrogate->ProcessObject::GetInput(1) );
  return deformationMesh;
  }

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::SetDestinationPoints( const DestinationPointsType * points )
{
  itkDebugMacro("setting input destination points to " << points);
  if( points != static_cast<const DestinationPointsType *>(this->ProcessObject::GetInput( 2 ) ) )
    {
    this->ProcessObject::SetNthInput(2, const_cast<DestinationPointsType *>( points ) );
    this->Modified();
    }
}

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
const typename
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::DestinationPointsType
* DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::GetDestinationPoints() const
  {
  Self *                        surrogate = const_cast<Self *>( this );
  const DestinationPointsType * destinationPoints =
    static_cast<const DestinationPointsType *>( surrogate->ProcessObject::GetInput(2) );
  return destinationPoints;
  }

template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>
::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  const ReferenceMeshType * referenceMesh = this->GetReferenceMesh();

  OutputMeshType * outputMesh = this->GetOutput();

  const DestinationPointsType * destinationPoints = this->GetDestinationPoints();

  const DestinationPointsContainerType * destinationPointsContainer = destinationPoints->GetPoints();

  typedef typename OutputMeshType::PointsContainer OutputPointsContainer;

  OutputPointsContainer * outputPoints = outputMesh->GetPoints();

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  const unsigned int referenceNumberOfPoints = referenceMesh->GetNumberOfPoints();

  const unsigned int destinationNumberOfPoints = destinationPoints->GetNumberOfPoints();

  if( destinationNumberOfPoints != referenceNumberOfPoints )
    {
    itkExceptionMacro("Reference Mesh and Destination Points have "
                      << "different number of points " << referenceNumberOfPoints
                      << " vs " << destinationNumberOfPoints );
    }

  ProgressReporter progress(this, 0, numberOfPoints);

  this->m_Interpolator->SetSphereCenter( this->m_SphereCenter );

  this->m_Interpolator->SetInputMesh( referenceMesh );
  this->m_Interpolator->Initialize();

  typedef typename OutputPointsContainer::Iterator OutputPointIterator;

  OutputPointIterator outputPointItr = outputPoints->Begin();
  OutputPointIterator outputPointEnd = outputPoints->End();

  typedef typename InterpolatorType::PointType PointType;
  typedef typename PointType::VectorType       VectorType;

  PointType evaluatedPoint;
  PointType inputPoint;

  while( outputPointItr != outputPointEnd )
    {
    inputPoint.CastFrom( outputPointItr.Value() );
    this->m_Interpolator->Evaluate( destinationPointsContainer, inputPoint, evaluatedPoint );

    //
    //  Project point to sphere surface
    //
    VectorType vectorToCenter( evaluatedPoint - this->m_SphereCenter );

    const double radialDistance = vectorToCenter.GetNorm();
    vectorToCenter *= this->m_SphereRadius / radialDistance;
    evaluatedPoint = this->m_SphereCenter + vectorToCenter;

    //
    // SetPoint() must be used here instead of a simple assignment in order to
    // maintain the topology of the QuadEdgeMesh.
    //
    outputPointItr.Value().SetPoint( evaluatedPoint );

    progress.CompletedPixel();

    ++outputPointItr;
    }
}
} // end namespace itk

#endif
