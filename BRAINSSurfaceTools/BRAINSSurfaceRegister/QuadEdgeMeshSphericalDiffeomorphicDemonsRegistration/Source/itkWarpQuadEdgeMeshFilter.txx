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
#ifndef __itkWarpQuadEdgeMeshFilter_txx
#define __itkWarpQuadEdgeMeshFilter_txx

#include "itkWarpQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TReferenceMesh, class TDeformationField>
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::WarpQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 3 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  // Setup default interpolator
  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  m_Interpolator =
    static_cast<InterpolatorType *>( interp.GetPointer() );
}

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::~WarpQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
void
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro("setting input mesh to " << mesh);
  if( mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
const typename
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>::InputMeshType
* WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::GetInputMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * inputMesh =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(0) );
  return inputMesh;
  }

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
void
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::SetReferenceMesh( const ReferenceMeshType * mesh )
{
  itkDebugMacro("setting input deformation mesh to " << mesh);
  if( mesh != static_cast<const ReferenceMeshType *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<ReferenceMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
const typename
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>::ReferenceMeshType
* WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::GetReferenceMesh() const
  {
  Self *                    surrogate = const_cast<Self *>( this );
  const ReferenceMeshType * referenceMesh =
    static_cast<const ReferenceMeshType *>( surrogate->ProcessObject::GetInput(1) );
  return referenceMesh;
  }

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
void
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::SetDeformationField( const DeformationFieldType * field )
{
  itkDebugMacro("setting input deformation field to " << field);
  if( field != static_cast<const DeformationFieldType *>(this->ProcessObject::GetInput( 2 ) ) )
    {
    this->ProcessObject::SetNthInput(2, const_cast<DeformationFieldType *>( field ) );
    this->Modified();
    }
}

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
const typename
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>::DeformationFieldType
* WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::GetDeformationField() const
  {
  Self *                       surrogate = const_cast<Self *>( this );
  const DeformationFieldType * deformationfield =
    static_cast<const DeformationFieldType *>( surrogate->ProcessObject::GetInput(2) );
  return deformationfield;
  }

template <class TInputMesh, class TReferenceMesh, class TDeformationField>
void
WarpQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDeformationField>
::GenerateData()
{
  // Copy the input mesh into the output mesh.

  this->CopyInputMeshToOutputMesh();

  OutputMeshType * outputMesh = this->GetOutput();

  //
  // Visit all nodes of the Mesh
  //

  OutputPointsContainerPointer points = outputMesh->GetPoints();

  if( points.IsNull() )
    {
    itkExceptionMacro("Mesh has NULL PointData");
    }

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  ProgressReporter progress(this, 0, numberOfPoints);

  OutputPointDataContainerPointer pointData = outputMesh->GetPointData();

  if( pointData.IsNull() )
    {
    pointData = OutputPointDataContainer::New();
    outputMesh->SetPointData( pointData );
    }

  pointData->Reserve( numberOfPoints );

  const ReferenceMeshType * referenceMesh = this->GetReferenceMesh();

  // Initialize the internal point locator structure
  this->m_Interpolator->SetInputMesh( referenceMesh );
  this->m_Interpolator->Initialize();

  //
  // deformation field: mesh with vector pixel value
  //
  const DeformationFieldType * deformationfield = this->GetDeformationField();

  const DisplacementVectorContainer * displacementVectors = deformationfield->GetPointData();

  typedef typename OutputMeshType::PointsContainer::ConstIterator PointIterator;
  typedef typename OutputMeshType::PointDataContainer::Iterator   PointDataIterator;

  PointIterator pointItr = points->Begin();
  PointIterator pointEnd = points->End();

  PointDataIterator pointDataItr = pointData->Begin();
  PointDataIterator pointDataEnd = pointData->End();

  typedef typename DeformationFieldType::PointDataContainer::ConstIterator DisplacementVectorIterator;

  DisplacementVectorIterator displacementVectorItr = displacementVectors->Begin();
  DisplacementVectorIterator displacementVectorEnd = displacementVectors->End();

  const unsigned int deformationfieldNumberOfPoints = deformationfield->GetNumberOfPoints();

  if( deformationfieldNumberOfPoints != numberOfPoints )
    {
    itkExceptionMacro("Input Fixed Mesh and Deformation Field have "
                      << "different number of points " << numberOfPoints
                      << " vs " << deformationfieldNumberOfPoints );
    }

  OutputPointType pointToEvaluate;
  while( pointItr != pointEnd && pointDataItr != pointDataEnd && displacementVectorItr != displacementVectorEnd )
    {
    pointToEvaluate.CastFrom( pointItr.Value() + displacementVectorItr.Value() );
    pointDataItr.Value() = this->m_Interpolator->Evaluate( pointToEvaluate );

    progress.CompletedPixel();

    ++pointItr;
    ++pointDataItr;
    ++displacementVectorItr;
    }
}
} // end namespace itk

#endif
