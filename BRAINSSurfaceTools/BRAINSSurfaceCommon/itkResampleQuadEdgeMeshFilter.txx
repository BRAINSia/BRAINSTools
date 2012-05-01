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
#ifndef __itkResampleQuadEdgeMeshFilter_txx
#define __itkResampleQuadEdgeMeshFilter_txx

#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkVersor.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TOutputMesh>
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::ResampleQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TOutputMesh>
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::~ResampleQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::SetReferenceMesh( const TOutputMesh * mesh )
{
  itkDebugMacro("setting input ReferenceMesh to " << mesh);
  if( mesh != static_cast<const TOutputMesh *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<TOutputMesh *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TOutputMesh>
const typename ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>::OutputMeshType
* ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::GetReferenceMesh() const
  {
  Self *                 surrogate = const_cast<Self *>( this );
  const OutputMeshType * referenceMesh =
    static_cast<const OutputMeshType *>(surrogate->ProcessObject::GetInput(1) );
  return referenceMesh;
  }

template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::GenerateData()
{
  // Copy the input mesh into the output mesh.
  this->CopyReferenceMeshToOutputMesh();

  OutputMeshPointer outputMesh = this->GetOutput();

  //
  // Visit all nodes of the Mesh
  //

  OutputPointsContainerPointer points = outputMesh->GetPoints();

  if( points.IsNull() )
    {
    itkExceptionMacro("Mesh has NULL Points");
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

  // Initialize the internal point locator structure
  this->m_Interpolator->SetInputMesh( this->GetInput() );
  this->m_Interpolator->Initialize();

  typedef typename OutputMeshType::PointsContainer::ConstIterator PointIterator;
  typedef typename OutputMeshType::PointDataContainer::Iterator   PointDataIterator;

  PointIterator pointItr = points->Begin();
  PointIterator pointEnd = points->End();

  PointDataIterator pointDataItr = pointData->Begin();
  PointDataIterator pointDataEnd = pointData->End();

  typedef typename TransformType::OutputPointType MappedPointType;

  OutputPointType inputPoint;
  OutputPointType pointToEvaluate;

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    inputPoint.CastFrom( pointItr.Value() );

    MappedPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    pointToEvaluate.CastFrom( transformedPoint );
    pointDataItr.Value() = this->m_Interpolator->Evaluate( pointToEvaluate );

    progress.CompletedPixel();

    ++pointItr;
    ++pointDataItr;
    }
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMesh()
{
  this->CopyReferenceMeshToOutputMeshGeometry();
  this->CopyReferenceMeshToOutputMeshFieldData();
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMeshGeometry()
{
  this->CopyReferenceMeshToOutputMeshPoints();
  this->CopyReferenceMeshToOutputMeshEdgeCells();
  this->CopyReferenceMeshToOutputMeshCells();
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMeshFieldData()
{
  this->CopyReferenceMeshToOutputMeshPointData();
  this->CopyReferenceMeshToOutputMeshCellData();
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMeshPoints()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshPoints( in, out );
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMeshEdgeCells()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshEdgeCells( in, out );
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMeshCells()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshCells( in, out );
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMeshPointData()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshPointData( in, out );
}

// ---------------------------------------------------------------------
template <class TInputMesh, class TOutputMesh>
void
ResampleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::CopyReferenceMeshToOutputMeshCellData()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshCellData( in, out );
}
} // end namespace itk

#endif
