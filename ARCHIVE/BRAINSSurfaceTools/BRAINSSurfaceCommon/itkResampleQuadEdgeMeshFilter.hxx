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
#ifndef __itkResampleQuadEdgeMeshFilter_hxx
#define __itkResampleQuadEdgeMeshFilter_hxx

#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkVersor.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template < typename TInputMesh, typename TOutputMesh >
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::ResampleQuadEdgeMeshFilter()
{}

template < typename TInputMesh, typename TOutputMesh >
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::~ResampleQuadEdgeMeshFilter()
{}

template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::SetReferenceMesh( const TOutputMesh * mesh )
{
  itkDebugMacro( "setting input ReferenceMesh to " << mesh );
  if ( mesh != static_cast< const TOutputMesh * >( this->ProcessObject::GetInput( 1 ) ) )
  {
    this->ProcessObject::SetNthInput( 1, const_cast< TOutputMesh * >( mesh ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TOutputMesh >
const typename ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::OutputMeshType *
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::GetReferenceMesh() const
{
  Self *                 surrogate = const_cast< Self * >( this );
  const OutputMeshType * referenceMesh =
    static_cast< const OutputMeshType * >( surrogate->ProcessObject::GetInput( 1 ) );
  return referenceMesh;
}

template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::GenerateData()
{
  // Copy the input mesh into the output mesh.
  this->CopyReferenceMeshToOutputMesh();

  OutputMeshPointer outputMesh = this->GetOutput();

  //
  // Visit all nodes of the Mesh
  //

  OutputPointsContainerPointer points = outputMesh->GetPoints();

  if ( points.IsNull() )
  {
    itkExceptionMacro( "Mesh has NULL Points" );
  }

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  ProgressReporter progress( this, 0, numberOfPoints );

  OutputPointDataContainerPointer pointData = outputMesh->GetPointData();

  if ( pointData.IsNull() )
  {
    pointData = OutputPointDataContainer::New();
    outputMesh->SetPointData( pointData );
  }

  pointData->Reserve( numberOfPoints );

  // Initialize the internal point locator structure
  this->m_Interpolator->SetInputMesh( this->GetInput() );
  this->m_Interpolator->Initialize();

  using PointIterator = typename OutputMeshType::PointsContainer::ConstIterator;
  using PointDataIterator = typename OutputMeshType::PointDataContainer::Iterator;

  PointIterator pointItr = points->Begin();
  PointIterator pointEnd = points->End();

  PointDataIterator pointDataItr = pointData->Begin();
  PointDataIterator pointDataEnd = pointData->End();

  using MappedPointType = typename TransformType::OutputPointType;

  OutputPointType inputPoint;
  OutputPointType pointToEvaluate;

  while ( pointItr != pointEnd && pointDataItr != pointDataEnd )
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
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMesh()
{
  this->CopyReferenceMeshToOutputMeshGeometry();
  this->CopyReferenceMeshToOutputMeshFieldData();
}

// ---------------------------------------------------------------------
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMeshGeometry()
{
  this->CopyReferenceMeshToOutputMeshPoints();
  this->CopyReferenceMeshToOutputMeshEdgeCells();
  this->CopyReferenceMeshToOutputMeshCells();
}

// ---------------------------------------------------------------------
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMeshFieldData()
{
  this->CopyReferenceMeshToOutputMeshPointData();
  this->CopyReferenceMeshToOutputMeshCellData();
}

// ---------------------------------------------------------------------
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMeshPoints()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshPoints( in, out );
}

// ---------------------------------------------------------------------
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMeshEdgeCells()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshEdgeCells( in, out );
}

// ---------------------------------------------------------------------
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMeshCells()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshCells( in, out );
}

// ---------------------------------------------------------------------
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMeshPointData()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshPointData( in, out );
}

// ---------------------------------------------------------------------
template < typename TInputMesh, typename TOutputMesh >
void
ResampleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::CopyReferenceMeshToOutputMeshCellData()
{
  const OutputMeshType * in = this->GetReferenceMesh();
  OutputMeshType *       out = this->GetOutput();

  CopyMeshToMeshCellData( in, out );
}
} // end namespace itk

#endif
