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
#ifndef __itkReplaceDestinationPointsQuadEdgeMeshFilter_hxx
#define __itkReplaceDestinationPointsQuadEdgeMeshFilter_hxx

#include "itkReplaceDestinationPointsQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template < typename TInputMesh, typename TInputPointSet >
ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::ReplaceDestinationPointsQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );
}

template < typename TInputMesh, typename TInputPointSet >
ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::~ReplaceDestinationPointsQuadEdgeMeshFilter()
{}

template < typename TInputMesh, typename TInputPointSet >
void
ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro( "setting input Mesh to " << mesh );
  if ( mesh != static_cast< const InputMeshType * >( this->ProcessObject::GetInput( 0 ) ) )
  {
    this->ProcessObject::SetNthInput( 0, const_cast< InputMeshType * >( mesh ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TInputPointSet >
const typename ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::InputMeshType *
ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::GetInputMesh() const
{
  Self *                surrogate = const_cast< Self * >( this );
  const InputMeshType * referenceMesh = static_cast< const InputMeshType * >( surrogate->ProcessObject::GetInput( 0 ) );
  return referenceMesh;
}

template < typename TInputMesh, typename TInputPointSet >
void
ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::SetDestinationPoints(
  const InputPointSetType * destinationPointSet )
{
  itkDebugMacro( "setting input ReferenceMesh to " << destinationPointSet );
  if ( destinationPointSet != static_cast< const InputPointSetType * >( this->ProcessObject::GetInput( 1 ) ) )
  {
    this->ProcessObject::SetNthInput( 1, const_cast< InputPointSetType * >( destinationPointSet ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TInputPointSet >
const typename ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::InputPointSetType *
ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::GetDestinationPoints() const
{
  Self *                    surrogate = const_cast< Self * >( this );
  const InputPointSetType * destinationPointSet =
    static_cast< const InputPointSetType * >( surrogate->ProcessObject::GetInput( 1 ) );
  return destinationPointSet;
}

template < typename TInputMesh, typename TInputPointSet >
void
ReplaceDestinationPointsQuadEdgeMeshFilter< TInputMesh, TInputPointSet >::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  const InputPointSetType * inputPointSet = this->GetDestinationPoints();

  if ( !inputPointSet )
  {
    itkExceptionMacro( "Input PointSet is missing" );
  }

  using DestinationPointsContainer = typename InputPointSetType::PointsContainer;

  const DestinationPointsContainer * destinationPoints = inputPointSet->GetPoints();

  if ( !destinationPoints )
  {
    itkExceptionMacro( "Input PointSet has no points" );
  }

  OutputMeshType * outputMesh = this->GetOutput();

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  if ( destinationPoints->Size() != numberOfPoints )
  {
    itkExceptionMacro( "The PointSet does not have the same number of points as the Mesh" );
  }

  ProgressReporter progress( this, 0, numberOfPoints );

  using DestinationPointIterator = typename DestinationPointsContainer::ConstIterator;

  DestinationPointIterator destinationPointItr = destinationPoints->Begin();
  DestinationPointIterator destinationPointEnd = destinationPoints->End();

  using OutputPointContainer = typename OutputMeshType::PointsContainer;
  using OutputPointIterator = typename OutputPointContainer::Iterator;

  OutputPointContainer * outputPoints = outputMesh->GetPoints();

  OutputPointIterator outputPointItr = outputPoints->Begin();

  while ( destinationPointItr != destinationPointEnd )
  {
    OutputPointType & outputPoint = outputPointItr.Value();
    outputPoint.SetPoint( destinationPointItr.Value() );
    ++outputPointItr;
    ++destinationPointItr;
  }
}
} // end namespace itk

#endif
