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
#ifndef __itkAssignScalarValuesQuadEdgeMeshFilter_hxx
#define __itkAssignScalarValuesQuadEdgeMeshFilter_hxx

#include "itkAssignScalarValuesQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template < typename TInputMesh, typename TSourceMesh, typename TOutputMesh >
AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::AssignScalarValuesQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );
}

template < typename TInputMesh, typename TSourceMesh, typename TOutputMesh >
AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::~AssignScalarValuesQuadEdgeMeshFilter()
{}

template < typename TInputMesh, typename TSourceMesh, typename TOutputMesh >
void
AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::SetSourceMesh(
  const SourceMeshType * mesh )
{
  itkDebugMacro( "setting source Mesh to " << mesh );
  if ( mesh != static_cast< const SourceMeshType * >( this->ProcessObject::GetInput( 0 ) ) )
  {
    this->ProcessObject::SetNthInput( 1, const_cast< SourceMeshType * >( mesh ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TSourceMesh, typename TOutputMesh >
const typename AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::SourceMeshType *
AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::GetSourceMesh() const
{
  Self *                 surrogate = const_cast< Self * >( this );
  const SourceMeshType * sourceMesh = static_cast< const SourceMeshType * >( surrogate->ProcessObject::GetInput( 1 ) );
  return sourceMesh;
}

template < typename TInputMesh, typename TSourceMesh, typename TOutputMesh >
void
AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro( "setting input Mesh to " << mesh );
  if ( mesh != static_cast< const InputMeshType * >( this->ProcessObject::GetInput( 1 ) ) )
  {
    this->ProcessObject::SetNthInput( 0, const_cast< InputMeshType * >( mesh ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TSourceMesh, typename TOutputMesh >
const typename AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::InputMeshType *
AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::GetInputMesh() const
{
  Self *                surrogate = const_cast< Self * >( this );
  const InputMeshType * inputMesh = static_cast< const InputMeshType * >( surrogate->ProcessObject::GetInput( 0 ) );
  return inputMesh;
}

template < typename TInputMesh, typename TSourceMesh, typename TOutputMesh >
void
AssignScalarValuesQuadEdgeMeshFilter< TInputMesh, TSourceMesh, TOutputMesh >::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  const SourceMeshType * sourceMesh = this->GetSourceMesh();

  const SourcePointDataContainer * sourcePointData = sourceMesh->GetPointData();

  OutputPointDataContainerPointer outputPointData = this->GetOutput()->GetPointData();

  if ( outputPointData.IsNull() )
  {
    outputPointData = OutputPointDataContainer::New();
    this->GetOutput()->SetPointData( outputPointData );
  }

  const unsigned int numberOfNodes = sourceMesh->GetNumberOfPoints();

  outputPointData->Reserve( numberOfNodes );

  if ( !sourcePointData )
  {
    itkExceptionMacro( "Source PointData is missing" );
  }

  OutputPointDataContainerIterator outputDataItr = outputPointData->Begin();

  using SourcePointDataIterator = typename SourcePointDataContainer::ConstIterator;
  SourcePointDataIterator sourceDataItr = sourcePointData->Begin();
  SourcePointDataIterator sourceDataEnd = sourcePointData->End();

  while ( sourceDataItr != sourceDataEnd )
  {
    outputDataItr.Value() = sourceDataItr.Value();

    ++outputDataItr;
    ++sourceDataItr;
  }
}
} // end namespace itk

#endif
