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
#ifndef __itkPiecewiseRescaleQuadEdgeMeshFilter_txx
#define __itkPiecewiseRescaleQuadEdgeMeshFilter_txx

#include "itkPiecewiseRescaleQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraits.h"

namespace itk
{
template < typename TInputMesh, typename TOutputMesh >
PiecewiseRescaleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::PiecewiseRescaleQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  this->m_OutputMaximum = NumericTraits< OutputPixelType >::max();
  this->m_OutputMinimum = NumericTraits< OutputPixelType >::NonpositiveMin();

  this->m_cValue = NumericTraits< InputPixelType >::ZeroValue();
  this->m_InputMaximum = NumericTraits< InputPixelType >::min();
  this->m_InputMinimum = NumericTraits< InputPixelType >::max();
}

template < typename TInputMesh, typename TOutputMesh >
PiecewiseRescaleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::~PiecewiseRescaleQuadEdgeMeshFilter()
{}

template < typename TInputMesh, typename TOutputMesh >
void
PiecewiseRescaleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro( "setting input mesh to " << mesh );
  if ( mesh != static_cast< const InputMeshType * >( this->ProcessObject::GetInput( 0 ) ) )
  {
    this->ProcessObject::SetNthInput( 0, const_cast< InputMeshType * >( mesh ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TOutputMesh >
const typename PiecewiseRescaleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::InputMeshType *
PiecewiseRescaleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::GetInputMesh() const
{
  Self *                surrogate = const_cast< Self * >( this );
  const InputMeshType * inputMesh = static_cast< const InputMeshType * >( surrogate->ProcessObject::GetInput( 0 ) );

  return inputMesh;
}

template < typename TInputMesh, typename TOutputMesh >
void
PiecewiseRescaleQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::GenerateData()
{
  // Copy the input mesh into the output mesh.
  this->CopyInputMeshToOutputMesh();

  OutputMeshType * outputMesh = this->GetOutput();

  const InputMeshType * inputMesh = this->GetInputMesh();

  //
  // Visit all nodes of the input mesh
  //
  InputPointDataContainerConstPointer inputPointData = inputMesh->GetPointData();
  using InputPointDataIterator = typename InputPointDataContainer::ConstIterator;

  InputPointDataIterator inputPointDataItr = inputPointData->Begin();
  InputPointDataIterator inputPointDataEnd = inputPointData->End();

  while ( inputPointDataItr != inputPointDataEnd )
  {
    if ( inputPointDataItr.Value() > this->m_InputMaximum )
    {
      this->m_InputMaximum = inputPointDataItr.Value();
    }

    if ( inputPointDataItr.Value() < this->m_InputMinimum )
    {
      this->m_InputMinimum = inputPointDataItr.Value();
    }

    ++inputPointDataItr;
  }

  // verify the cValue
  if ( m_cValue > m_InputMaximum || m_cValue < m_InputMinimum )
  {
    itkExceptionMacro( "cValue should be in [InputMinimum, InputMaximum]!" );
  }

  if ( ( this->m_InputMinimum != this->m_cValue ) && ( this->m_InputMaximum != this->m_cValue ) )
  {
    this->m_Scale_a = ( static_cast< double >( this->m_cValue ) - static_cast< double >( this->m_OutputMinimum ) ) /
                      ( static_cast< double >( this->m_cValue ) - static_cast< double >( this->m_InputMinimum ) );

    this->m_Scale_b = ( static_cast< double >( this->m_OutputMaximum ) - static_cast< double >( this->m_cValue ) ) /
                      ( static_cast< double >( this->m_InputMaximum ) - static_cast< double >( this->m_cValue ) );
  }
  else if ( this->m_InputMinimum != this->m_InputMaximum )
  {
    this->m_Scale_a =
      ( static_cast< double >( this->m_OutputMaximum ) - static_cast< double >( this->m_OutputMinimum ) ) /
      ( static_cast< double >( this->m_InputMaximum ) - static_cast< double >( this->m_InputMinimum ) );
    this->m_Scale_b = this->m_Scale_a;
  }
  else if ( this->m_InputMaximum != NumericTraits< InputPixelType >::ZeroValue() )
  {
    this->m_Scale_a =
      ( static_cast< double >( this->m_OutputMaximum ) - static_cast< double >( this->m_OutputMinimum ) ) /
      static_cast< double >( this->m_InputMaximum );
    this->m_Scale_b = this->m_Scale_a;
  }
  else
  {
    this->m_Scale_a = 0.0;
    this->m_Scale_b = 0.0;
  }

  //
  // Visit all nodes of the output mesh
  // and change the scalars
  OutputPointDataContainer * pointData = outputMesh->GetPointData();

  using PointDataIterator = typename OutputPointDataContainer::Iterator;

  PointDataIterator pointDataItr = pointData->Begin();
  PointDataIterator pointDataEnd = pointData->End();

  while ( pointDataItr != pointDataEnd )
  {
    if ( pointDataItr.Value() < static_cast< OutputPixelType >( m_cValue ) )
    {
      pointDataItr.Value() = static_cast< OutputPixelType >(
        this->m_OutputMinimum + ( ( pointDataItr.Value() - this->m_InputMinimum ) * this->m_Scale_a ) );
    }
    else if ( pointDataItr.Value() > static_cast< OutputPixelType >( m_cValue ) )
    {
      pointDataItr.Value() = static_cast< OutputPixelType >(
        this->m_cValue + ( ( pointDataItr.Value() - this->m_cValue ) * this->m_Scale_b ) );
    }
    else
    {
      pointDataItr.Value() = static_cast< OutputPixelType >( m_cValue );
    }

    ++pointDataItr;
  }
}
} // end namespace itk

#endif
