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
#ifndef __itkQuadEdgeMeshAddScalarsFilter_txx
#define __itkQuadEdgeMeshAddScalarsFilter_txx

#include "itkQuadEdgeMeshAddScalarsFilter.h"
#include "itkProgressReporter.h"

namespace itk
{
template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::QuadEdgeMeshAddScalarsFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );
}

template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::~QuadEdgeMeshAddScalarsFilter()
{}

template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
void
QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::SetInput1( const InputMeshType1 * mesh )
{
  itkDebugMacro( "setting input1 to " << mesh );
  if ( mesh != static_cast< const InputMeshType1 * >( this->ProcessObject::GetInput( 0 ) ) )
  {
    this->ProcessObject::SetNthInput( 0, const_cast< InputMeshType1 * >( mesh ) );
    this->Modified();
  }
}

template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
void
QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::SetInput2( const InputMeshType2 * mesh )
{
  itkDebugMacro( "setting input2 to " << mesh );
  if ( mesh != static_cast< const InputMeshType2 * >( this->ProcessObject::GetInput( 1 ) ) )
  {
    this->ProcessObject::SetNthInput( 1, const_cast< InputMeshType2 * >( mesh ) );
    this->Modified();
  }
}

template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
const typename QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::InputMeshType1 *
QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::GetInput1() const
{
  Self *                 surrogate = const_cast< Self * >( this );
  const InputMeshType1 * inputMesh = static_cast< const InputMeshType1 * >( surrogate->ProcessObject::GetInput( 0 ) );
  return inputMesh;
}

template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
const typename QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::InputMeshType2 *
QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::GetInput2() const
{
  Self *                 surrogate = const_cast< Self * >( this );
  const InputMeshType2 * inputMesh = static_cast< const InputMeshType2 * >( surrogate->ProcessObject::GetInput( 1 ) );
  return inputMesh;
}

template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
void
QuadEdgeMeshAddScalarsFilter< TInputMesh1, TInputMesh2, TOutputMesh >::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  const InputMeshType1 * inputMesh1 = this->GetInput1();
  const InputMeshType2 * inputMesh2 = this->GetInput2();

  OutputMeshType * outputMesh = this->GetOutput();

  const unsigned int numberOfPoints1 = inputMesh1->GetNumberOfPoints();
  const unsigned int numberOfPoints2 = inputMesh2->GetNumberOfPoints();

  if ( numberOfPoints1 != numberOfPoints2 )
  {
    itkExceptionMacro( "The inputs does not have the same number of points" );
  }

  using DisplacementVectorContainer1 = typename InputMeshType1::PointDataContainer;
  using DisplacementVectorContainer2 = typename InputMeshType2::PointDataContainer;
  using OutputDisplacementVectorContainer = typename OutputMeshType::PointDataContainer;

  const DisplacementVectorContainer1 * displacementVectors1 = inputMesh1->GetPointData();
  const DisplacementVectorContainer2 * displacementVectors2 = inputMesh2->GetPointData();

  OutputDisplacementVectorContainer * outputDisplacementVectors = outputMesh->GetPointData();

  using OutputDisplacementVectorIterator = typename OutputDisplacementVectorContainer::Iterator;

  OutputDisplacementVectorIterator outputItr = outputDisplacementVectors->Begin();
  OutputDisplacementVectorIterator outputEnd = outputDisplacementVectors->End();

  using InputIterator1 = typename DisplacementVectorContainer1::ConstIterator;
  using InputIterator2 = typename DisplacementVectorContainer2::ConstIterator;

  InputIterator1 inputItr1 = displacementVectors1->Begin();
  InputIterator2 inputItr2 = displacementVectors2->Begin();

  while ( outputItr != outputEnd )
  {
    outputItr.Value() = inputItr1.Value() + inputItr2.Value();

    ++outputItr;
    ++inputItr1;
    ++inputItr2;
  }
}
} // end namespace itk

#endif
