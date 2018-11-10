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
#ifndef __itkQuadEdgeMeshGenerateDeformationFieldFilter_hxx
#define __itkQuadEdgeMeshGenerateDeformationFieldFilter_hxx

#include "itkQuadEdgeMeshGenerateDeformationFieldFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <typename TInputMesh, typename TInputPointSet, typename TOutputMesh>
QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>
::QuadEdgeMeshGenerateDeformationFieldFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  this->m_SphereCenter.Fill( 0.0 );
  this->m_SphereRadius = 1.0;
}

template <typename TInputMesh, typename TInputPointSet, typename TOutputMesh>
QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>
::~QuadEdgeMeshGenerateDeformationFieldFilter()
{
}

template <typename TInputMesh, typename TInputPointSet, typename TOutputMesh>
void
QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>
::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro("setting input Mesh to " << mesh);
  if( mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
    this->Modified();
    }
}

template <typename TInputMesh, typename TInputPointSet, typename TOutputMesh>
const typename
QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>::InputMeshType
* QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>
::GetInputMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * referenceMesh =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(0) );
  return referenceMesh;
  }

template <typename TInputMesh, typename TInputPointSet, typename TOutputMesh>
void
QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>
::SetDestinationPoints( const InputPointSetType * destinationPointSet )
{
  itkDebugMacro("setting input ReferenceMesh to " << destinationPointSet);
  if( destinationPointSet != static_cast<const InputPointSetType *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<InputPointSetType *>( destinationPointSet ) );
    this->Modified();
    }
}

template <typename TInputMesh, typename TInputPointSet, typename TOutputMesh>
const typename
QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>::InputPointSetType
* QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>
::GetDestinationPoints() const
  {
  Self *                    surrogate = const_cast<Self *>( this );
  const InputPointSetType * destinationPointSet =
    static_cast<const InputPointSetType *>( surrogate->ProcessObject::GetInput(1) );
  return destinationPointSet;
  }

template <typename TInputMesh, typename TInputPointSet, typename TOutputMesh>
void
QuadEdgeMeshGenerateDeformationFieldFilter<TInputMesh, TInputPointSet, TOutputMesh>
::GenerateData()
{
  const InputMeshType * inputMesh = this->GetInputMesh();

  if( !inputMesh )
    {
    itkExceptionMacro("Input mesh is missing");
    }

  this->CopyInputMeshToOutputMesh();

  const InputPointSetType * inputPointSet = this->GetDestinationPoints();

  if( !inputPointSet )
    {
    itkExceptionMacro("Input PointSet is missing");
    }

  using DestinationPointsContainer = typename InputPointSetType::PointsContainer;

  const DestinationPointsContainer * destinationPoints = inputPointSet->GetPoints();

  if( !destinationPoints )
    {
    itkExceptionMacro("Input PointSet has no points");
    }

  const unsigned int numberOfPoints = inputMesh->GetNumberOfPoints();

  if( destinationPoints->Size() != numberOfPoints )
    {
    itkExceptionMacro("The PointSet does not have the same number of points as the Mesh");
    }

  OutputMeshType * outputMesh = this->GetOutput();

  ProgressReporter progress(this, 0, numberOfPoints);

  using DestinationPointIterator = typename DestinationPointsContainer::ConstIterator;

  DestinationPointIterator destinationPointItr = destinationPoints->Begin();
  DestinationPointIterator destinationPointEnd = destinationPoints->End();

  using SourcePointContainer = typename InputMeshType::PointsContainer;
  using SourcePointIterator = typename SourcePointContainer::ConstIterator;

  const SourcePointContainer * sourcePoints = inputMesh->GetPoints();

  SourcePointIterator sourcePointItr = sourcePoints->Begin();

  using DisplacementVectorContainer = typename OutputMeshType::PointDataContainer;
  using DisplacementVectorContainerPointer = typename DisplacementVectorContainer::Pointer;
  using DisplacementVectorIterator = typename DisplacementVectorContainer::Iterator;

  DisplacementVectorContainerPointer displacementVectors = DisplacementVectorContainer::New();

  displacementVectors->Reserve( numberOfPoints );

  DisplacementVectorIterator displacementVectorItr = displacementVectors->Begin();

  using VectorType = typename InputMeshType::PointType::VectorType;

  while( destinationPointItr != destinationPointEnd )
    {
    VectorType vectorToCenter = sourcePointItr.Value() - this->m_SphereCenter;

    VectorType dstVectorToCenter = destinationPointItr.Value() - this->m_SphereCenter;

    vectorToCenter.Normalize();

    displacementVectorItr.Value() =
      CrossProduct( vectorToCenter,
                    CrossProduct( vectorToCenter, dstVectorToCenter ) );

    displacementVectorItr.Value() *= -1.0;

    ++displacementVectorItr;
    ++sourcePointItr;
    ++destinationPointItr;
    }

  outputMesh->SetPointData( displacementVectors );
}
} // end namespace itk

#endif
