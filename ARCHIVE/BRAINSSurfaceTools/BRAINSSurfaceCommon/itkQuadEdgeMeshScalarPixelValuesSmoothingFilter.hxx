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
#ifndef __itkQuadEdgeMeshScalarPixelValuesSmoothingFilter_hxx
#define __itkQuadEdgeMeshScalarPixelValuesSmoothingFilter_hxx

#include "itkQuadEdgeMeshScalarPixelValuesSmoothingFilter.h"
#include "itkProgressReporter.h"
#include "itkVersor.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <typename TInputMesh, typename TOutputMesh>
QuadEdgeMeshScalarPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::QuadEdgeMeshScalarPixelValuesSmoothingFilter()
{
  this->m_Lambda = 1.0;
  this->m_MaximumNumberOfIterations = 10;
}

template <typename TInputMesh, typename TOutputMesh>
QuadEdgeMeshScalarPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::~QuadEdgeMeshScalarPixelValuesSmoothingFilter()
{
}

template <typename TInputMesh, typename TOutputMesh>
void
QuadEdgeMeshScalarPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::GenerateData()
{
  // Copy the input mesh into the output mesh.
  this->CopyInputMeshToOutputMesh();

  OutputMeshPointer outputMesh = this->GetOutput();

  //
  // Visit all nodes of the Mesh
  //
  OutputPointsContainerPointer points = outputMesh->GetPoints();

  if( points.IsNull() )
    {
    itkExceptionMacro("Mesh has NULL PointData");
    }

  OutputPointDataContainerPointer pointData = outputMesh->GetPointData();

  if( pointData.IsNull() )
    {
    itkExceptionMacro("Output Mesh has NULL PointData");
    }

  const double weightFactor = std::exp( -1.0 / ( 2.0 * this->m_Lambda ) );

  ProgressReporter progress(this, 0, this->m_MaximumNumberOfIterations);

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  std::cout << "Output Mesh numberOfPoints " << numberOfPoints << std::endl;
  for( unsigned int iter = 0; iter < this->m_MaximumNumberOfIterations; ++iter )
    {
    std::cout << " Smoothing Iteration " << iter << std::endl;

    using EdgeType = typename OutputMeshType::QEPrimal;

    OutputPointDataContainerPointer newPointDataContainer = OutputPointDataContainer::New();

    newPointDataContainer->Reserve( pointData->Size() );

    using AccumulatePixelType = typename NumericTraits<OutputPixelType>::AccumulateType;
    for( unsigned int pointId = 0; pointId < numberOfPoints; pointId++ )
      {
      const OutputPixelType & centralPixelValue = pointData->GetElement( pointId );

      const EdgeType * edgeToFirstNeighborPoint = outputMesh->FindEdge( pointId );
      const EdgeType * edgeToNeighborPoint = edgeToFirstNeighborPoint;

      AccumulatePixelType pixelSum = centralPixelValue;

      unsigned int numberOfNeighbors = 0;

      do
        {
        const OutputPointIdentifier neighborPointId = edgeToNeighborPoint->GetDestination();
        const OutputPixelType &     neighborPixelValue = pointData->GetElement( neighborPointId );

        pixelSum += weightFactor * neighborPixelValue;

        numberOfNeighbors++;

        edgeToNeighborPoint = edgeToNeighborPoint->GetOnext();
        }
      while( edgeToNeighborPoint != edgeToFirstNeighborPoint );

      const double normalizationFactor = 1.0 / ( 1.0 + numberOfNeighbors * weightFactor );

      OutputPixelType smoothedPixelValue = pixelSum * normalizationFactor;

      newPointDataContainer->SetElement( pointId, smoothedPixelValue );
      }

    outputMesh->SetPointData( newPointDataContainer );

    pointData = newPointDataContainer;

    progress.CompletedPixel();  // potential exception thrown here
    }
}
} // end namespace itk

#endif
