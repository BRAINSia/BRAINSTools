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
#ifndef __itkQuadEdgeMeshScalarPixelValuesSmoothingFilter_txx
#define __itkQuadEdgeMeshScalarPixelValuesSmoothingFilter_txx

#include "itkQuadEdgeMeshScalarPixelValuesSmoothingFilter.h"
#include "itkProgressReporter.h"
#include "itkVersor.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshScalarPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::QuadEdgeMeshScalarPixelValuesSmoothingFilter()
{
  this->m_Lambda = 1.0;
  this->m_MaximumNumberOfIterations = 10;
}

template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshScalarPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::~QuadEdgeMeshScalarPixelValuesSmoothingFilter()
{
}

template <class TInputMesh, class TOutputMesh>
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
  typedef typename OutputPointDataContainer::ConstIterator OutputPointDataIterator;

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

  const double weightFactor = vcl_exp( -1.0 / ( 2.0 * this->m_Lambda ) );

  ProgressReporter progress(this, 0, this->m_MaximumNumberOfIterations);

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  std::cout << "Output Mesh numberOfPoints " << numberOfPoints << std::endl;
  for( unsigned int iter = 0; iter < this->m_MaximumNumberOfIterations; ++iter )
    {
    std::cout << " Smoothing Iteration " << iter << std::endl;

    typedef typename OutputMeshType::QEPrimal EdgeType;

    OutputPointDataContainerPointer newPointDataContainer = OutputPointDataContainer::New();

    newPointDataContainer->Reserve( pointData->Size() );

    typedef typename NumericTraits<OutputPixelType>::AccumulateType AccumulatePixelType;
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
