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
#ifndef __itkQuadEdgeMeshVectorPixelValuesSmoothingFilter_txx
#define __itkQuadEdgeMeshVectorPixelValuesSmoothingFilter_txx

#include "itkQuadEdgeMeshVectorPixelValuesSmoothingFilter.h"
#include "itkProgressReporter.h"
#include "itkVersor.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshVectorPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::QuadEdgeMeshVectorPixelValuesSmoothingFilter()
{
  this->m_Lambda = 1.0;
  this->m_MaximumNumberOfIterations = 10;
  this->m_SphereCenter.Fill( 0.0 );
  this->m_SphereRadius = 1.0;
}

template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshVectorPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::~QuadEdgeMeshVectorPixelValuesSmoothingFilter()
{
}

template <class TInputMesh, class TOutputMesh>
void
QuadEdgeMeshVectorPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
::ParalelTransport(
  const OutputPointType sourcePoint, const OutputPointType destinationPoint,
  const InputPixelType & inputPixelValue,
  InputPixelType & transportedPixelValue ) const
{
  OutputVectorType vsrc = sourcePoint - this->m_SphereCenter;
  OutputVectorType vdst = destinationPoint - this->m_SphereCenter;

  OutputVectorType axis = CrossProduct( vsrc, vdst );

  const double scaledSinus   = axis.GetNorm();
  const double scaledCosinus = vsrc * vdst;

  double angle = vcl_atan2( scaledSinus, scaledCosinus );

  // if angle is very close to zero
  if( (angle - 0.0) < 0.000001 )
    {
    transportedPixelValue = inputPixelValue;
    }
  else
    {
    typedef Versor<double> VersorType;

    VersorType versor;
    versor.Set( axis, angle );

    transportedPixelValue = versor.Transform( inputPixelValue );
    }
}

template <class TInputMesh, class TOutputMesh>
void
QuadEdgeMeshVectorPixelValuesSmoothingFilter<TInputMesh, TOutputMesh>
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

  OutputPixelType smoothedPixelValue;

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  std::cout << "Output Mesh numberOfPoints " << numberOfPoints << std::endl;
  for( unsigned int iter = 0; iter < this->m_MaximumNumberOfIterations; ++iter )
    {
    std::cout << " Smoothing Iteration " << iter << std::endl;

    typedef typename OutputMeshType::QEPrimal EdgeType;

    OutputPointDataContainerPointer newPointDataContainer = OutputPointDataContainer::New();

    newPointDataContainer->Reserve( pointData->Size() );

    OutputPixelType transportedPixelValue;

    typedef typename NumericTraits<OutputPixelType>::AccumulateType AccumulatePixelType;
    for( unsigned int pointId = 0; pointId < numberOfPoints; pointId++ )
      {
      const OutputPointType & centralPoint = points->GetElement( pointId );
      const OutputPixelType & centralPixelValue = pointData->GetElement( pointId );

      const EdgeType * edgeToFirstNeighborPoint = outputMesh->FindEdge( pointId );
      const EdgeType * edgeToNeighborPoint = edgeToFirstNeighborPoint;

      AccumulatePixelType pixelSum;
      for( unsigned int k = 0; k < PointDimension; k++ )
        {
        pixelSum[k] = centralPixelValue[k];
        }

      unsigned int numberOfNeighbors = 0;

      do
        {
        const OutputPointIdentifier neighborPointId = edgeToNeighborPoint->GetDestination();
        const OutputPointType &     neighborPoint = points->GetElement( neighborPointId );
        const OutputPixelType &     neighborPixelValue = pointData->GetElement( neighborPointId );

        this->ParalelTransport( neighborPoint, centralPoint, neighborPixelValue, transportedPixelValue );
        for( unsigned int k = 0; k < PointDimension; k++ )
          {
          pixelSum[k] += weightFactor * transportedPixelValue[k];
          }

        numberOfNeighbors++;

        edgeToNeighborPoint = edgeToNeighborPoint->GetOnext();
        }
      while( edgeToNeighborPoint != edgeToFirstNeighborPoint );

      const double normalizationFactor = 1.0 / ( 1.0 + numberOfNeighbors * weightFactor );
      for( unsigned int k = 0; k < PointDimension; k++ )
        {
        smoothedPixelValue[k] = pixelSum[k] * normalizationFactor;
        }

      newPointDataContainer->SetElement( pointId, smoothedPixelValue );
      }

    outputMesh->SetPointData( newPointDataContainer );

    pointData = newPointDataContainer;

    progress.CompletedPixel();  // potential exception thrown here
    }
}
} // end namespace itk

#endif
