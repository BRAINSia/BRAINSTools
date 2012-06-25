/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2006-04-25 23:20:16 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter_hxx
#define __itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter_hxx

#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkProgressReporter.h"

namespace itk
{
template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::QuadEdgeMeshSphericalDiffeomorphicDemonsFilter()
{
  this->SetNumberOfIndexedInputs( 2 );
  this->SetNumberOfIndexedOutputs( 3 );

  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 3 );

  this->SetNthOutput( 0, OutputMeshType::New() );
  this->SetNthOutput( 1, FixedMeshType::New() );
  this->SetNthOutput( 2, DestinationPointSetType::New() );

  this->m_BasisSystemAtNode = BasisSystemContainerType::New();
  this->m_DestinationPoints = DestinationPointContainerType::New();
  this->m_DestinationPointsSwap = DestinationPointContainerType::New();
  this->m_UserProvidedInitialDestinationPoints = false;

  this->m_TriangleListBasisSystemCalculator = TriangleListBasisSystemCalculatorType::New();

  this->m_NodeScalarGradientCalculator = NodeScalarGradientCalculatorType::New();

  this->m_NodeVectorJacobianCalculator = NodeVectorJacobianCalculatorType::New();

  // this->m_ResampledMovingValuesContainer = ResampledMovingValuesContainerType::New();

  this->m_ScalarInterpolator = ScalarInterpolatorType::New();
  this->m_ScalarInterpolator->SetUseNearestNeighborInterpolationAsBackup(false);

  this->m_DeformationInterpolator = DeformationInterpolatorType::New();
  this->m_DeformationInterpolator->SetUseNearestNeighborInterpolationAsBackup(true);

  this->m_MaximumNumberOfIterations = 50;

  this->m_SphereCenter.Fill( 0.0 );
  this->m_SphereRadius = 1.0;

  this->m_SigmaX = 1.0;
  this->m_Epsilon = 1.0 / (this->m_SigmaX * this->m_SigmaX);

  this->m_Lambda = 1.0;
  this->m_MaximumNumberOfSmoothingIterations = 1;

  this->m_ShortestEdgeLength = 1.0;
  this->m_ScalingAndSquaringNumberOfIterations = 2;

  this->m_MetricValue = 0.0;
  this->m_MetricChange = 0.0;
  this->m_MetricSignificance = 1.0; // 1% change

  this->m_SelfRegulatedMode = false;

  this->m_SelfStopMode = false;

  this->m_FixedMeshAtInitialDestinationPoints = FixedMeshType::New();
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::~QuadEdgeMeshSphericalDiffeomorphicDemonsFilter()
{
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
DataObject::Pointer
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::MakeOutput(size_t idx)
{
  DataObject::Pointer output;

  switch( idx )
    {
    case 0:
      {
      output = (OutputMeshType::New() ).GetPointer();
      }
      break;
    case 1:
      {
      output = (FixedMeshType::New() ).GetPointer();
      }
      break;
    case 2:
      {
      output = (DestinationPointSetType::New() ).GetPointer();
      }
      break;
    }
  return output.GetPointer();
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::SetFixedMesh( const FixedMeshType * fixedMesh )
{
  itkDebugMacro("setting Fixed Mesh to " << fixedMesh );

  if( this->m_FixedMesh.GetPointer() != fixedMesh )
    {
    this->m_FixedMesh = fixedMesh;

    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(0, const_cast<FixedMeshType *>( fixedMesh ) );

    this->Modified();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::SetMovingMesh( const MovingMeshType * movingMesh )
{
  itkDebugMacro("setting Moving Mesh to " << movingMesh );

  if( this->m_MovingMesh.GetPointer() != movingMesh )
    {
    this->m_MovingMesh = movingMesh;

    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(1, const_cast<MovingMeshType *>( movingMesh ) );

    this->Modified();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::SetInitialDestinationPoints( const DestinationPointSetType * destinationPointSet )
{
  itkDebugMacro("setting Destination PointSet to " << destinationPointSet );

  if( this->GetInitialDestinationPoints() != destinationPointSet )
    {
    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(2, const_cast<DestinationPointSetType *>( destinationPointSet ) );

    this->m_UserProvidedInitialDestinationPoints = true;

    this->Modified();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
const typename QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                                              TOutputMesh>::DestinationPointSetType
* QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::GetInitialDestinationPoints() const
  {
  if( this->GetNumberOfInputs() < 3 )
    {
    return NULL;
    }

  return static_cast<const DestinationPointSetType *>(this->ProcessObject::GetInput(2) );
  }

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
typename QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::DestinationPointSetType
* QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::
GetFinalDestinationPoints() const
  {
  if( this->GetNumberOfOutputs() < 3 )
    {
    return 0;
    }

  const DestinationPointSetType * pointSet =
    dynamic_cast<const DestinationPointSetType *>(this->ProcessObject::GetOutput(2) );

  return const_cast<DestinationPointSetType *>( pointSet );
  }

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::CopyInitialDestinationPoints()
{
  this->m_FixedMeshAtInitialDestinationPoints = FixedMeshType::New();

  this->CopyMeshToMesh( this->m_FixedMesh, this->m_FixedMeshAtInitialDestinationPoints  );

  const DestinationPointSetType *       destinationPointSet = this->GetInitialDestinationPoints();
  const DestinationPointContainerType * destinationPoints = destinationPointSet->GetPoints();

  itkDebugMacro("setting Destination Points to " << destinationPoints );

  if( destinationPoints == NULL )
    {
    itkExceptionMacro("Pointer to DestinationPoints was NULL");
    }

  if( destinationPoints->Size() != this->m_FixedMesh->GetNumberOfPoints() )
    {
    itkExceptionMacro(
      "Number of destination points " << destinationPoints->Size()
                                      << " does not match fixed mesh number of points "
                                      << this->m_FixedMesh->GetNumberOfPoints() );
    }

  FixedPointsContainer * fixedPoints = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();
  FixedPointsIterator    fixedPointItr = fixedPoints->Begin();

  this->m_DestinationPoints = DestinationPointContainerType::New();
  this->m_DestinationPoints->Reserve( destinationPoints->Size() );

  DestinationPointConstIterator srcPointItr = destinationPoints->Begin();

  DestinationPointIterator dstPointItr = this->m_DestinationPoints->Begin();
  DestinationPointIterator dstPointEnd = this->m_DestinationPoints->End();

  PointType point;
  while( dstPointItr != dstPointEnd )
    {
    point = srcPointItr.Value();

    this->ProjectPointToSphereSurface( point );

    dstPointItr.Value() = point;

    fixedPointItr.Value().SetPoint( point );

    ++srcPointItr;
    ++dstPointItr;
    ++fixedPointItr;
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::GenerateData()
{
  this->m_Chronometer.Start("DataPreProcessing");
  // Prepare data
  this->CopyInputMeshToOutputMesh();
  this->AllocateInternalArrays();
  this->ComputeInitialArrayOfDestinationPoints();
  this->InitializeFixedNodesSigmas();
  this->ComputeBasisSystemAtEveryNode();
  this->ComputeShortestEdgeLength();
  this->ComposeDestinationPointsOutputPointSet();
  this->InitializeInterpolators();
  this->InitializeGradientCalculators();
  this->m_Chronometer.Stop("DataPreProcessing");

  // Compute deformations
  this->m_Chronometer.Start("RunIterations");
  this->RunIterations();
  this->m_Chronometer.Stop("RunIterations");

  // Gathering outputs
  this->m_Chronometer.Start("DataPostProcessing");
  this->ComputeMappedMovingValueAtEveryNode();
  this->AssignResampledMovingValuesToOutputMesh();
  this->ComposeFixedMeshOutputDisplacedToMovingMesh();
  this->m_Chronometer.Stop("DataPostProcessing");
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::ChronometerReport(
  std::ofstream & os ) const
{
  this->m_Chronometer.Report( os );
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::AllocateInternalArrays()
{
  const PointIdentifier numberOfNodes = this->m_FixedMesh->GetNumberOfPoints();

  //
  // create new containers and allocate memory for them, in case the filter has
  // been run previously with a mesh having a larger number of nodes than the
  // current mesh.
  //

  this->m_BasisSystemAtNode = BasisSystemContainerType::New();
  this->m_BasisSystemAtNode->Reserve( numberOfNodes );

  if( !this->m_UserProvidedInitialDestinationPoints )
    {
    this->m_DestinationPoints = DestinationPointContainerType::New();
    this->m_DestinationPoints->Reserve( numberOfNodes );
    }

  this->m_DestinationPointsSwap = DestinationPointContainerType::New();
  this->m_DestinationPointsSwap->Reserve( numberOfNodes );

  this->m_DisplacementField = DestinationPointContainerType::New();
  this->m_DisplacementField->Reserve( numberOfNodes );

  this->m_DisplacementFieldSwap = DestinationPointContainerType::New();
  this->m_DisplacementFieldSwap->Reserve( numberOfNodes );

  this->m_ResampledMovingValuesContainer = ResampledMovingValuesContainerType::New();
  this->m_ResampledMovingValuesContainer->Reserve( numberOfNodes );

  this->m_VelocityField = VelocityVectorContainer::New();
  this->m_VelocityField->Reserve( numberOfNodes );

  this->m_TangentVectorField = TangentVectorContainer::New();
  this->m_TangentVectorField->Reserve( numberOfNodes );

  this->m_TangentVectorFieldSwap = TangentVectorContainer::New();
  this->m_TangentVectorFieldSwap->Reserve( numberOfNodes );

  this->m_ShortestEdgeLengthPerPoint = ShortestLengthContainerType::New();
  this->m_ShortestEdgeLengthPerPoint->Reserve( numberOfNodes );
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::InitializeFixedNodesSigmas()
{
  const PointIdentifier numberOfNodes =
    this->m_FixedMeshAtInitialDestinationPoints->GetNumberOfPoints();

  if( this->m_FixedNodesSigmas.IsNull() || this->m_FixedNodesSigmas->Size() != numberOfNodes )
    {
    NodeSigmaContainerPointer sigmas = NodeSigmaContainerType::New();
    sigmas->Reserve( numberOfNodes );
    NodeSigmaContainerIterator sigmaItr = sigmas->Begin();
    NodeSigmaContainerIterator sigmaEnd = sigmas->End();
    while( sigmaItr != sigmaEnd )
      {
      sigmaItr.Value() = 1.0;
      ++sigmaItr;
      }

    this->SetFixedNodesSigmas( sigmas );
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::ComputeBasisSystemAtEveryNode()
{
  const PointIdentifier numberOfNodes = this->m_FixedMeshAtInitialDestinationPoints->GetNumberOfPoints();

  typedef typename TFixedMesh::PointsContainer PointsContainer;
  const PointsContainer * points = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();

  typedef typename TFixedMesh::QEPrimal EdgeType;
  for( PointIdentifier pointId1 = 0; pointId1 < numberOfNodes; pointId1++ )
    {
    const EdgeType * edge = this->m_FixedMeshAtInitialDestinationPoints->FindEdge( pointId1 );

    if( !edge )
      {
      itkExceptionMacro("FindEdge() returned NULL for pointId " << pointId1 );
      }

    PointIdentifier pointId2 = edge->GetDestination();

    const PointType point1 = points->GetElement( pointId1 );
    const PointType point2 = points->GetElement( pointId2 );

    const VectorType v12    = point1 - point2;

    // v12 is not necessarily tangent to the sphere, therefore we must use
    // cross products in order to find an orthogonal system.

    const VectorType radial = point1.GetVectorFromOrigin();

    VectorType u12 = CrossProduct( v12, radial );
    VectorType w12 = CrossProduct( radial, u12 );

    w12.Normalize();
    u12.Normalize();

    BasisSystemType basis;
    basis.SetVector( 0, w12 );
    basis.SetVector( 1, u12 );

    this->m_BasisSystemAtNode->SetElement( pointId1, basis );
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeInitialArrayOfDestinationPoints()
{
  if( this->m_UserProvidedInitialDestinationPoints )
    {
    this->CopyInitialDestinationPoints();
    }
  else
    {
    this->CopySourcePoinstAsDestinationPoints();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::CopySourcePoinstAsDestinationPoints()
{
  this->m_FixedMeshAtInitialDestinationPoints = FixedMeshType::New();

  this->CopyMeshToMesh( this->m_FixedMesh, this->m_FixedMeshAtInitialDestinationPoints  );

  const FixedPointsContainer * points = this->m_FixedMesh->GetPoints();

  FixedPointsConstIterator srcPointItr = points->Begin();

  DestinationPointIterator dstPointItr = this->m_DestinationPoints->Begin();
  DestinationPointIterator dstPointEnd = this->m_DestinationPoints->End();

  while( dstPointItr != dstPointEnd )
    {
    dstPointItr.Value() = srcPointItr.Value();
    ++srcPointItr;
    ++dstPointItr;
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::InitializeInterpolators()
{
  this->m_ScalarInterpolator->SetInputMesh( this->m_MovingMesh );
  this->m_ScalarInterpolator->Initialize();

  this->m_DeformationInterpolator->SetInputMesh( this->m_FixedMeshAtInitialDestinationPoints );
  this->m_DeformationInterpolator->Initialize();
  this->m_DeformationInterpolator->SetSphereCenter( this->m_SphereCenter );
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::InitializeGradientCalculators()
{
  this->m_TriangleListBasisSystemCalculator->SetInputMesh( this->m_FixedMeshAtInitialDestinationPoints );
  this->m_TriangleListBasisSystemCalculator->Calculate();

  this->m_NodeScalarGradientCalculator->SetInputMesh( this->m_FixedMeshAtInitialDestinationPoints );
  this->m_NodeScalarGradientCalculator->SetDataContainer( this->m_ResampledMovingValuesContainer );

  this->m_NodeScalarGradientCalculator->SetBasisSystemList(
    this->m_TriangleListBasisSystemCalculator->GetBasisSystemList() );

  this->m_NodeScalarGradientCalculator->SetSphereCenter( this->m_SphereCenter );
  this->m_NodeScalarGradientCalculator->SetSphereRadius( this->m_SphereRadius );

  this->m_NodeVectorJacobianCalculator->SetInputMesh( this->m_FixedMeshAtInitialDestinationPoints );
  this->m_NodeVectorJacobianCalculator->SetBasisSystemList(
    this->m_TriangleListBasisSystemCalculator->GetBasisSystemList() );

  this->m_NodeVectorJacobianCalculator->SetSphereCenter( this->m_SphereCenter );
  this->m_NodeVectorJacobianCalculator->SetSphereRadius( this->m_SphereRadius );
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::RunIterations()
{
  // Report the progress
  ProgressReporter progress( this, 0, this->m_MaximumNumberOfIterations );

  double pre_Metric = 0.0; // to save the metric value of previous iteration

  unsigned int count = 0;
  bool         timeOn = false;

  for( unsigned int i = 0; i < this->m_MaximumNumberOfIterations; i++ )
    {
    this->m_Chronometer.Start("ComputeMappedMovingValueAtEveryNode");
    this->ComputeMappedMovingValueAtEveryNode();
    this->m_Chronometer.Stop("ComputeMappedMovingValueAtEveryNode");

    this->m_Chronometer.Start("ComputeGradientsOfMappedMovingValueAtEveryNode");
    this->ComputeGradientsOfMappedMovingValueAtEveryNode();
    this->m_Chronometer.Stop("ComputeGradientsOfMappedMovingValueAtEveryNode");

    std::cout << "Iteration: " << i << std::endl;
    this->m_Chronometer.Start("ComputeSelfRegulatedVelocityField");
    this->ComputeSelfRegulatedVelocityField();
    this->m_Chronometer.Stop("ComputeSelfRegulatedVelocityField");

    // metric calculation
    // iteration stops when there is no significant change in continuous 5 iterations
    if( m_SelfStopMode )
      {
      if( (i > 15) && (pre_Metric != 0.0) )
        {
        m_MetricChange = fabs( this->GetMetricValue() - pre_Metric )
          / pre_Metric * 100.0;

        if( m_MetricChange < m_MetricSignificance )
          {
          count += 1;
          timeOn = true;
          }
        else
          {
          timeOn = false;
          count = 0;
          }
        if( count == 5 )
          {
          break;
          }
        }
      pre_Metric = this->GetMetricValue();
      }

    this->m_Chronometer.Start("ComputeScalingAndSquaringNumberOfIterations");
    this->ComputeScalingAndSquaringNumberOfIterations();
    this->m_Chronometer.Stop("ComputeScalingAndSquaringNumberOfIterations");

    this->m_Chronometer.Start("ComputeDeformationByScalingAndSquaring");
    this->ComputeDeformationByScalingAndSquaring();
    this->m_Chronometer.Stop("ComputeDeformationByScalingAndSquaring");

    this->m_Chronometer.Start("ComposeDeformationUpdateWithPreviousDeformation");
    this->ComposeDeformationUpdateWithPreviousDeformation();
    this->m_Chronometer.Stop("ComposeDeformationUpdateWithPreviousDeformation");

    this->m_Chronometer.Start("SmoothDeformationField");
    this->SmoothDeformationField();
    this->m_Chronometer.Stop("SmoothDeformationField");

    this->m_Chronometer.Start("ComposeDestinationPointsOutputPointSet");
    this->ComposeDestinationPointsOutputPointSet();
    this->m_Chronometer.Stop("ComposeDestinationPointsOutputPointSet");

    this->m_Chronometer.Start("ComputeSelfRegulatedSigmaXandEpsilon");
    this->ComputeSelfRegulatedSigmaXandEpsilon();
    this->m_Chronometer.Stop("ComputeSelfRegulatedSigmaXandEpsilon");

    // Report progress via Events
    progress.CompletedPixel();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeGradientsOfMappedMovingValueAtEveryNode()
{
  this->m_NodeScalarGradientCalculator->Initialize();
  this->m_NodeScalarGradientCalculator->Compute();

  this->m_NodeVectorJacobianCalculator->SetVectorContainer( this->m_DestinationPoints );

  this->m_NodeVectorJacobianCalculator->Initialize();
  this->m_NodeVectorJacobianCalculator->Compute();
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeMappedMovingValueAtEveryNode()
{
  DestinationPointIterator pointItr = this->m_DestinationPoints->Begin();
  DestinationPointIterator pointEnd = this->m_DestinationPoints->End();

  ResampledMovingValuesContainerIterator resampledArrayItr = this->m_ResampledMovingValuesContainer->Begin();

  while( pointItr != pointEnd )
    {
    resampledArrayItr.Value() = this->m_ScalarInterpolator->Evaluate( pointItr.Value() );

    ++pointItr;
    ++resampledArrayItr;
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeSelfRegulatedVelocityField()
{
  if( this->m_SelfRegulatedMode )
    {
    unsigned int iterations = 0;

    do
      {
      this->m_Chronometer.Start("ComputeVelocityField");
      this->ComputeVelocityField();
      this->m_Chronometer.Stop("ComputeVelocityField");
      this->ComputeSelfRegulatedSigmaXandEpsilon();
      // std::cout<<this->ComputeLargestVelocityMagnitude()<<std::endl;
      // std::cout<<this->m_Epsilon<<std::endl;
      // std::cout<<this->m_SigmaX<<std::endl;
      }
    while( ( this->m_LargestVelocityToEdgeLengthRatio > 1.5 ) && ( iterations++ < 10 ) );
    }
  else
    {
    this->ComputeVelocityField();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::ComputeVelocityField()
{
  const PointIdentifier numberOfNodes = this->m_FixedMeshAtInitialDestinationPoints->GetNumberOfPoints();

  const FixedPointsContainer * points = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();

  FixedPointsConstIterator pointItr = points->Begin();

  const FixedPointDataContainer * pointData = this->m_FixedMeshAtInitialDestinationPoints->GetPointData();

  FixedPointDataConstIterator fixedPointDataItr = pointData->Begin();

  BasisSystemContainerIterator basisItr = this->m_BasisSystemAtNode->Begin();

  NodeSigmaContainerConstIterator sigmaItr = this->m_FixedNodesSigmas->Begin();

  VelocityVectorIterator velocityItr = this->m_VelocityField->Begin();

  ResampledMovingValuesContainerIterator resampledArrayItr =
    this->m_ResampledMovingValuesContainer->Begin();

  typedef vnl_matrix_fixed<double, 3, 3> VnlMatrix33Type;
  typedef vnl_vector_fixed<double, 2>    VnlVector2Type;
  typedef vnl_vector_fixed<double, 3>    VnlVector3Type;
  typedef vnl_matrix_fixed<double, 3, 2> VnlMatrix32Type;
  typedef vnl_matrix_fixed<double, 2, 3> VnlMatrix23Type;
  typedef vnl_matrix_fixed<double, 2, 2> VnlMatrix22Type;

  VnlMatrix33Type Gn;
  VnlMatrix33Type Gn2;
  VnlMatrix33Type mn2;
  VnlMatrix32Type En;
  VnlMatrix23Type EnT;
  VnlMatrix22Type EpsilonI22;
  VnlMatrix33Type Gn2Sn2;
  VnlMatrix33Type Gn2Sn2m2;
  VnlMatrix22Type EnTGn2Sn2m2En;
  VnlMatrix22Type EnTGn2Sn2m2EnGI22;
  VnlMatrix22Type EnTGn2Sn2m2EnGI22I;
  VnlMatrix33Type SnT;
  VnlMatrix33Type Gn2Sn;

  VnlVector3Type     mn;
  VnlVector2Type     EnTmn;
  VnlVector3Type     intensitySlope;
  VelocityVectorType Vn;

  VectorType vectorToCenter;

  typedef typename NodeVectorJacobianCalculatorType::OutputType JacobianType;

  JacobianType destinationJacobian;

  EpsilonI22.set_identity();
  EpsilonI22 *= ( this->m_Epsilon );

  double sumOfSquaredDifferences = 0.0;

  const double sigmaX2 = ( this->m_SigmaX * this->m_SigmaX );
  for( PointIdentifier pointId = 0; pointId < numberOfNodes; pointId++ )
    {
    vectorToCenter = pointItr.Value() - this->m_SphereCenter;

    vectorToCenter.Normalize();

    Gn(0, 0) = 0.0;
    Gn(1, 1) = 0.0;
    Gn(2, 2) = 0.0;

    Gn(0, 1) = -vectorToCenter[2];
    Gn(0, 2) =  vectorToCenter[1];
    Gn(1, 2) = -vectorToCenter[0];

    Gn(1, 0) =  vectorToCenter[2];
    Gn(2, 0) = -vectorToCenter[1];
    Gn(2, 1) =  vectorToCenter[0];

    Gn2 = Gn * Gn;

    typedef typename NodeScalarGradientCalculatorType::DerivativeType DerivativeType;
    DerivativeType derivative = this->m_NodeScalarGradientCalculator->Evaluate( pointId );

    destinationJacobian = this->m_NodeVectorJacobianCalculator->Evaluate( pointId );

    const BasisSystemType &   basis = basisItr.Value();
    const VectorType &        v0 = basis.GetVector(0);
    const VectorType &        v1 = basis.GetVector(1);
    const MovingPixelRealType Mv = resampledArrayItr.Value();
    const FixedPixelRealType  Fv = fixedPointDataItr.Value();
    for( unsigned int i = 0; i < 3; i++ )
      {
      En(i, 0) = v0[i];
      En(i, 1) = v1[i];
      EnT(0, i) = v0[i];
      EnT(1, i) = v1[i];
      mn[i] = derivative[i];
      }
    for( unsigned int r = 0; r < 3; r++ )
      {
      for( unsigned int c = 0; c < 3; c++ )
        {
        mn2(r, c) = mn[r] * mn[c];
        SnT(r, c) = destinationJacobian(c, r);  // FIXME : Check for potential transposition here...
        }
      }

    Gn2Sn = Gn2 * SnT;

    Gn2Sn2 = Gn2Sn.transpose() * Gn2Sn;

    //
    // The general form of this addition would involve two weights,
    // representing the variance of each term at this node.
    //
    const double sigmaN2 = sigmaItr.Value() * sigmaItr.Value();

    Gn2Sn2m2 = mn2 / sigmaN2 + Gn2Sn2 / sigmaX2;

    EnTGn2Sn2m2En = EnT * Gn2Sn2m2 * En;

    EnTGn2Sn2m2EnGI22 = EnTGn2Sn2m2En + EpsilonI22;

    EnTGn2Sn2m2EnGI22I = vnl_matrix_inverse<double>( EnTGn2Sn2m2EnGI22 );

    EnTmn = EnT * mn;

    intensitySlope = En * EnTGn2Sn2m2EnGI22I * EnTmn;

    Vn.SetVnlVector( intensitySlope * ( Fv - Mv ) );

    sumOfSquaredDifferences += ( Fv - Mv ) * ( Fv - Mv ) / sigmaN2;

    velocityItr.Value() = Vn;

    ++velocityItr;
    ++sigmaItr;
    ++basisItr;
    ++resampledArrayItr;
    ++fixedPointDataItr;
    ++pointItr;
    }

  const double averageOfSquaredDifferences = sumOfSquaredDifferences / numberOfNodes;

  this->m_MetricValue = averageOfSquaredDifferences;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeScalingAndSquaringNumberOfIterations()
{
  const unsigned int minimumNumberOfIterations = 1;

  if( this->m_SelfRegulatedMode )
    {
    if( this->m_LargestVelocityToEdgeLengthRatio < 1.0 )
      {
      this->m_ScalingAndSquaringNumberOfIterations = minimumNumberOfIterations;
      }
    else
      {
      unsigned int iterations =
        static_cast<unsigned int>(
          vcl_log( this->m_LargestVelocityToEdgeLengthRatio ) / vcl_log( 2.0 ) ) + 2;

      if( iterations < minimumNumberOfIterations )
        {
        iterations = minimumNumberOfIterations;
        }

      this->m_ScalingAndSquaringNumberOfIterations = iterations;
      }
    }
  else
    {
    const double ratio = this->ComputeLargestVelocityMagnitude() / this->m_ShortestEdgeLength;

    if( ratio < 1.0 )
      {
      this->m_ScalingAndSquaringNumberOfIterations = minimumNumberOfIterations;
      }
    else
      {
      unsigned int iterations =
        static_cast<unsigned int>( vcl_log( ratio ) / vcl_log( 2.0 ) ) + 2;

      if( iterations < minimumNumberOfIterations )
        {
        iterations = minimumNumberOfIterations;
        }

      this->m_ScalingAndSquaringNumberOfIterations = iterations;
      }
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::ComputeShortestEdgeLength()
{
  double shortestLength = NumericTraits<double>::max();

  typedef typename FixedMeshType::QEPrimal EdgeType;

  const FixedPointsContainer * points = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();

  FixedPointsConstIterator pointItr = points->Begin();
  FixedPointsConstIterator pointEnd = points->End();

  ShortestLengthContainerIterator shortestEdgeItr = this->m_ShortestEdgeLengthPerPoint->Begin();

  while( pointItr != pointEnd )
    {
    EdgeType * edge1 = this->m_FixedMeshAtInitialDestinationPoints->FindEdge( pointItr.Index() );

    EdgeType * temp1 = NULL;
    EdgeType * temp2 = edge1;

    const PointType & point = pointItr.Value();

    double localShortestLength = NumericTraits<double>::max();

    do
      {
      temp1 = temp2;
      temp2 = temp1->GetOnext();

      const PointIdentifier neighborPointId = temp1->GetDestination();

      const PointType & neighborPoint = points->GetElement( neighborPointId );

      const double distance = point.EuclideanDistanceTo( neighborPoint );

      if( distance < localShortestLength )
        {
        localShortestLength = distance;
        }
      }
    while( temp2 != edge1 );

    shortestEdgeItr.Value() = localShortestLength;

    if( localShortestLength < shortestLength )
      {
      shortestLength = localShortestLength;
      }

    ++pointItr;
    ++shortestEdgeItr;
    }

  this->m_ShortestEdgeLength = shortestLength;
  // std::cout << "m_ShortestEdgeLength = " << this->m_ShortestEdgeLength << std::endl;

  if( this->m_ShortestEdgeLength < vnl_math::eps )
    {
    itkExceptionMacro("The shortest edge length is too close to zero = " << shortestLength );
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
double
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeLargestVelocityMagnitude() const
{
  double largestVelocityMagnitude = NumericTraits<double>::Zero;

  VelocityVectorConstIterator velocityItr = this->m_VelocityField->Begin();
  VelocityVectorConstIterator velocityEnd = this->m_VelocityField->End();

  while( velocityItr != velocityEnd )
    {
    const double velocityMagnitude = velocityItr.Value().GetNorm();

    if( velocityMagnitude > largestVelocityMagnitude )
      {
      largestVelocityMagnitude = velocityMagnitude;
      }

    ++velocityItr;
    }

  return largestVelocityMagnitude;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeDeformationByScalingAndSquaring()
{
  unsigned long powerOfTwo = 1;

  powerOfTwo <<= this->m_ScalingAndSquaringNumberOfIterations;

  const double scalingFactor = 1.0 / powerOfTwo;

  DestinationPointIterator displacementItr = this->m_DisplacementField->Begin();
  DestinationPointIterator displacementEnd = this->m_DisplacementField->End();

  VelocityVectorConstIterator velocityItr = this->m_VelocityField->Begin();

  const FixedPointsContainer * points = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();

  FixedPointsConstIterator pointItr = points->Begin();

  PointType destinationPoint;

  while( displacementItr != displacementEnd )
    {
    destinationPoint = pointItr.Value() +  velocityItr.Value() * scalingFactor;

    this->ProjectPointToSphereSurface( destinationPoint );

    displacementItr.Value() = destinationPoint;

    ++displacementItr;
    ++pointItr;
    ++velocityItr;
    }
  for( unsigned int i = 0; i < this->m_ScalingAndSquaringNumberOfIterations; i++ )
    {
    DestinationPointConstIterator oldDisplacementItr = this->m_DisplacementField->Begin();
    DestinationPointConstIterator oldDisplacementEnd = this->m_DisplacementField->End();

    DestinationPointIterator newDisplacementItr = this->m_DisplacementFieldSwap->Begin();

    while( oldDisplacementItr != oldDisplacementEnd )
      {
      destinationPoint = oldDisplacementItr.Value();

      this->ProjectPointToSphereSurface( destinationPoint );

      newDisplacementItr.Value() =
        this->InterpolateDestinationFieldAtPoint(
          this->m_DisplacementField, destinationPoint );

      ++newDisplacementItr;
      ++oldDisplacementItr;
      }

    this->SwapOldAndNewDisplacementFieldContainers();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComputeSelfRegulatedSigmaXandEpsilon()
{
  if( this->m_SelfRegulatedMode )
    {
    //
    // Largest velocity vector Vn  magnitude  / 2^(N-2) < 1/2 Vertex distance
    //
    //  const double largestVelocityMagnitude = this->ComputeLargestVelocityMagnitude();
    //

    this->ComputeLargestVelocityMagnitudeToShortestEdgeLengthRatio();

    this->m_SigmaX /= vcl_sqrt( this->m_LargestVelocityToEdgeLengthRatio );
    this->m_Epsilon =  1.0 / ( this->m_SigmaX * this->m_SigmaX );
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComposeDeformationUpdateWithPreviousDeformation()
{
  DestinationPointConstIterator displacementItr = this->m_DisplacementField->Begin();
  DestinationPointConstIterator displacementEnd = this->m_DisplacementField->End();

  DestinationPointIterator newDestinationPointItr = this->m_DestinationPointsSwap->Begin();

  while( displacementItr != displacementEnd )
    {
    PointType point = displacementItr.Value();

    this->ProjectPointToSphereSurface( point );

    PointType destinationPoint =
      this->InterpolateDestinationFieldAtPoint( this->m_DestinationPoints, point );

    this->ProjectPointToSphereSurface( destinationPoint );

    newDestinationPointItr.Value() = destinationPoint;

    ++newDestinationPointItr;
    ++displacementItr;
    }

  this->SwapOldAndNewDestinationPointContainers();
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
typename QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::PointType
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::InterpolateDestinationFieldAtPoint(
  const DestinationPointContainerType * destinationField,
  const PointType & point )
{
  PointType interpolatedDestinationPoint;

  const bool found = this->m_DeformationInterpolator->Evaluate( destinationField, point, interpolatedDestinationPoint );

  if( !found )
    {
    itkExceptionMacro("Point not found in the interpolation" << point );
    }

  this->ProjectPointToSphereSurface( interpolatedDestinationPoint );

  return interpolatedDestinationPoint;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::ProjectPointToSphereSurface(
  PointType & point ) const
{
  VectorType vectorToCenter = point - this->m_SphereCenter;

  const double radialDistance = vectorToCenter.GetNorm();

  vectorToCenter *= this->m_SphereRadius / radialDistance;
  point = this->m_SphereCenter + vectorToCenter;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::SwapOldAndNewDisplacementFieldContainers()
{
  DestinationPointContainerPointer temp = this->m_DisplacementField;

  this->m_DisplacementField = this->m_DisplacementFieldSwap;
  this->m_DisplacementFieldSwap = temp;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::SwapOldAndNewDestinationPointContainers()
{
  DestinationPointContainerPointer temp = this->m_DestinationPoints;

  this->m_DestinationPoints = this->m_DestinationPointsSwap;
  this->m_DestinationPointsSwap = temp;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::SwapOldAndNewTangetFieldContainers()
{
  TangentVectorPointer temp = this->m_TangentVectorField;

  this->m_TangentVectorField = this->m_TangentVectorFieldSwap;
  this->m_TangentVectorFieldSwap = temp;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::SmoothDeformationField()
{
  this->ConvertDeformationFieldToTangentVectorField();
  this->SmoothTangentVectorField();
  this->ConvertTangentVectorFieldToDeformationField();
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ConvertDeformationFieldToTangentVectorField()
{
  DestinationPointIterator dstPointItr = this->m_DestinationPoints->Begin();
  DestinationPointIterator dstPointEnd = this->m_DestinationPoints->End();

  const FixedPointsContainer * points = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();

  FixedPointsConstIterator pointItr = points->Begin();

  TangentVectorIterator tangentItr = this->m_TangentVectorField->Begin();

  const double factor = -1.0 / this->m_SphereRadius;

  while( dstPointItr != dstPointEnd )
    {
    VectorType vectorToCenter = pointItr.Value() - this->m_SphereCenter;

    vectorToCenter.Normalize();

    tangentItr.Value() =
      CrossProduct( vectorToCenter,
                    CrossProduct( vectorToCenter, dstPointItr.Value().GetVectorFromOrigin() ) );

    tangentItr.Value() *= factor;

    ++dstPointItr;
    ++tangentItr;
    ++pointItr;
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::SmoothTangentVectorField()
{
  const double weightFactor = vcl_exp( -1.0 / ( 2.0 * this->m_Lambda ) );

  const FixedPointsContainer * points = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();

  TangentVectorType smoothedVector;
  TangentVectorType transportedTangentVector;

  for( unsigned int iter = 0; iter < this->m_MaximumNumberOfSmoothingIterations; ++iter )
    {
    typedef typename OutputMeshType::QEPrimal EdgeType;

    TangentVectorIterator tangentItr = this->m_TangentVectorField->Begin();
    TangentVectorIterator tangentEnd = this->m_TangentVectorField->End();

    TangentVectorIterator smoothedTangentItr = this->m_TangentVectorFieldSwap->Begin();

    FixedPointsConstIterator pointItr = points->Begin();

    typedef typename NumericTraits<TangentVectorType>::AccumulateType AccumulatePixelType;

    while( tangentItr != tangentEnd )
      {
      const TangentVectorType & centralTangentVector = tangentItr.Value();

      const EdgeType * edgeToFirstNeighborPoint =
        this->m_FixedMeshAtInitialDestinationPoints->FindEdge( tangentItr.Index() );
      const EdgeType * edgeToNeighborPoint = edgeToFirstNeighborPoint;

      AccumulatePixelType tangentVectorSum;
      for( unsigned int k = 0; k < PointDimension; k++ )
        {
        tangentVectorSum[k] = centralTangentVector[k];
        }

      unsigned int numberOfNeighbors = 0;

      do
        {
        const PointIdentifier     neighborPointId = edgeToNeighborPoint->GetDestination();
        const PointType &         neighborPoint = points->GetElement( neighborPointId );
        const TangentVectorType & neighborTangentVector =
          this->m_TangentVectorField->GetElement( neighborPointId );

        this->ParalelTransport( neighborPoint, pointItr.Value(),
                                neighborTangentVector, transportedTangentVector );
        for( unsigned int k = 0; k < PointDimension; k++ )
          {
          tangentVectorSum[k] += weightFactor * transportedTangentVector[k];
          }

        numberOfNeighbors++;

        edgeToNeighborPoint = edgeToNeighborPoint->GetOnext();
        }
      while( edgeToNeighborPoint != edgeToFirstNeighborPoint );

      const double normalizationFactor = 1.0 / ( 1.0 + numberOfNeighbors * weightFactor );
      for( unsigned int k = 0; k < PointDimension; k++ )
        {
        smoothedVector[k] = tangentVectorSum[k] * normalizationFactor;
        }

      smoothedTangentItr.Value() = smoothedVector;

      ++tangentItr;
      ++smoothedTangentItr;
      ++pointItr;
      }

    this->SwapOldAndNewTangetFieldContainers();
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::ParalelTransport(
  const PointType sourcePoint, const PointType destinationPoint,
  const TangentVectorType & inputVector,
  TangentVectorType & transportedVector ) const
{
  VectorType vsrc = sourcePoint - this->m_SphereCenter;
  VectorType vdst = destinationPoint - this->m_SphereCenter;

  VectorType axis = CrossProduct( vsrc, vdst );

  const double scaledSinus   = axis.GetNorm();
  const double scaledCosinus = vsrc * vdst;

  double angle = vcl_atan2( scaledSinus, scaledCosinus );

  typedef Versor<double> VersorType;

  VersorType versor;
  versor.Set( axis, angle );

  transportedVector = versor.Transform( inputVector );
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ConvertTangentVectorFieldToDeformationField()
{
  TangentVectorIterator tangentItr = this->m_TangentVectorField->Begin();
  TangentVectorIterator tangentEnd = this->m_TangentVectorField->End();

  DestinationPointIterator dstPointItr = this->m_DestinationPoints->Begin();

  const FixedPointsContainer * points = this->m_FixedMeshAtInitialDestinationPoints->GetPoints();

  FixedPointsConstIterator pointItr = points->Begin();

  typedef Versor<double> VersorType;
  VersorType versor;

  const double normEpsilon = itk::NumericTraits<double>::min();

  while( tangentItr != tangentEnd )
    {
    VectorType vectorToCenter = pointItr.Value() - this->m_SphereCenter;

    vectorToCenter.Normalize();

    const VectorType & tangent = tangentItr.Value();

    const double sinTheta = tangent.GetNorm();

    if( sinTheta > normEpsilon )
      {
      const double theta = vcl_asin( sinTheta );

      const VectorType axis = CrossProduct( vectorToCenter, tangent );

      versor.Set( axis, theta );

      dstPointItr.Value() = versor.Transform( pointItr.Value() );
      }
    else
      {
      dstPointItr.Value() = pointItr.Value();
      }

    ++dstPointItr;
    ++tangentItr;
    ++pointItr;
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::AssignResampledMovingValuesToOutputMesh()
{
  OutputMeshPointer outputMesh = this->GetOutput();

  OutputPointDataContainerPointer outputPointData = outputMesh->GetPointData();

  const PointIdentifier numberOfNodes = this->m_FixedMesh->GetNumberOfPoints();

  outputPointData->Reserve( numberOfNodes );

  OutputPointDataContainerIterator outputDataItr = outputPointData->Begin();

  ResampledMovingValuesContainerIterator resampledArrayItr = this->m_ResampledMovingValuesContainer->Begin();
  ResampledMovingValuesContainerIterator resampledArrayEnd = this->m_ResampledMovingValuesContainer->End();

  while( resampledArrayItr != resampledArrayEnd )
    {
    outputDataItr.Value() = resampledArrayItr.Value();

    ++outputDataItr;
    ++resampledArrayItr;
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
const
typename QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::FixedMeshType
* QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::
GetDeformedFixedMesh() const
  {
  if( this->GetNumberOfOutputs() < 2 )
    {
    return 0;
    }

  return dynamic_cast<const FixedMeshType *>(this->ProcessObject::GetOutput(1) );
  }

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComposeFixedMeshOutputDisplacedToMovingMesh()
{
  const FixedMeshType * in = this->m_FixedMesh.GetPointer();

  FixedMeshType * out =
    dynamic_cast<FixedMeshType *>(this->ProcessObject::GetOutput(1) );

  CopyMeshToMeshPoints( in, out );
  CopyMeshToMeshEdgeCells( in, out );
  CopyMeshToMeshCells( in, out );
  CopyMeshToMeshPointData( in, out );
  CopyMeshToMeshCellData( in, out );

  this->CopyDestinationPointsToDeformedFixedMesh();
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::CopyDestinationPointsToDeformedFixedMesh()
{
  DestinationPointConstIterator srcPointItr = this->m_DestinationPoints->Begin();
  DestinationPointConstIterator srcPointEnd = this->m_DestinationPoints->End();

  FixedMeshType * deformedFixedMesh =
    dynamic_cast<FixedMeshType *>(this->ProcessObject::GetOutput(1) );

  FixedPointsContainer * points = deformedFixedMesh->GetPoints();

  FixedPointsIterator dstPointItr = points->Begin();

  while( srcPointItr != srcPointEnd )
    {
    dstPointItr.Value().SetPoint( srcPointItr.Value() );
    ++dstPointItr;
    ++srcPointItr;
    }
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh,
                                               TOutputMesh>::ComposeDestinationPointsOutputPointSet()
{
  DestinationPointSetType * destinationPointSet =
    dynamic_cast<DestinationPointSetType *>(this->ProcessObject::GetOutput(2) );

  if( !destinationPointSet )
    {
    itkExceptionMacro("Problem found while composing the destination PointSet");
    }

  destinationPointSet->SetPoints( this->m_DestinationPoints );
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::PrintOutDeformationVectors(
  std::ostream & os )
{
  os << std::endl;
  os << "Deformation Vectors at every node " <<  std::endl;
  DestinationPointIterator dstPointItr = this->m_DestinationPoints->Begin();
  DestinationPointIterator dstPointEnd = this->m_DestinationPoints->End();

  const FixedPointsContainer * points = this->m_FixedMesh->GetPoints();
  FixedPointsConstIterator     srcPointItr = points->Begin();

  while( dstPointItr != dstPointEnd )
    {
    os <<  dstPointItr.Value() - srcPointItr.Value() << std::endl;

    ++dstPointItr;
    ++srcPointItr;
    }

  os << std::endl;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>
::ComputeLargestVelocityMagnitudeToShortestEdgeLengthRatio()
{
  double largestRatio = NumericTraits<double>::Zero;

  VelocityVectorConstIterator velocityItr = this->m_VelocityField->Begin();
  VelocityVectorConstIterator velocityEnd = this->m_VelocityField->End();

  ShortestLengthContainerIterator shortestEdgeItr = this->m_ShortestEdgeLengthPerPoint->Begin();

  while( velocityItr != velocityEnd )
    {
    const double velocityMagnitude = velocityItr.Value().GetNorm();

    const double ratio = velocityMagnitude / ( 2.0 * shortestEdgeItr.Value() );

    if( ratio > largestRatio )
      {
      largestRatio = ratio;
      }

    ++velocityItr;
    ++shortestEdgeItr;
    }

  this->m_LargestVelocityToEdgeLengthRatio = largestRatio;
}

template <class TFixedMesh, class TMovingMesh, class TOutputMesh>
void
QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TFixedMesh, TMovingMesh, TOutputMesh>::PrintSelf(std::ostream& os,
                                                                                                Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "User provided initial destination points : "
     << (this->m_UserProvidedInitialDestinationPoints ? " true " : " false " ) << std::endl;
  os << "Sphere center: " << this->m_SphereCenter << std::endl;
  os << "Sphere radius: " << this->m_SphereRadius << std::endl;
  os << "Epsilon : " << this->m_Epsilon << std::endl;
  os << "Sigma X : " << this->m_SigmaX << std::endl;
  os << "Lambda : " << this->m_Lambda << std::endl;
  os << "Maximum number of iterations : " << this->m_MaximumNumberOfIterations << std::endl;
  os << "Maximum number of smoothing iterations : " << this->m_MaximumNumberOfSmoothingIterations << std::endl;
  os << "Fixed nodes sigmas container pointer : " << this->m_FixedNodesSigmas.GetPointer() << std::endl;
  os << "Shortest edge length : " << this->m_ShortestEdgeLength << std::endl;
  os << "Scaling and squaring number of iterations : " << this->m_ScalingAndSquaringNumberOfIterations << std::endl;
}
}

#endif
