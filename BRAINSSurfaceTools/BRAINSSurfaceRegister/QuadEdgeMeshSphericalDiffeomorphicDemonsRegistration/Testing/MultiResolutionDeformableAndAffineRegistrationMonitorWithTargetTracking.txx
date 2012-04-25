/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTestingMacros.h,v $
  Language:  C++
  Date:      $Date: 2009-05-09 17:40:20 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking.h"

#include "itkMacro.h"
#include "itkDeformQuadEdgeMeshFilter.h"

#include "vtkPoints.h"

/** Constructor */
template <class TMultiResolutionDeformationFilter, class TPointSet>
MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<TMultiResolutionDeformationFilter, TPointSet>
::MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking()
{
  this->m_EvaluateDistanceToTarget = false;

  this->m_MultiResolutionDemonsRegistrationFilter = NULL;
}

/** Destructor */
template <class TMultiResolutionDeformationFilter, class TPointSet>
MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<TMultiResolutionDeformationFilter, TPointSet>
::~MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking()
{
}

template <class TMultiResolutionDeformationFilter, class TPointSet>
void
MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<TMultiResolutionDeformationFilter, TPointSet>
::SetMultiResolutionDemonsFilter( TMultiResolutionDeformationFilter * filter )
{
  this->m_MultiResolutionDemonsRegistrationFilter = filter;
  this->m_MultiResolutionDemonsRegistrationFilter->SetRegistrationMonitor( this );
}

template <class TMultiResolutionDeformationFilter, class TPointSet>
void
MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<TMultiResolutionDeformationFilter, TPointSet>
::SetEvaluateDistanceToTarget( bool value )
{
  this->m_EvaluateDistanceToTarget = value;
}

template <class TMultiResolutionDeformationFilter, class TPointSet>
bool
MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<TMultiResolutionDeformationFilter, TPointSet>
::GetEvaluateDistanceToTarget() const
{
  return this->m_EvaluateDistanceToTarget;
}

template <class TMultiResolutionDeformationFilter, class TPointSet>
void
MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<TMultiResolutionDeformationFilter, TPointSet>
::PrintOutUpdateMessage()
{
  this->Superclass::PrintOutUpdateMessage();

  if( this->m_MultiResolutionDemonsRegistrationFilter->GetRegistrationMode() ==
      TMultiResolutionDeformationFilter::DEFORMABLE )
    {
    this->EvaluateDistanceToTarget();
    }
}

template <class TMultiResolutionDeformationFilter, class TPointSet>
void
MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<TMultiResolutionDeformationFilter, TPointSet>
::EvaluateDistanceToTarget() const
{
  if( !this->m_EvaluateDistanceToTarget )
    {
    return;
    }

  typedef typename DeformationFilterType::FixedMeshType           FixedMeshType;
  typedef typename DeformationFilterType::MovingMeshType          MovingMeshType;
  typedef typename DeformationFilterType::DestinationPointSetType DestinationPointSetType;

  //
  //   Deform source fixed new mesh using current fixed mesh and destination points
  //
  typedef itk::DeformQuadEdgeMeshFilter<FixedMeshType, FixedMeshType, DestinationPointSetType> DeformFilterType;

  typename DeformFilterType::Pointer deformFilter = DeformFilterType::New();

  deformFilter->SetInput( this->m_FixedMeshSource );

  const FixedMeshType * referenceMesh =
    this->m_MultiResolutionDemonsRegistrationFilter->GetCurrentLevelInitialFixedMesh();

  const DestinationPointSetType * destinationPointSet =
    this->m_MultiResolutionDemonsRegistrationFilter->GetCurrentDestinationPoints();

  if( referenceMesh->GetNumberOfPoints() != destinationPointSet->GetNumberOfPoints() )
    {
    return;
    }

  deformFilter->SetReferenceMesh( referenceMesh );
  deformFilter->SetDestinationPoints( destinationPointSet );

  deformFilter->SetSphereRadius( this->m_MultiResolutionDemonsRegistrationFilter->GetSphereRadius() );
  deformFilter->SetSphereCenter( this->m_MultiResolutionDemonsRegistrationFilter->GetSphereCenter() );

  try
    {
    deformFilter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    throw excp;
    }

  const FixedMeshType * mappedFixedSource = deformFilter->GetOutput();

  typedef typename DeformationFilterType::FixedPointsConstIterator FixedPointsConstIterator;

  FixedPointsConstIterator sourcePointItr   = this->m_FixedMeshSource->GetPoints()->Begin();
  FixedPointsConstIterator sourcePointEnd   = this->m_FixedMeshSource->GetPoints()->End();

  FixedPointsConstIterator mappedPointItr   = mappedFixedSource->GetPoints()->Begin();

  FixedPointsConstIterator targetPointItr   = this->m_FixedMeshTarget->GetPoints()->Begin();

  double sumOfDistances1 = 0.0;
  double sumOfDistances2 = 0.0;

  const unsigned long numberOfPoints = mappedFixedSource->GetNumberOfPoints();

  while( sourcePointItr != sourcePointEnd )
    {
    const double distanceSquared1 =
      sourcePointItr.Value().SquaredEuclideanDistanceTo( mappedPointItr.Value() );

    sumOfDistances1 += distanceSquared1;

    const double distanceSquared2 =
      mappedPointItr.Value().SquaredEuclideanDistanceTo( targetPointItr.Value() );

    sumOfDistances2 += distanceSquared2;

    ++sourcePointItr;
    ++mappedPointItr;
    ++targetPointItr;
    }

  const double distancesRMS1 = vcl_sqrt( sumOfDistances1 / numberOfPoints );
  const double distancesRMS2 = vcl_sqrt( sumOfDistances2 / numberOfPoints );

  std::cout << " RMS to source " << distancesRMS1 << "  RMS to target " << distancesRMS2 << std::endl;
}
