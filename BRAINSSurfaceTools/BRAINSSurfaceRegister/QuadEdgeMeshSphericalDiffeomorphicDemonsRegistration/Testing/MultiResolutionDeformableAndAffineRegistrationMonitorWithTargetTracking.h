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

#ifndef __MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking_h
#define __MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking_h

#include "MultiResolutionDeformableAndAffineRegistrationMonitor.h"

/** \class MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking
 *  This class provides a VTK visualization pipeline configured for monitoring
 *  the progress of a deformable registration process. It take multiple resolution
 *  surfaces and use them for the visualization of the proper resolution level.
 *  It also takes two surfaces as source and target and track the progress of the
 *  destination points from the source to the target with an RMS measure of distances.
 */
template <class TMultiResolutionDeformationFilter, class TPointSet>
class MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking :
  public MultiResolutionDeformableAndAffineRegistrationMonitor<
    TMultiResolutionDeformationFilter, TPointSet>
{
public:

  typedef MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking Self;
  typedef MultiResolutionDeformableAndAffineRegistrationMonitor<
      TMultiResolutionDeformationFilter, TPointSet> Superclass;

  MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking();
  virtual ~MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking();

  typedef typename TMultiResolutionDeformationFilter::MeshType MeshType;

  void SetFixedMeshSource( const MeshType * mesh )
  {
    this->m_FixedMeshSource = mesh;
  }

  void SetFixedMeshTarget( const MeshType * mesh )
  {
    this->m_FixedMeshTarget = mesh;
  }

  void SetEvaluateDistanceToTarget( bool );

  bool GetEvaluateDistanceToTarget() const;

  void EvaluateDistanceToTarget() const;

  void PrintOutUpdateMessage();

  void SetMultiResolutionDemonsFilter( TMultiResolutionDeformationFilter * filter );

private:

  typename MeshType::ConstPointer     m_FixedMeshSource;
  typename MeshType::ConstPointer     m_FixedMeshTarget;
  bool m_EvaluateDistanceToTarget;

  typedef typename TMultiResolutionDeformationFilter::DeformationFilterType DeformationFilterType;

  typename TMultiResolutionDeformationFilter::Pointer   m_MultiResolutionDemonsRegistrationFilter;
};

#include "MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking.txx"

#endif
