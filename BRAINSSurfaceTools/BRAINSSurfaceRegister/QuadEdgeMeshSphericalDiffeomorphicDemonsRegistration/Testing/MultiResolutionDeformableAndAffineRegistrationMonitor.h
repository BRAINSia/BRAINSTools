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

#ifndef __MultiResolutionDeformableAndAffineRegistrationMonitor_h
#define __MultiResolutionDeformableAndAffineRegistrationMonitor_h

#include "DeformableAndAffineRegistrationMonitor.h"
#include <vector>

/** \class MultiResolutionDeformableAndAffineRegistrationMonitor
 *  This class provides a VTK visualization pipeline configured for monitoring
 *  the progress of a deformable registration process. It take multiple resolution
 *  surfaces and use them for the visualization of the proper resolution level.
 */
template <class TMultiResolutionDeformationFilter, class TPointSet>
class MultiResolutionDeformableAndAffineRegistrationMonitor :
  public DeformableAndAffineRegistrationMonitor<
    typename TMultiResolutionDeformationFilter::DeformationFilterType, TPointSet>
{
public:

  typedef MultiResolutionDeformableAndAffineRegistrationMonitor Self;
  typedef DeformableAndAffineRegistrationMonitor<
      typename TMultiResolutionDeformationFilter::DeformationFilterType, TPointSet> Superclass;

  MultiResolutionDeformableAndAffineRegistrationMonitor();
  virtual ~MultiResolutionDeformableAndAffineRegistrationMonitor();

  void SetNumberOfResolutionLevels( unsigned int number );

  void SetFixedSurface(unsigned int level, vtkPolyData * surface );

  void SetMovingSurface(unsigned int level, vtkPolyData * surface );

  void SetResolutionLevel(unsigned int level);

private:
  typedef std::vector<vtkPolyData *> PolyDataArray;

  PolyDataArray m_FixedSurfaces;
  PolyDataArray m_MovingSurfaces;
};

#include "MultiResolutionDeformableAndAffineRegistrationMonitor.txx"

#endif
