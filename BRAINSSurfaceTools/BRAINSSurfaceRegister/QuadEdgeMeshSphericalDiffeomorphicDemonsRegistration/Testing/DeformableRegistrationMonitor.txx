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

#include "DeformableRegistrationMonitor.h"

#include "itkMacro.h"
#include "vtkPoints.h"

/** Constructor */
template <class TDeformationFilter, class TPointSet>
DeformableRegistrationMonitor<TDeformationFilter, TPointSet>
::DeformableRegistrationMonitor()
{
}

/** Destructor */
template <class TDeformationFilter, class TPointSet>
DeformableRegistrationMonitor<TDeformationFilter, TPointSet>
::~DeformableRegistrationMonitor()
{
}

/** Set the objects to be observed: optimizer and transform */
template <class TDeformationFilter, class TPointSet>
void
DeformableRegistrationMonitor<TDeformationFilter, TPointSet>
::ObserveData( const PointSetType * destinationPointSet )
{
  this->ObservedPointSet = destinationPointSet;
}

/** Update the Visualization */
template <class TDeformationFilter, class TPointSet>
void
DeformableRegistrationMonitor<TDeformationFilter, TPointSet>
::UpdateDataBeforeRendering()
{
  typedef typename PointSetType::PointsContainer  PointsContainer;
  typedef typename PointsContainer::ConstIterator PointsIterator;
  typedef typename PointSetType::PointType        PointType;

  const PointsContainer * displacedPoints = this->ObservedPointSet->GetPoints();

  if( !displacedPoints )
    {
    std::cerr << "displacedPoints = NULL" << std::endl;
    return;
    }

  vtkPoints * vtkpoints = this->GetFixedSurfacePoints();

  if( !vtkpoints )
    {
    std::cerr << "vtkpoints = NULL" << std::endl;
    return;
    }

  unsigned int numberOfPoints =
    static_cast<unsigned int>( vtkpoints->GetNumberOfPoints() );

  if( numberOfPoints != displacedPoints->Size() )
    {
    std::cerr << "vtkPolyData numberOfPoints  = " << numberOfPoints << std::endl;
    std::cerr << "number of deformable points = " << displacedPoints->Size() << std::endl;
    itkGenericExceptionMacro("Fixed Mesh and Deformed Points have different size");
    }

  PointsIterator pointItr = displacedPoints->Begin();
  PointsIterator pointEnd = displacedPoints->End();

  vtkFloatingPointType displacedPoint[3];

  int pointId = 0;

  while( pointItr != pointEnd )
    {
    const PointType & observedDisplacedPoint = pointItr.Value();

    displacedPoint[0] = observedDisplacedPoint[0];
    displacedPoint[1] = observedDisplacedPoint[1];
    displacedPoint[2] = observedDisplacedPoint[2];

    vtkpoints->SetPoint( pointId, displacedPoint );

    ++pointItr;
    ++pointId;
    }

  this->MarkFixedSurfaceAsModified();
}

/** Print out iteration information */
template <class TDeformationFilter, class TPointSet>
void
DeformableRegistrationMonitor<TDeformationFilter, TPointSet>
::PrintOutUpdateMessage()
{
  itk::Object * observedObject = this->GetObservedObject();

  TDeformationFilter * deformationFilter =
    dynamic_cast<TDeformationFilter *>( observedObject );

  if( deformationFilter != NULL )
    {
    std::cout << " Metric = " << deformationFilter->GetMetricValue() << std::endl;
    }
}
