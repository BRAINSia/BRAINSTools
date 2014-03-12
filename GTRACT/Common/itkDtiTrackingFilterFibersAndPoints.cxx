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
/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <iostream>
#include <fstream>

#include "itkDtiTrackingFilterBase.h"
#include "algo.h"
#include "GtractTypes.h"
#include "gtractCommonWin32.h"

#include <map>
#include <string>

void MinimumDistanceBetweenFiberGroups(VTKFiberListType fiberList1,
                                       VTKFiberListType fiberList2,
                                       TVector & dis)
{
  typedef itk::IdentityTransform<double> TransformType;
  TransformType::Pointer transform = TransformType::New();
  typedef itk::EuclideanDistancePointMetric<PointSetType, PointSetType> DistanceType;
  DistanceType::TransformParametersType parameter( transform->GetNumberOfParameters() );
  parameter.fill(0);

  VTKFiberListType::iterator fiberIt1, fiberIt2;
  int                        index = 0;
  for( fiberIt1 = fiberList1.begin(); fiberIt1 != fiberList1.end(); ++index, ++fiberIt1 )
    {
    vtkPolyData *      fiber1 = *fiberIt1;
    std::vector<float> avrDis;
    avrDis.clear();
    for( fiberIt2 = fiberList2.begin(); fiberIt2 != fiberList2.end(); ++fiberIt2 )
      {
      vtkPolyData *         fiber2 = *fiberIt2;
      PointSetType::Pointer pSet1 = PolyDataToPointSet(fiber1);
      PointSetType::Pointer pSet2 = PolyDataToPointSet(fiber2);
      DistanceType::Pointer distance = DistanceType::New();
      distance->SetMovingPointSet(pSet1);
      distance->SetFixedPointSet(pSet2);
      distance->SetTransform(transform);
      DistanceType::MeasureType value = distance->GetValue(parameter);
      int                       n = distance->GetNumberOfValues();
      float                     avr = 0;
      for( int i = 0; i < n; ++i )
        {
        avr += value[i];
        }
      avr /= n;
      avrDis.push_back(avr);
      }
    float min = 999;
    for( unsigned int i = 0; i < avrDis.size(); ++i )
      {
      if( min > avrDis[i] )
        {
        min = avrDis[i];
        }
      }
    dis[index] = min;
    }
}

PointSetType::Pointer PolyDataToPointSet(vtkPolyData *fiber)
{
  const int               npts = fiber->GetNumberOfPoints();
  PointSetType::Pointer   pSet = PointSetType::New();
  vtkFloatingPointType    p[3];
  PointSetType::PointType point;

  for( int i = 0; i < npts; ++i )
    {
    fiber->GetPoint(i, p);
    point[0] = p[0]; point[1] = p[1]; point[2] = p[2];
    pSet->SetPoint(i, point);
    }
  return pSet;
}
