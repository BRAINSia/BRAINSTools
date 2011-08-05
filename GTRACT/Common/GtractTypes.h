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

#ifndef _GtractTypes_h
#define _GtractTypes_h

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <itkContinuousIndex.h>

typedef vnl_matrix<float> TMatrix;
typedef vnl_vector<float> TVector;

class BranchPointType
{
public:
  TVector m_Direction;
  int     m_DivergePoint;
};

class FiberPointType
{
public:
  itk::ContinuousIndex<float, 3> m_Point;
  float                          m_AI;
  float                          m_AISum;
};

#endif
