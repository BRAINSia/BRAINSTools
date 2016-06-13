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

#ifndef __itkDtiGuidedTrackingFilter_h
#define __itkDtiGuidedTrackingFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkBlobSpatialObject.h"

#include "itkDtiTrackingFilterBase.h"
#include "algo.h"
#include "GtractTypes.h"
#include "gtractCommonWin32.h"

#include <map>
#include <string>

// #include "vtkPoints.h"

// ////////////////////////////////////////////////////////////////////////

namespace itk
{
/** \class DtiGuidedTrackingFilter
 */

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
class DtiGuidedTrackingFilter : public itk::DtiTrackingFilterBase<TTensorImageType,
                                                                  TAnisotropyImageType,
                                                                  TMaskImageType>
{
public:
  /** Standard class typedefs. */
  typedef DtiGuidedTrackingFilter                                                            Self;
  typedef itk::DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType> Superclass;
  typedef SmartPointer<Self>                                                                 Pointer;
  typedef SmartPointer<const Self>                                                           ConstPointer;

  typedef vtkPolyData *GuideFiberType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DtiGuidedTrackingFilter, itk::DtiTrackingFilterBase);

  itkSetMacro(CurvatureThreshold, double);
  itkGetMacro(CurvatureThreshold, double);
  itkSetMacro(GuidedCurvatureThreshold, double);
  itkGetMacro(GuidedCurvatureThreshold, double);
  itkSetMacro(MaximumGuideDistance, double);
  itkGetMacro(MaximumGuideDistance, double);

  // void InitializeSeeds();
  void SetGuideFiber( GuideFiberType );

  void Update();

protected:
  DtiGuidedTrackingFilter();
  ~DtiGuidedTrackingFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DtiGuidedTrackingFilter);

  bool GuideDirection(typename Self::ContinuousIndexType, GuideFiberType, const float, TVector &);

  GuideFiberType m_GuideFiber;
  double         m_CurvatureThreshold;
  double         m_GuidedCurvatureThreshold;
  double         m_MaximumGuideDistance;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDtiGuidedTrackingFilter.hxx"
#endif

#endif
