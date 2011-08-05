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

// ///////////// VTK Version Compatibility   //////////////////////////////
#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif
// ////////////////////////////////////////////////////////////////////////

namespace itk
{
/** \class DtiGuidedTrackingFilter
 */

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
GTRACT_COMMON_EXPORT class DtiGuidedTrackingFilter : public itk::DtiTrackingFilterBase<TTensorImageType,
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
  DtiGuidedTrackingFilter(const Self &); // purposely not implemented
  void operator=(const Self &);          // purposely not implemented

  bool GuideDirection(typename Self::ContinuousIndexType, GuideFiberType, const float, TVector &);

  GuideFiberType m_GuideFiber;
  double         m_CurvatureThreshold;
  double         m_GuidedCurvatureThreshold;
  double         m_MaximumGuideDistance;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDtiGuidedTrackingFilter.txx"
#endif

#endif
