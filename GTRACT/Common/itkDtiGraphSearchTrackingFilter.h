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

#ifndef __itkDtiGraphSearchTrackingFilter_h
#define __itkDtiGraphSearchTrackingFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkBlobSpatialObject.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "itkDtiTrackingFilterBase.h"
#include "algo.h"
#include "GtractTypes.h"
#include "gtractCommonWin32.h"

#include <map>
#include <string>

namespace itk
{
/** \class DtiGraphSearchTrackingFilter
 */

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
GTRACT_COMMON_EXPORT class DtiGraphSearchTrackingFilter : public itk::DtiTrackingFilterBase<TTensorImageType,
                                                                                            TAnisotropyImageType,
                                                                                            TMaskImageType>
{
public:
  /** Standard class typedefs. */
  typedef DtiGraphSearchTrackingFilter                                                       Self;
  typedef itk::DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType> Superclass;
  typedef SmartPointer<Self>                                                                 Pointer;
  typedef SmartPointer<const Self>                                                           ConstPointer;

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomGeneratorType;
  typedef RandomGeneratorType::Pointer                           RandomGeneratorPointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DtiGraphSearchTrackingFilter, itk::DtiTrackingFilterBase);

  itkSetMacro(AnisotropyBranchingValue, float);
  itkSetMacro(RandomSeed, int);
  itkSetMacro(MaximumBranches, unsigned int);
  itkSetMacro(UseRandomWalk, bool);

  itkSetMacro(RandomWalkAngle, double);
  itkGetMacro(RandomWalkAngle, double);
  itkSetMacro(CurvatureBranchAngle, double);
  itkGetMacro(CurvatureBranchAngle, double);

  typename itk::Point<double, 3> InitializeCenterOfMask();

  void Update();

protected:
  DtiGraphSearchTrackingFilter();
  ~DtiGraphSearchTrackingFilter()
  {
  }

private:
  DtiGraphSearchTrackingFilter(const Self &); // purposely not implemented
  void operator=(const Self &);               // purposely not implemented

  RandomGeneratorPointer m_RandomGenerator;

  float        m_AnisotropyBranchingValue;
  double       m_CurvatureBranchAngle;
  unsigned int m_MaximumBranches;
  bool         m_UseRandomWalk;
  double       m_RandomWalkAngle;
  int          m_RandomSeed;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDtiGraphSearchTrackingFilter.txx"
#endif

#endif
