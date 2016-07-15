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

#ifndef __itkDtiStreamlineTrackingFilter_h
#define __itkDtiStreamlineTrackingFilter_h

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

namespace itk
{
/** \class DtiStreamlineTrackingFilter
 */

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
class DtiStreamlineTrackingFilter : public itk::DtiTrackingFilterBase<TTensorImageType,
                                                                      TAnisotropyImageType,
                                                                      TMaskImageType>
{
public:
  /** Standard class typedefs. */
  typedef DtiStreamlineTrackingFilter                                                        Self;
  typedef itk::DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType> Superclass;
  typedef SmartPointer<Self>                                                                 Pointer;
  typedef SmartPointer<const Self>                                                           ConstPointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DtiStreamlineTrackingFilter, itk::DtiTrackingFilterBase);

  itkSetMacro(CurvatureThreshold, double);
  itkGetMacro(CurvatureThreshold, double);

  // void SetSeeds(SeedListType);
  // void SetTrackingDirections(DirectionListType);

  void Update();

protected:
  DtiStreamlineTrackingFilter();
  ~DtiStreamlineTrackingFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DtiStreamlineTrackingFilter);

  double m_CurvatureThreshold;
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDtiStreamlineTrackingFilter.hxx"
#endif

#endif
