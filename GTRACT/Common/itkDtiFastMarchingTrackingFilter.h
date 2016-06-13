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

#ifndef __itkDtiFastMarchingTrackingFilter_h
#define __itkDtiFastMarchingTrackingFilter_h

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <queue>

#include <itkImage.h>
#include <metaCommand.h>
#include <itkObject.h>
#include <itkImageToImageFilter.h>
#include <itkIOCommon.h>
#include <itkIndex.h>
#include <vnl/vnl_math.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkSingleValuedCostFunction.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkPointSet.h>
#include <vtkPoints.h>

#include "gtractCommonWin32.h"
#include "GtractTypes.h"
#include "itkFastMarchingCostFunction.h"
#include "itkDtiTrackingFilterBase.h"

namespace itk
{
/** \class DtiFastMarchingTrackingFilter
 */

template <class TTensorImageType, class TAnisotropyImageType, class TCostImageType, class TMaskImageType>
class GTRACT_COMMON_EXPORT DtiFastMarchingTrackingFilter : public itk::DtiTrackingFilterBase<TTensorImageType,
                                                                                             TAnisotropyImageType,
                                                                                             TMaskImageType>
{
public:
  /** Standard class typdedefs. */
  typedef DtiFastMarchingTrackingFilter                                                      Self;
  typedef itk::DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType> Superclass;
  typedef SmartPointer<Self>                                                                 Pointer;
  typedef SmartPointer<const Self>                                                           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DtiFastMarchingTrackingFilter, itk::DtiTrackingFilterBase);

  /** Typedef support of input Cost Image Type */
  typedef TCostImageType                        CostImageType;
  typedef typename CostImageType::Pointer       CostImagePointer;
  typedef typename CostImageType::ConstPointer  CostImageConstPointer;
  typedef typename CostImageType::RegionType    CostImageRegionType;
  typedef typename CostImageType::SizeType      CostImageSizeType;
  typedef typename CostImageType::SpacingType   CostImageSpacingType;
  typedef typename CostImageType::PointType     CostImagePointType;
  typedef typename CostImageType::PixelType     CostImagePixelType;
  typedef typename CostImageType::DirectionType CostImageDirectionType;

  typedef typename itk::LinearInterpolateImageFunction<CostImageType, double> CostIPType;
  typedef typename CostIPType::ContinuousIndexType                            ContinuousIndexType;

  // Setup the CostFunction

  typedef itk::FastMarchingCostFunction    CostFunctionType;
  typedef CostFunctionType::Pointer        CostFunctionPointer;
  typedef CostFunctionType::ParametersType ParametersType;

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef OptimizerType::ScalesType                ScalesType;
  typedef OptimizerType::DerivativeType            DerivativeType;

  typedef typename std::list<ContinuousIndexType>
    StartPointsListType;

  itkSetObjectMacro( CostFN, CostFunctionType );
  itkGetConstObjectMacro( CostFN, CostFunctionType );
  itkSetObjectMacro(CostImage,  CostImageType);
  itkGetConstObjectMacro(CostImage,  CostImageType);

  itkSetMacro(MaxStepSize, double);               // for Gradient Descent with
                                                  // default set to 1.0
  itkSetMacro(MinStepSize, double);               // for Gradient Descent with
                                                  // default set to 0.01
  itkSetMacro(NumberOfIterations, unsigned long)  // for Gradient Descent with
  // default set to 150;

  /* For Cost Funtion neighborhood iterator; Size 1 default with smaller size
    for partial voxel */
  itkSetMacro(CostFunctionStepSize, float);

  void Update(); // InitializeSeeds() and starts data generation

protected:
  DtiFastMarchingTrackingFilter();
  ~DtiFastMarchingTrackingFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DtiFastMarchingTrackingFilter);

  void PrintSelf( std::ostream & os, Indent indent ) const ITK_OVERRIDE;

  void InitializeSeeds();

  void GradientDescent( ContinuousIndexType & index);

  // Input and Output Image
  CostImagePointer    m_CostImage;
  CostFunctionPointer m_CostFN;
  StartPointsListType m_StartPoints;

  typename CostIPType::Pointer m_CostIP;
  OptimizerType::Pointer m_GradientOP;

  double        m_MaxStepSize;
  double        m_MinStepSize;
  float         m_CostFunctionStepSize;
  unsigned long m_NumberOfIterations;
}; // end class
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDtiFastMarchingTrackingFilter.hxx"
#endif

#endif
