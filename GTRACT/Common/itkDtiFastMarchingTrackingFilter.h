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
#include <itkMath.h>
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

template <typename TTensorImageType, typename TAnisotropyImageType, typename TCostImageType, typename TMaskImageType>
class GTRACT_COMMON_EXPORT DtiFastMarchingTrackingFilter
  : public itk::DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(DtiFastMarchingTrackingFilter);

  /** Standard class typdedefs. */
  using Self = DtiFastMarchingTrackingFilter;
  using Superclass = itk::DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DtiFastMarchingTrackingFilter, itk::DtiTrackingFilterBase);

  /** Typedef support of input Cost Image Type */
  using CostImageType = TCostImageType;
  using CostImagePointer = typename CostImageType::Pointer;
  using CostImageConstPointer = typename CostImageType::ConstPointer;
  using CostImageRegionType = typename CostImageType::RegionType;
  using CostImageSizeType = typename CostImageType::SizeType;
  using CostImageSpacingType = typename CostImageType::SpacingType;
  using CostImagePointType = typename CostImageType::PointType;
  using CostImagePixelType = typename CostImageType::PixelType;
  using CostImageDirectionType = typename CostImageType::DirectionType;

  using CostIPType = typename itk::LinearInterpolateImageFunction<CostImageType, double>;
  using ContinuousIndexType = typename CostIPType::ContinuousIndexType;

  // Setup the CostFunction

  using CostFunctionType = itk::FastMarchingCostFunction;
  using CostFunctionPointer = CostFunctionType::Pointer;
  using ParametersType = CostFunctionType::ParametersType;

  using OptimizerType = itk::RegularStepGradientDescentOptimizer;
  using ScalesType = OptimizerType::ScalesType;
  using DerivativeType = OptimizerType::DerivativeType;

  using StartPointsListType = typename std::list<ContinuousIndexType>;

  itkSetObjectMacro(CostFN, CostFunctionType);
  itkGetConstObjectMacro(CostFN, CostFunctionType);
  itkSetObjectMacro(CostImage, CostImageType);
  itkGetConstObjectMacro(CostImage, CostImageType);

  itkSetMacro(MaxStepSize, double);              // for Gradient Descent with
                                                 // default set to 1.0
  itkSetMacro(MinStepSize, double);              // for Gradient Descent with
                                                 // default set to 0.01
  itkSetMacro(NumberOfIterations, unsigned long) // for Gradient Descent with
                                                 // default set to 150;

    /* For Cost Funtion neighborhood iterator; Size 1 default with smaller size
      for partial voxel */
    itkSetMacro(CostFunctionStepSize, float);

  void
  Update(); // InitializeSeeds() and starts data generation

protected:
  DtiFastMarchingTrackingFilter();
  ~DtiFastMarchingTrackingFilter() override = default;

private:
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  InitializeSeeds();

  void
  GradientDescent(ContinuousIndexType & index);

  // Input and Output Image
  CostImagePointer    m_CostImage;
  CostFunctionPointer m_CostFN;
  StartPointsListType m_StartPoints;

  typename CostIPType::Pointer m_CostIP;
  OptimizerType::Pointer       m_GradientOP;

  double        m_MaxStepSize;
  double        m_MinStepSize;
  float         m_CostFunctionStepSize;
  unsigned long m_NumberOfIterations;
}; // end class
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDtiFastMarchingTrackingFilter.hxx"
#endif

#endif
