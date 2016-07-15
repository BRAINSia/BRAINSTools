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

#ifndef __itkAnatomicalBSplineFilter_h
#define __itkAnatomicalBSplineFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include <itkExtractImageFilter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkImageRegistrationMethod.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineDeformableTransform.h>
#include <itkLBFGSBOptimizer.h>
#include <itkCenteredTransformInitializer.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFactory.h>
#include "gtractCommonWin32.h"

#include <map>
#include <string>

namespace itk
{
/** \class DtiToAnatomicalBSplineRegistrationFilter
 * \brief This is a convience class to perform non linear image registration
 *  between a 3D image set and a 4D image set. The input
 *  images must be signed short. In future revisions, the
 *  class will be modified to support vector images instead
 *  of 4D images.
 *
 * For image registration the B-Spline Registration is used with the
 * MattesMutualInformationImageToImageMetric. A bulk transform can be
 * specified proving an initial starting point for the deformable
 * registration.
 */

class GTRACT_COMMON_EXPORT AnatomicalBSplineFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef AnatomicalBSplineFilter  Self;
  typedef itk::Object              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Fixed Image typedefs. */
  typedef itk::Image<signed short, 3>      RegisterImageType;
  typedef RegisterImageType::Pointer       RegisterImagePointer;
  typedef RegisterImageType::ConstPointer  RegisterImageConstPointer;
  typedef RegisterImageType::RegionType    RegisterImageRegionType;
  typedef RegisterImageType::SizeType      RegisterImageSizeType;
  typedef RegisterImageType::SpacingType   RegisterImageSpacingType;
  typedef RegisterImageType::PointType     RegisterImagePointType;
  typedef RegisterImageType::PixelType     RegisterImagePixelType;
  typedef RegisterImageType::DirectionType RegisterImageDirectionType;
  typedef RegisterImageType::IndexType     RegisterImageIndexType;

  static const unsigned int SpaceDimension = 3;
  static const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;
  typedef itk::BSplineDeformableTransform<
      CoordinateRepType,
      SpaceDimension,
      SplineOrder>               TransformType;
  typedef TransformType::RegionType     TransformRegionType;
  typedef TransformRegionType::SizeType TransformSizeType;
  typedef TransformType::SpacingType    TransformSpacingType;
  typedef TransformType::OriginType     TransformOriginType;
  typedef TransformType::ParametersType TransformParametersType;

  typedef itk::LBFGSBOptimizer OptimizerType;

  typedef itk::MattesMutualInformationImageToImageMetric<
      RegisterImageType,
      RegisterImageType>          MetricType;

  typedef itk::LinearInterpolateImageFunction<
      RegisterImageType,
      double>        InterpolatorType;

  typedef itk::ImageRegistrationMethod<
      RegisterImageType,
      RegisterImageType>          RegistrationType;

  typedef TransformType::Pointer            TransformTypePointer;
  typedef MetricType::Pointer               MetricTypePointer;
  typedef OptimizerType::Pointer            OptimizerTypePointer;
  typedef OptimizerType::ParametersType     OptimizerParameterType;
  typedef OptimizerType::ScalesType         OptimizerScalesType;
  typedef OptimizerType::BoundSelectionType OptimizerBoundSelectionType;
  typedef OptimizerType::BoundValueType     OptimizerBoundValueType;

  typedef InterpolatorType::Pointer InterpolatorTypePointer;
  typedef RegistrationType::Pointer RegistrationTypePointer;

  /** Typedef of the bulk transform. */
  /*
  typedef itk::VersorRigid3DTransform< double >     BulkTransformType;
  typedef BulkTransformType::Pointer           BulkTransformPointer;
   */
  typedef Transform<CoordinateRepType, itkGetStaticConstMacro(SpaceDimension),
                    itkGetStaticConstMacro(SpaceDimension)> BulkTransformType;
  typedef BulkTransformType::ConstPointer BulkTransformPointer;
  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(AnatomicalBSplineFilter, itk::Object);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(FixedImage, RegisterImageType);
  itkSetObjectMacro(MovingImage, RegisterImageType);
  itkSetConstObjectMacro(BulkTransform, BulkTransformType);
  itkGetConstObjectMacro(Output, TransformType);

  itkSetMacro(SpatialSampleScale, int);
  itkSetMacro(MaximumNumberOfIterations, int);
  itkSetMacro(MaximumNumberOfEvaluations, int);
  itkSetMacro(MaximumNumberOfCorrections, int);
  itkSetMacro(BSplineHistogramBins, int);
  itkSetMacro(GridSize, TransformSizeType);
  itkSetMacro(GridBorderSize, int);
  itkSetMacro(CostFunctionConvergenceFactor, float);
  itkSetMacro(ProjectedGradientTolerance, float);

  itkSetMacro(BoundTypeX, int);
  itkGetMacro(BoundTypeX, int);
  itkSetMacro(BoundTypeY, int);
  itkGetMacro(BoundTypeY, int);
  itkSetMacro(BoundTypeZ, int);
  itkGetMacro(BoundTypeZ, int);
  itkSetMacro(LowerBoundX, float);
  itkGetMacro(LowerBoundX, float);
  itkSetMacro(LowerBoundY, float);
  itkGetMacro(LowerBoundY, float);
  itkSetMacro(LowerBoundZ, float);
  itkGetMacro(LowerBoundZ, float);
  itkSetMacro(UpperBoundX, float);
  itkGetMacro(UpperBoundX, float);
  itkSetMacro(UpperBoundY, float);
  itkGetMacro(UpperBoundY, float);
  itkSetMacro(UpperBoundZ, float);
  itkGetMacro(UpperBoundZ, float);

  void Update();

protected:
  AnatomicalBSplineFilter();
  ~AnatomicalBSplineFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(AnatomicalBSplineFilter);

  /*** Input and Output Objects ***/
  RegisterImagePointer m_FixedImage;
  RegisterImagePointer m_MovingImage;
  BulkTransformPointer m_BulkTransform;
  TransformTypePointer m_Output;

  // Parameters for the image registration
  int               m_SpatialSampleScale;
  int               m_MaximumNumberOfIterations;
  int               m_MaximumNumberOfEvaluations;
  int               m_MaximumNumberOfCorrections;
  int               m_BSplineHistogramBins;
  TransformSizeType m_GridSize;
  int               m_GridBorderSize;
  float             m_CostFunctionConvergenceFactor;
  float             m_ProjectedGradientTolerance;

  int   m_BoundTypeX;
  int   m_BoundTypeY;
  int   m_BoundTypeZ;
  float m_LowerBoundX;
  float m_LowerBoundY;
  float m_LowerBoundZ;
  float m_UpperBoundX;
  float m_UpperBoundY;
  float m_UpperBoundZ;
};  // end of class
} // end namespace itk

#endif
