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

#ifndef __itkAnatomicalVersorRigidFilter_h
#define __itkAnatomicalVersorRigidFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include <itkExtractImageFilter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkVersorRigid3DTransformOptimizer.h>
#include <itkImageRegistrationMethod.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkVersorRigid3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFactory.h>
#include "gtractCommonWin32.h"

#include <map>
#include <string>

namespace itk
{
/** \class DtiToAnatomicalRigidRegistrationFilter
 * \brief Rigid registration convience class between
 * a 4D image time series and a 3D anatomical image.
 *
 * The Versor Rigid registration is used along with the
 * Mattes mutual information registration. The output
 * of the filter is the resulting registration.
 *
 */

class GTRACT_COMMON_EXPORT AnatomicalVersorRigidFilter : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( AnatomicalVersorRigidFilter );

  /** Standard class type alias. */
  using Self = AnatomicalVersorRigidFilter;
  using Superclass = itk::Object;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Fixed Image type alias. */
  using RegisterImageType = itk::Image< signed short, 3 >;
  using RegisterImagePointer = RegisterImageType::Pointer;
  using RegisterImageConstPointer = RegisterImageType::ConstPointer;
  using RegisterImageRegionType = RegisterImageType::RegionType;
  using RegisterImageSizeType = RegisterImageType::SizeType;
  using RegisterImageSpacingType = RegisterImageType::SpacingType;
  using RegisterImagePointType = RegisterImageType::PointType;
  using RegisterImagePixelType = RegisterImageType::PixelType;
  using RegisterImageDirectionType = RegisterImageType::DirectionType;
  using RegisterImageIndexType = RegisterImageType::IndexType;

  /** Output Transform type alias. */
  using TransformType = itk::VersorRigid3DTransform< double >;
  using OptimizerType = itk::VersorRigid3DTransformOptimizer;
  using MetricType = itk::MattesMutualInformationImageToImageMetric< RegisterImageType, RegisterImageType >;

  using InterpolatorType = itk::LinearInterpolateImageFunction< RegisterImageType, double >;

  using RegistrationType = itk::ImageRegistrationMethod< RegisterImageType, RegisterImageType >;

  typedef itk::CenteredTransformInitializer< TransformType, RegisterImageType, RegisterImageType >
    TransformInitializerType;
  using TransformTypePointer = TransformType::Pointer;
  using VersorType = TransformType::VersorType;
  using VectorType = VersorType::VectorType;
  using MetricTypePointer = MetricType::Pointer;
  using OptimizerTypePointer = OptimizerType::Pointer;
  using OptimizerParameterType = OptimizerType::ParametersType;
  using OptimizerScalesType = OptimizerType::ScalesType;
  using InterpolatorTypePointer = InterpolatorType::Pointer;
  using RegistrationTypePointer = RegistrationType::Pointer;
  using TransformInitializerTypePointer = TransformInitializerType::Pointer;

  /** Standard New method. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( AnatomicalVersorRigidFilter, itk::Object );

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro( FixedImage, RegisterImageType );
  itkSetObjectMacro( MovingImage, RegisterImageType );
  itkGetConstObjectMacro( Output, TransformType );

  itkSetMacro( NumberOfSpatialSamples, int );
  itkSetMacro( NumberOfIterations, int );
  itkSetMacro( TranslationScale, float );
  itkSetMacro( MaximumStepLength, float );
  itkSetMacro( MinimumStepLength, float );
  itkSetMacro( RelaxationFactor, float );
  itkSetMacro( InitialRotationAngle, double );
  itkSetMacro( InitialRotationAxis, int );

  void
  Update();

protected:
  AnatomicalVersorRigidFilter();
  ~AnatomicalVersorRigidFilter() override {}

private:
  // Input and Output Image
  RegisterImagePointer m_FixedImage;
  RegisterImagePointer m_MovingImage;
  TransformTypePointer m_Output;

  // Registration Parameters
  float  m_TranslationScale;
  float  m_MaximumStepLength;
  float  m_MinimumStepLength;
  float  m_RelaxationFactor;
  int    m_NumberOfSpatialSamples;
  int    m_NumberOfIterations;
  int    m_InitialRotationAxis;
  double m_InitialRotationAngle;
}; // end of class
} // end namespace itk

#endif
