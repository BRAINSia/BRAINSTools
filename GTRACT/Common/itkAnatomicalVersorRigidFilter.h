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
  /** Standard class typedefs. */
  typedef AnatomicalVersorRigidFilter Self;
  typedef itk::Object                 Superclass;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer<const Self>    ConstPointer;

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

  /** Output Transform typedefs. */
  typedef itk::VersorRigid3DTransform<double>  TransformType;
  typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetric<
      RegisterImageType,
      RegisterImageType>        MetricType;

  typedef itk::LinearInterpolateImageFunction<
      RegisterImageType,
      double>         InterpolatorType;

  typedef itk::ImageRegistrationMethod<
      RegisterImageType,
      RegisterImageType>        RegistrationType;

  typedef itk::CenteredTransformInitializer<TransformType,
                                            RegisterImageType,
                                            RegisterImageType
                                            >  TransformInitializerType;
  typedef TransformType::Pointer            TransformTypePointer;
  typedef TransformType::VersorType         VersorType;
  typedef VersorType::VectorType            VectorType;
  typedef MetricType::Pointer               MetricTypePointer;
  typedef OptimizerType::Pointer            OptimizerTypePointer;
  typedef OptimizerType::ParametersType     OptimizerParameterType;
  typedef OptimizerType::ScalesType         OptimizerScalesType;
  typedef InterpolatorType::Pointer         InterpolatorTypePointer;
  typedef RegistrationType::Pointer         RegistrationTypePointer;
  typedef TransformInitializerType::Pointer TransformInitializerTypePointer;

  /** Standard New method. */
  itkNewMacro(Self);

/** Runtime information support. */
  itkTypeMacro(AnatomicalVersorRigidFilter, itk::Object);

/* SetInput and GetOutput Macros */
  itkSetObjectMacro(FixedImage,  RegisterImageType);
  itkSetObjectMacro(MovingImage, RegisterImageType);
  itkGetConstObjectMacro(Output,      TransformType);

  itkSetMacro(NumberOfSpatialSamples, int);
  itkSetMacro(NumberOfIterations,     int);
  itkSetMacro(TranslationScale,       float);
  itkSetMacro(MaximumStepLength,      float);
  itkSetMacro(MinimumStepLength,      float);
  itkSetMacro(RelaxationFactor,       float);
  itkSetMacro(InitialRotationAngle,   double);
  itkSetMacro(InitialRotationAxis,    int);

  void Update();

protected:
  AnatomicalVersorRigidFilter();
  ~AnatomicalVersorRigidFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(AnatomicalVersorRigidFilter);

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
};  // end of class
} // end namespace itk

#endif
