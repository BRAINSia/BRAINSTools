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

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkOrientImageFilter.h"
#include "itkAnatomicalVersorRigidFilter.h"
#include <itkIOCommon.h>
#include <itkCastImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkExtractImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include <iostream>

namespace itk
{
AnatomicalVersorRigidFilter::AnatomicalVersorRigidFilter()
{
  m_InitialRotationAxis = 0;
  m_InitialRotationAngle = 0.0;
  m_NumberOfSpatialSamples = 100000;
  m_NumberOfIterations = 1000;
  m_TranslationScale = 1000.0;
  m_MaximumStepLength = 0.2;
  m_MinimumStepLength = 0.0001;
  m_RelaxationFactor = 0.5;
}

void AnatomicalVersorRigidFilter::Update()
{
  std::cout << "AnatomicalVersorRigidFilter()...." << std::endl;

  MetricTypePointer       metric        = MetricType::New();
  OptimizerTypePointer    optimizer     = OptimizerType::New();
  InterpolatorTypePointer interpolator  = InterpolatorType::New();
  RegistrationTypePointer registration  = RegistrationType::New();
  TransformType::Pointer  transform     = TransformType::New();

  /*** Set up the Registration ***/
  metric->SetNumberOfSpatialSamples( m_NumberOfSpatialSamples );
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  registration->SetTransform(     transform     );

  registration->SetFixedImage(   m_FixedImage   );
  registration->SetMovingImage(   m_MovingImage   );
  registration->SetFixedImageRegion( m_FixedImage->GetBufferedRegion() );

  TransformInitializerTypePointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  m_FixedImage );
  initializer->SetMovingImage( m_MovingImage );
  initializer->MomentsOn();
  initializer->InitializeTransform();

  std::cout << "Initializer, center: " << transform->GetCenter()
            << ", offset: " << transform->GetOffset()
            << "." << std::endl;

  VersorType rotation;
  VectorType axis;

  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 0.0;
  axis[m_InitialRotationAxis] = 1.0;

  const double pi = 3.14159265358979323846;
  const double inRadians = pi / 180.0;
  const double angle = m_InitialRotationAngle * inRadians;

  rotation.Set(  axis, angle  );
  transform->SetRotation( rotation );
  registration->SetInitialTransformParameters( transform->GetParameters() );

  const double translationScale = 1.0 / m_TranslationScale;

  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales( optimizerScales );

  optimizer->SetMaximumStepLength( m_MaximumStepLength );
  optimizer->SetMinimumStepLength( m_MinimumStepLength );

  optimizer->SetRelaxationFactor( m_RelaxationFactor );

  optimizer->SetNumberOfIterations( m_NumberOfIterations );

  std::cout << "Before Rigid Registration, center: " << transform->GetCenter()
            << ", offset: " << transform->GetOffset()
            << "." << std::endl;

  try
    {
    registration->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    }

  OptimizerParameterType finalParameters = registration->GetLastTransformParameters();

  const double       versorX              = finalParameters[0];
  const double       versorY              = finalParameters[1];
  const double       versorZ              = finalParameters[2];
  const double       finalTranslationX    = finalParameters[3];
  const double       finalTranslationY    = finalParameters[4];
  const double       finalTranslationZ    = finalParameters[5];
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double       bestValue = optimizer->GetValue();

  // Print out results
  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << versorX  << std::endl;
  std::cout << " versor Y      = " << versorY  << std::endl;
  std::cout << " versor Z      = " << versorZ  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  if( numberOfIterations == static_cast<unsigned int>( m_NumberOfIterations ) )
    {
    std::cout << std::endl << std::endl
              << "|##>> Iterations maxed out!  A solution was not found in the numberOfIterations permitted. "
              << std::endl
              << std::endl;
    }

  transform->SetParameters( finalParameters );

  std::cout << "After Rigid Registration, center: " << transform->GetCenter()
            << ", offset: " << transform->GetOffset()
            << "." << std::endl;

  TransformType::MatrixType matrix = transform->GetMatrix();
  TransformType::OffsetType offset = transform->GetOffset();

  std::cout << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << std::endl << offset << std::endl;

  m_Output = TransformType::New();
  m_Output->SetCenter( transform->GetCenter() );
  m_Output->SetParameters( transform->GetParameters() );
}
} // end namespace itk
