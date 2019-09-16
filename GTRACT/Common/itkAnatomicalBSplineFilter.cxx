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
#include "itkAnatomicalBSplineFilter.h"
#include <itkIOCommon.h>
#include <itkCastImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkExtractImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"
#include <iostream>

namespace itk
{
AnatomicalBSplineFilter::AnatomicalBSplineFilter()
{
  m_SpatialSampleScale = 100;
  m_MaximumNumberOfIterations = 500;
  m_MaximumNumberOfEvaluations = 500;
  m_MaximumNumberOfCorrections = 12;
  m_BSplineHistogramBins = 50;
  m_GridSize[0] = 12;
  m_GridSize[1] = 4;
  m_GridSize[2] = 12;
  m_GridBorderSize = 3;
  m_CostFunctionConvergenceFactor = 1e+7;
  m_ProjectedGradientTolerance = 1e-4;

  m_BoundTypeX = 0;
  m_BoundTypeY = 0;
  m_BoundTypeZ = 0;
  m_LowerBoundX = 0.0;
  m_LowerBoundY = 0.0;
  m_LowerBoundZ = 0.0;
  m_UpperBoundX = 0.0;
  m_UpperBoundY = 0.0;
  m_UpperBoundZ = 0.0;

  m_BulkTransform = nullptr;
  m_Output = TransformType::New();
}

void
AnatomicalBSplineFilter::Update()
{
  std::cout << "AnatomicalBSplineFilter()...." << std::endl;

  MetricTypePointer       metric = MetricType::New();
  OptimizerTypePointer    optimizer = OptimizerType::New();
  InterpolatorTypePointer interpolator = InterpolatorType::New();
  RegistrationTypePointer registration = RegistrationType::New();
  // INFO: Needed for ITKv4 registration registration->InPlaceOn();

  /*** Set up the Registration ***/
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);
  registration->SetTransform(m_Output);

  /*** Setup the Registration ***/
  registration->SetFixedImage(m_FixedImage);
  registration->SetMovingImage(m_MovingImage);

  RegisterImageRegionType fixedImageRegion = m_FixedImage->GetBufferedRegion();

  registration->SetFixedImageRegion(fixedImageRegion);

  /*** Setup the B-SPline Parameters ***/
  TransformRegionType bsplineRegion;
  TransformSizeType   gridSizeOnImage;
  TransformSizeType   gridBorderSize;
  TransformSizeType   totalGridSize;

  gridSizeOnImage.SetSize(m_GridSize.GetSize());
  gridBorderSize.Fill(m_GridBorderSize); // Border for spline order = 3 ( 1
                                         // lower, 2 upper )
  totalGridSize = gridSizeOnImage + gridBorderSize;

  // bsplineRegion.SetSize( totalGridSize );

  if (m_BulkTransform.IsNotNull())
  {
    std::cout << "Using Bulk Transform" << std::endl;
    m_Output->SetBulkTransform(m_BulkTransform);
  }

  const unsigned int numberOfParameters = m_Output->GetNumberOfParameters();

  TransformParametersType parameters(numberOfParameters);

  parameters.Fill(0.0);

  m_Output->SetParameters(parameters);

  registration->SetInitialTransformParameters(m_Output->GetParameters());

  OptimizerBoundSelectionType boundSelect(m_Output->GetNumberOfParameters());
  OptimizerBoundValueType     upperBound(m_Output->GetNumberOfParameters());
  OptimizerBoundValueType     lowerBound(m_Output->GetNumberOfParameters());
  /* Old Method - Unbounded Deformations in X,Y,Z */
  // boundSelect.Fill( 0 );
  // upperBound.Fill( 0.0 );
  // lowerBound.Fill( 0.0 );
  /* New Method - User Specifies the Displacement Bounds in X,Y,Z */
  /*    Default is the same as the Old method          */
  for (unsigned int i = 0; i < boundSelect.size(); i += 3)
  {
    boundSelect[i + 0] = m_BoundTypeX;
    boundSelect[i + 1] = m_BoundTypeY;
    boundSelect[i + 2] = m_BoundTypeZ;
    lowerBound[i + 0] = m_LowerBoundX;
    lowerBound[i + 1] = m_LowerBoundY;
    lowerBound[i + 2] = m_LowerBoundZ;
    upperBound[i + 0] = m_UpperBoundX;
    upperBound[i + 1] = m_UpperBoundY;
    upperBound[i + 2] = m_UpperBoundZ;
  }

  optimizer->SetBoundSelection(boundSelect);
  optimizer->SetUpperBound(upperBound);
  optimizer->SetLowerBound(lowerBound);

  optimizer->SetCostFunctionConvergenceFactor(m_CostFunctionConvergenceFactor);
  optimizer->SetProjectedGradientTolerance(m_ProjectedGradientTolerance);
  optimizer->SetMaximumNumberOfIterations(m_MaximumNumberOfIterations);
  optimizer->SetMaximumNumberOfEvaluations(m_MaximumNumberOfEvaluations);
  optimizer->SetMaximumNumberOfCorrections(m_MaximumNumberOfCorrections);

  metric->SetNumberOfHistogramBins(m_BSplineHistogramBins);

  /*** Make this a Parameter ***/
  const unsigned int numberOfSamples = fixedImageRegion.GetNumberOfPixels() / m_SpatialSampleScale;

  metric->SetNumberOfSpatialSamples(numberOfSamples);
  metric->ReinitializeSeed(76926294);

  // Add a time probe
  itk::TimeProbesCollectorBase collector;

  std::cout << std::endl << "Starting Registration" << std::endl;

  try
  {
    collector.Start("Registration");
    registration->Update();
    collector.Stop("Registration");
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return;
  }

  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

  collector.Report();

  /* This call is required to copy the parameters */
  m_Output->SetParametersByValue(finalParameters);

  std::cout << m_Output << std::endl;
}
} // end namespace itk
