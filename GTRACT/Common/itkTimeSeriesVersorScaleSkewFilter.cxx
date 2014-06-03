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
#include "itkTimeSeriesVersorScaleSkewFilter.h"
#include <itkIOCommon.h>
#include <itkCastImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkExtractImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include <iostream>

namespace itk
{
TimeSeriesVersorScaleSkewFilter::TimeSeriesVersorScaleSkewFilter()
{
  m_NumberOfSpatialSamples = 100000;
  m_NumberOfIterations = 1000;
  m_NumberOfHistogramBins = 0;
  m_TranslationScale = 1000.0;
  m_ScalingScale = 1.0;
  m_SkewScale = 1.0;
  m_MaximumStepLength = 0.2;
  m_MinimumStepLength = 0.0001;
  m_RelaxationFactor = 0.5;
  m_BaseImage = 0;
  m_Output = InputImageType::New();
}

void TimeSeriesVersorScaleSkewFilter::Update()
{
  std::cout << "RigidRegisterDtiImages()...." << std::endl;

  MetricTypePointer       metric        = MetricType::New();
  OptimizerTypePointer    optimizer     = OptimizerType::New();
  InterpolatorTypePointer interpolator  = InterpolatorType::New();
  RegistrationTypePointer registration  = RegistrationType::New();
  registration->InPlaceOn();

  ExtractFilterTypePointer extractBaseImageFilter = ExtractFilterType::New();
  ExtractFilterTypePointer extractImageFilter     = ExtractFilterType::New();

  InputImageRegionType  fixedRegion  = m_Input->GetLargestPossibleRegion();
  InputImageSizeType    fixedSize    = fixedRegion.GetSize();
  InputImageSpacingType fixedSpacing = m_Input->GetSpacing();
  InputImagePointType   fixedOrigin  = m_Input->GetOrigin();
  int                   numVolumes = fixedSize[3];

  m_Output->SetRegions(fixedRegion);
  m_Output->SetSpacing(fixedSpacing);
  m_Output->SetOrigin(fixedOrigin);
  m_Output->Allocate();

  std::cout << "Region: " << fixedRegion << std::endl;
  std::cout << "Spacing: " << fixedSpacing << std::endl;
  std::cout << "Origin: " << fixedOrigin << std::endl;

  fixedSize[3] = m_BaseImage;
  fixedRegion.SetSize(fixedSize);

  InputImageIndexType fixedIndex = fixedRegion.GetIndex();
  fixedIndex[0] = 0;
  fixedIndex[1] = 0;
  fixedIndex[2] = 0;
  fixedIndex[3] = 0;
  fixedRegion.SetIndex(fixedIndex);
  std::cout << "Region: " << fixedRegion << std::endl;

  extractBaseImageFilter->SetExtractionRegion( fixedRegion );
  extractBaseImageFilter->SetInput(m_Input);
  extractBaseImageFilter->Update();
  for( int i = 0; i < numVolumes; i++ )
    {
    std::cout << "\tVolume: " << i << std::endl;

    /*** Set the Extraction Region ***/
    fixedIndex[3] = i;
    fixedRegion.SetIndex(fixedIndex);
    extractImageFilter->SetExtractionRegion( fixedRegion );
    extractImageFilter->SetInput(m_Input);

    /*** Set up the Registration ***/
    metric->SetNumberOfHistogramBins( m_NumberOfHistogramBins );
    metric->SetNumberOfSpatialSamples( m_NumberOfSpatialSamples );
    registration->SetMetric(        metric        );
    registration->SetOptimizer(     optimizer     );
    registration->SetInterpolator(  interpolator  );

    TransformType::Pointer transform = TransformType::New();
    registration->SetTransform( transform );

    extractImageFilter->Update();

    std::cout << extractBaseImageFilter->GetOutput()->GetBufferedRegion() << std::endl;
    std::cout << extractImageFilter->GetOutput()->GetBufferedRegion() << std::endl;

    itk::Point<double, 3> zeroOrigin;
    zeroOrigin.GetVnlVector().fill(0.0);
    extractBaseImageFilter->GetOutput()->SetOrigin(zeroOrigin);
    extractImageFilter->GetOutput()->SetOrigin(zeroOrigin);

    registration->SetFixedImage( extractBaseImageFilter->GetOutput() );
    registration->SetMovingImage( extractImageFilter->GetOutput() );
    registration->SetFixedImageRegion( extractBaseImageFilter->GetOutput()->GetBufferedRegion() );

    TransformInitializerTypePointer initializer = TransformInitializerType::New();
    initializer->SetTransform(   transform );
    initializer->SetFixedImage( extractBaseImageFilter->GetOutput() );
    initializer->SetMovingImage( extractImageFilter->GetOutput() );
    initializer->MomentsOn();
    initializer->InitializeTransform();

    std::cout << "Initializer, center: " << transform->GetCenter()
              << ", offset: " << transform->GetOffset()
              << "." << std::endl;

    VersorType rotation;
    VectorType axis;

    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;

    const double angle = 0;

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

    /*** Scaling Scales ***/
    const double scaleScale = 1.0 / m_ScalingScale;
    optimizerScales[6] = scaleScale;
    optimizerScales[7] = scaleScale;
    optimizerScales[8] = scaleScale;

    /*** Skew Scale ***/
    const double skewScale = 1.0 / m_SkewScale;
    optimizerScales[9] = skewScale;
    optimizerScales[10] = skewScale;
    optimizerScales[11] = skewScale;
    optimizerScales[12] = skewScale;
    optimizerScales[13] = skewScale;
    optimizerScales[14] = skewScale;

    optimizer->SetScales( optimizerScales );
    optimizer->SetMaximumStepLength( m_MaximumStepLength );
    optimizer->SetMinimumStepLength( m_MinimumStepLength );
    optimizer->SetRelaxationFactor( m_RelaxationFactor );
    optimizer->SetNumberOfIterations( m_NumberOfIterations );

    std::cout << "Before Scale-Skew Registration, center: " << transform->GetCenter()
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

    // const double versorX              = finalParameters[0];
    // const double versorY              = finalParameters[1];
    // const double versorZ              = finalParameters[2];
    // const double finalTranslationX    = finalParameters[3];
    // const double finalTranslationY    = finalParameters[4];
    // const double finalTranslationZ    = finalParameters[5];
    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    const double       bestValue = optimizer->GetValue();

    // Print out results
    std::cout << std::endl << std::endl;
    std::cout << "Result = " << std::endl;
    std::cout << " versor X      = " << finalParameters[0]  << std::endl;
    std::cout << " versor Y      = " << finalParameters[1]  << std::endl;
    std::cout << " versor Z      = " << finalParameters[2]  << std::endl;
    std::cout << " Translation X = " << finalParameters[3]  << std::endl;
    std::cout << " Translation Y = " << finalParameters[4]  << std::endl;
    std::cout << " Translation Z = " << finalParameters[5]  << std::endl;
    std::cout << " Scale X       = " << finalParameters[6]  << std::endl;
    std::cout << " Scale Y       = " << finalParameters[7]  << std::endl;
    std::cout << " Scale Z       = " << finalParameters[8]  << std::endl;
    std::cout << " Skew X        = " << finalParameters[9]  << std::endl;
    std::cout << " Skew Y        = " << finalParameters[10]  << std::endl;
    std::cout << " Skew Z        = " << finalParameters[11]  << std::endl;
    std::cout << " Skew XX       = " << finalParameters[12]  << std::endl;
    std::cout << " Skew YY       = " << finalParameters[13]  << std::endl;
    std::cout << " Skew ZZ       = " << finalParameters[14]  << std::endl;
    std::cout << " Iterations    = " << numberOfIterations << std::endl;
    std::cout << " Metric value  = " << bestValue          << std::endl;

    transform->SetParameters( finalParameters );

    std::cout << "After Scale-Skew Registration, center: " << transform->GetCenter()
              << ", offset: " << transform->GetOffset()
              << "." << std::endl;

    TransformType::MatrixType matrix = transform->GetMatrix();
    TransformType::OffsetType offset = transform->GetOffset();

    std::cout << "Matrix = " << std::endl << matrix << std::endl;
    std::cout << "Offset = " << std::endl << offset << std::endl;

    TransformTypePointer finalTransform = TransformType::New();
    finalTransform->SetCenter( transform->GetCenter() );
    finalTransform->SetParameters( transform->GetParameters() );
    /* Add Transform Writer */

    /* Resample the Image */
    ResampleFilterTypePointer resampler = ResampleFilterType::New();
    resampler->SetTransform( finalTransform );
    resampler->SetInput( extractImageFilter->GetOutput() );
    resampler->SetSize( extractBaseImageFilter->GetOutput()->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin( extractBaseImageFilter->GetOutput()->GetOrigin() );
    resampler->SetOutputSpacing( extractBaseImageFilter->GetOutput()->GetSpacing() );
    resampler->SetOutputDirection( extractBaseImageFilter->GetOutput()->GetDirection() );
    resampler->SetDefaultPixelValue( 0 );
    resampler->Update();

    ExtractImagePointer resampleImage = resampler->GetOutput();

    /*** Iterate and Update Image ***/
    std::cout << "Write Resampled Image" << std::endl;
    for( unsigned int z = 0; z < fixedSize[2]; z++ )
      {
      for( unsigned int y = 0; y < fixedSize[1]; y++ )
        {
        for( unsigned int x = 0; x < fixedSize[0]; x++ )
          {
          InputImageIndexType   index4D;
          ExtractImageIndexType index3D;
          index3D[0] = index4D[0] = x;
          index3D[1] = index4D[1] = y;
          index3D[2] = index4D[2] = z;
          index4D[3] = i;
          m_Output->SetPixel( index4D, resampleImage->GetPixel(index3D) );
          }
        }
      }
    }

  m_Output->SetMetaDataDictionary( m_Input->GetMetaDataDictionary() );
  m_Output->SetDirection( m_Input->GetDirection() );
  std::cout << "REGISTERED IMAGE: " << m_Output << std::endl;
}
} // end namespace itk
