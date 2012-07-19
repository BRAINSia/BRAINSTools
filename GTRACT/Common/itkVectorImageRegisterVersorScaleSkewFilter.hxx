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

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkIOCommon.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkCastImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkMetaDataObject.h>
#include <itkProgressAccumulator.h>

#include "itkVectorImageRegisterVersorScaleSkewFilter.h"

#include <iostream>

namespace itk
{
template <class TInputImage, class TOutputImage>
VectorImageRegisterVersorScaleSkewFilter<TInputImage, TOutputImage>
::VectorImageRegisterVersorScaleSkewFilter()
{
  m_NumberOfSpatialSamples = 100000;
  m_NumberOfIterations = 1000;
  m_NumberOfHistogramBins = 32;
  m_TranslationScale = 1000.0;
  m_ScalingScale = 1.0;
  m_SkewScale = 1.0;
  m_MaximumStepLength = 0.2;
  m_MinimumStepLength = 0.0001;
  m_RelaxationFactor = 0.5;
  m_OutputParameterFile = "";
}

template <class TInputImage, class TOutputImage>
void VectorImageRegisterVersorScaleSkewFilter<TInputImage, TOutputImage>
::GenerateData()
{
  itkDebugMacro("VectorImageRegisterVersorRigidFilter::GenerateData()");

  // Create a process accumulator for tracking the progress of this minipipeline
  typename ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  // No need to allocate the output vcl_since the minipipeline does it
  // this->AllocateOutputs();

  MetricTypePointer       metric        = MetricType::New();
  OptimizerTypePointer    optimizer     = OptimizerType::New();
  InterpolatorTypePointer interpolator  = InterpolatorType::New();
  RegistrationTypePointer registration  = RegistrationType::New();

  /* Allocate Output Image*/
  m_Output = OutputImageType::New();
  m_Output->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
  m_Output->SetSpacing( this->GetInput()->GetSpacing() );
  m_Output->SetOrigin( this->GetInput()->GetOrigin() );
  m_Output->SetDirection( this->GetInput()->GetDirection() );
  m_Output->SetVectorLength( this->GetInput()->GetVectorLength() );
  m_Output->SetMetaDataDictionary( this->GetInput()->GetMetaDataDictionary() );
  m_Output->Allocate();

  /* Create the Vector Index Extraction / Cast Filter */
  VectorIndexFilterPointer extractImageFilter = VectorIndexFilterType::New();
  extractImageFilter->SetInput( this->GetInput() );

  /* Create the Transform */
  TransformType::Pointer          transform = TransformType::New();
  TransformInitializerTypePointer initializer = TransformInitializerType::New();

  /* Create the Image Resampling Filter */
  ResampleFilterTypePointer resampler = ResampleFilterType::New();
  resampler->SetSize( m_FixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( m_FixedImage->GetOrigin() );
  resampler->SetOutputSpacing( m_FixedImage->GetSpacing() );
  resampler->SetOutputDirection( m_FixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 0 );

  /* Create the Cast Image Filter */
  CastFilterTypePointer castImageFilter = CastFilterType::New();
  for( int i = 0; i < this->GetInput()->GetVectorLength(); i++ )
    {
    itkDebugMacro("\tRegister Volume: " << i);

    extractImageFilter->SetIndex( i );
    extractImageFilter->Update();

    /*** Set up the Registration ***/
    metric->SetNumberOfSpatialSamples( m_NumberOfSpatialSamples );
    metric->SetNumberOfHistogramBins( m_NumberOfHistogramBins );
    registration->SetMetric(        metric        );
    registration->SetOptimizer(     optimizer     );
    registration->SetInterpolator(  interpolator  );
    registration->SetTransform( transform );

    /* Do we need this ???
    itk::Point<double, 3> zeroOrigin;
    zeroOrigin.GetVnlVector().fill(0.0);
    extractImageFilter->GetOutput()->SetOrigin(zeroOrigin);
    */

    registration->SetFixedImage(    m_FixedImage    );
    registration->SetMovingImage( extractImageFilter->GetOutput() );
    registration->SetFixedImageRegion( m_FixedImage->GetBufferedRegion() );

    initializer->SetTransform(   transform );
    initializer->SetFixedImage(  m_FixedImage );
    initializer->SetMovingImage( extractImageFilter->GetOutput() );
    initializer->MomentsOn();
    initializer->InitializeTransform();

    VersorType rotation;
    VectorType axis;

    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;

    const double angle = 0;

    rotation.Set(  axis, angle  );
    transform->SetRotation( rotation );
    registration->SetInitialTransformParameters( transform->GetParameters() );

    OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

    /*** Rotation Scales ***/
    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;

    /*** Translation Scales ***/
    const double translationScale = 1.0 / m_TranslationScale;
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

    registration->Update();

    OptimizerParameterType finalParameters = registration->GetLastTransformParameters();

    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    const double       bestValue = optimizer->GetValue();

    itkDebugMacro("\tResult");
    itkDebugMacro("\t\tVersor X      = " << finalParameters[0] );
    itkDebugMacro("\t\tVersor Y      = " << finalParameters[1] );
    itkDebugMacro("\t\tVersor Z      = " << finalParameters[2] );
    itkDebugMacro("\t\tTranslation X = " << finalParameters[3] );
    itkDebugMacro("\t\tTranslation Y = " << finalParameters[4] );
    itkDebugMacro("\t\tTranslation Z = " << finalParameters[5] );
    itkDebugMacro("\t\tScale X       = " << finalParameters[6] );
    itkDebugMacro("\t\tScale Y       = " << finalParameters[7] );
    itkDebugMacro("\t\tScale Z       = " << finalParameters[8] );
    itkDebugMacro("\t\tSkew X        = " << finalParameters[9] );
    itkDebugMacro("\t\tSkew Y        = " << finalParameters[10] );
    itkDebugMacro("\t\tSkew Z        = " << finalParameters[11] );
    itkDebugMacro("\t\tSkew XX       = " << finalParameters[12] );
    itkDebugMacro("\t\tSkew YY       = " << finalParameters[13] );
    itkDebugMacro("\t\tSkew ZZ       = " << finalParameters[14] );
    itkDebugMacro("\t\tIterations    = " << numberOfIterations);
    itkDebugMacro("\t\tMetric value  = " << bestValue         );

    transform->SetParameters( finalParameters );

    /* This step can probably be removed */
    TransformTypePointer finalTransform = TransformType::New();
    finalTransform->SetCenter( transform->GetCenter() );
    finalTransform->SetParameters( transform->GetParameters() );

    /* Add Transform Writer */
    if( m_OutputParameterFile.length() != 0 )
      {
      typedef itk::TransformFileWriter TransformWriterType;
      TransformWriterType::Pointer transformWriter =  TransformWriterType::New();
      transformWriter->SetFileName( m_OutputParameterFile.c_str() );
      transformWriter->SetInput( finalTransform );
      if( index == 0 )
        {
        transformWriter->SetAppendOff();
        }
      else
        {
        transformWriter->SetAppendOn();
        }
      transformWriter->Update();
      }

    /* Resample the Image */
    resampler->SetTransform( finalTransform );
    resampler->SetInput( extractImageFilter->GetOutput() );

    /* Cast the Image */
    castImageFilter->SetInput( resampler->GetOutput() );
    castImageFilter->Update();

    /* Insert the Registered Vector Index Image into the Output Vecot Image */
    typedef ImageRegionConstIterator<FixedImageType> ConstIteratorType;
    typedef ImageRegionIterator<OutputImageType>     IteratorType;

    ConstIteratorType it( castImageFilter->GetOutput(), castImageFilter->GetOutput()->GetRequestedRegion() );

    IteratorType ot( m_Output, m_Output->GetRequestedRegion() );

    OutputImagePixelType vectorImagePixel;
    for( ot.GoToBegin(), it.GoToBegin(); !ot.IsAtEnd(); ++ot, ++it )
      {
      vectorImagePixel = ot.Get();
      vectorImagePixel[i] = it.Value();
      ot.Set( vectorImagePixel );
      }

    /* Update Progress */
    this->UpdateProgress( (float)( i + 1 ) / (float)( this->GetInput()->GetVectorLength() ) );
    }
}
} // end namespace itk
