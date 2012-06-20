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
#include <vnl/vnl_inverse.h>

#include "itkVectorImageRegisterVersorRigidFilter.h"

#include <iostream>

namespace itk
{
template <class TInputImage, class TOutputImage>
VectorImageRegisterVersorRigidFilter<TInputImage, TOutputImage>
::VectorImageRegisterVersorRigidFilter()
{
  m_NumberOfSpatialSamples = 100000;
  m_NumberOfIterations = 1000;
  m_TranslationScale = 1000.0;
  m_MaximumStepLength = 0.2;
  m_MinimumStepLength = 0.0001;
  m_RelaxationFactor = 0.5;
  m_RegisterB0Only = false;
  m_OutputParameterFile = "";
}

template <class TInputImage, class TOutputImage>
void VectorImageRegisterVersorRigidFilter<TInputImage, TOutputImage>
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

  /* Allocate the Writer - if requested */
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = NULL;
  if( m_OutputParameterFile.length() != 0 )
    {
    transformWriter =  TransformWriterType::New();
    transformWriter->SetFileName( m_OutputParameterFile.c_str() );
    }

  itkDebugMacro( "\tVector Length: " << this->GetInput()->GetVectorLength() );
  for( int i = 0; i < static_cast<int>( this->GetInput()->GetVectorLength() ); i++ )
    {
    itkDebugMacro("\tRegister Volume: " << i);

    // Get Current Gradient Direction
    vnl_vector<double> curGradientDirection(3);

    char tmpStr[64];
    sprintf(tmpStr, "DWMRI_gradient_%04d", i);
    std::string KeyString(tmpStr);

    std::string NrrdValue;

    itk::MetaDataDictionary inputMetaDataDictionary = this->GetInput()->GetMetaDataDictionary();
    itk::ExposeMetaData<std::string>(inputMetaDataDictionary, KeyString, NrrdValue);
    /* %lf is 'long float', i.e., double. */
    sscanf(
      NrrdValue.c_str(), " %lf %lf %lf", &curGradientDirection[0], &curGradientDirection[1], &curGradientDirection[2]);

    extractImageFilter->SetIndex( i );
    extractImageFilter->Update();

    bool registerVolumeFlag = true;
    if( ( m_RegisterB0Only )
        && ( ( curGradientDirection[0] != 0.0 ) || ( curGradientDirection[1] != 0.0 )
             && ( curGradientDirection[2] != 0.0 ) ) )
      {
      registerVolumeFlag = false;
      }

    if( registerVolumeFlag )
      {
      /*** Set up the Registration ***/
      metric->SetNumberOfSpatialSamples( m_NumberOfSpatialSamples );
      registration->SetMetric(        metric        );
      registration->SetOptimizer(     optimizer     );
      registration->SetInterpolator(  interpolator  );
      registration->SetTransform(     transform     );
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

      registration->StartRegistration();

      OptimizerParameterType finalParameters = registration->GetLastTransformParameters();

#ifndef NDEBUG
      const double       versorX              = finalParameters[0];
      const double       versorY              = finalParameters[1];
      const double       versorZ              = finalParameters[2];
      const double       finalTranslationX    = finalParameters[3];
      const double       finalTranslationY    = finalParameters[4];
      const double       finalTranslationZ    = finalParameters[5];
      const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
      const double       bestValue = optimizer->GetValue();

      itkDebugMacro("\tResult = ");
      itkDebugMacro("\t\tVersor X      = " << versorX);
      itkDebugMacro("\t\tVersor Y      = " << versorY);
      itkDebugMacro("\t\tVersor Z      = " << versorZ);
      itkDebugMacro("\t\tTranslation X = " << finalTranslationX);
      itkDebugMacro("\t\tTranslation Y = " << finalTranslationY);
      itkDebugMacro("\t\tTranslation Z = " << finalTranslationZ);
      itkDebugMacro("\t\tIterations    = " << numberOfIterations);
      itkDebugMacro("\t\tMetric value  = " << bestValue);
#endif  // NDEBUG

      transform->SetParameters( finalParameters );

      /* This step can probably be removed */
      TransformTypePointer finalTransform = TransformType::New();
      finalTransform->SetCenter( transform->GetCenter() );
      finalTransform->SetParameters( transform->GetParameters() );

      /* Add Transform Writer */
      if( m_OutputParameterFile.length() != 0 )
        {
        transformWriter->AddTransform( finalTransform );
        }

      /* Resample the Image */
      resampler->SetTransform( finalTransform );
      resampler->SetInput( extractImageFilter->GetOutput() );

      castImageFilter->SetInput( resampler->GetOutput() );

      /* Update the Gradient direction */
      curGradientDirection = finalTransform->GetMatrix().GetVnlMatrix() * curGradientDirection;
      }
    else
      {
      itkDebugMacro("\tIgnore Volume: " << i);
      castImageFilter->SetInput( extractImageFilter->GetOutput() );
      }

    /* Cast the Image */

    castImageFilter->Update();

    /* Insert the Registered Vector Index Image into the Output Vector Image */
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

    /* Set the gradient direction in the resulting image */
    sprintf(tmpStr,
            " %18.15lf %18.15lf %18.15lf",
            curGradientDirection[0],
            curGradientDirection[1],
            curGradientDirection[2]);
    NrrdValue = tmpStr;
    itk::EncapsulateMetaData<std::string>(m_Output->GetMetaDataDictionary(), KeyString, NrrdValue);

    /* Update Progress */
    this->UpdateProgress( (float)( i + 1 ) / (float)( this->GetInput()->GetVectorLength() ) );
    }

  if( m_OutputParameterFile.length() != 0 )
    {
    transformWriter->Update();
    }
}
} // end namespace itk
