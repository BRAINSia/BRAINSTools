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
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile$
 *  Language:  C++
 *  Date:      $Date: 2007-08-02 14:58:12 -0500 (Thu, 02 Aug 2007) $
 *  Version:   $Revision: 10282 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __genericRegistrationHelper_hxx
#define __genericRegistrationHelper_hxx

#include "genericRegistrationHelper.h"
#include "GenericTransformImage.h"

#include "itkVersor.h"
#include "itkMatrix.h"
#include "ConvertToRigidAffine.h"

#include "itkWindowedSincInterpolateImageFunction.h"

#include "itkStatisticsImageFilter.h"
#include "itkImageDuplicator.h"

#include "itkTransformFileWriter.h"

#include "itkImageMomentsCalculator.h"

extern void debug_catch(void);

namespace itk
{

/*
  * Constructor
  */
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::MultiModal3DMutualRegistrationHelper() :
  m_FixedImage(nullptr),                         // has to be provided by the user.
  m_MovingImage(nullptr),                        // has to be provided by the user.
  m_FixedImage2(nullptr),                        // has to be provided by the user.
  m_MovingImage2(nullptr),                       // has to be provided by the user.
  m_CompositeTransform(nullptr),              /* It is set by initial moving transform and
                                              integrates that with registration output transform.*/
  m_Transform(nullptr),                          // has to be provided by
                                           // this->Initialize().
  m_Registration(nullptr),                       // has to be provided by
                                           // this->Initialize().
  m_CostMetricObject(nullptr),
  m_SamplingPercentage(1),
  m_NumberOfHistogramBins(200),
  m_NumberOfIterations(0),
  m_RelaxationFactor(0.5),
  m_MaximumStepLength(0.2000),
  m_MinimumStepLength(0.0001),
  m_TranslationScale(1000.0),
  m_ReproportionScale(25.0),
  m_SkewScale(25.0),
  m_BackgroundFillValue(0.0),
  m_ActualNumberOfIterations(0),
  m_DisplayDeformedImage(false),
  m_PromptUserAfterDisplay(false),
  m_FinalMetricValue(0),
  m_ObserveIterations(true),
  m_SamplingStrategy(AffineRegistrationType::NONE),
  m_InternalTransformTime(0)
{
  this->SetNumberOfRequiredOutputs(1);    // for the Transform

  TransformOutputPointer transformDecorator = static_cast<TransformOutputType *>( this->MakeOutput(0U).GetPointer() );

  this->ProcessObject::SetNthOutput( 0, transformDecorator.GetPointer() );

  this->SetTransform( TransformType::New() );
  this->m_Transform->SetIdentity();
}

/*
 * Initialize by setting the interconnects between components.
 */
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::Initialize(void) // throw ( ExceptionObject )
{
  if( !m_FixedImage )
    {
    itkExceptionMacro(<< "FixedImage is not present");
    }

  if( !m_MovingImage )
    {
    itkExceptionMacro(<< "MovingImage is not present");
    }

  if( this->m_CompositeTransform.IsNull() )
    {
    itkExceptionMacro(<< "Input composite transform should include at least one Identity initial transform.");
    }

  m_Registration = RegistrationType::New();
  m_Registration->InPlaceOn();

  m_Transform = TransformType::New();

  if( this->m_CompositeTransform->GetNumberOfTransforms() == 1 )
    {
    const itk::Transform<double, 3, 3>::ConstPointer genericInit = this->m_CompositeTransform->GetFrontTransform();
    const typename TransformType::ConstPointer tempInitializerITKTransform =
                                      dynamic_cast<TransformType const *>( genericInit.GetPointer() );
    if( tempInitializerITKTransform.IsNull() )
      {
      std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
      }
    AssignRigid::AssignConvertedTransform(m_Transform, tempInitializerITKTransform);
    }
  else
    {
    itkExceptionMacro(<< "We should have just one initial transform in the input composite transform.");
    }

  //================================= SET SCALES AND OPTIMIZERS =============================================
  GenericOptimizerType::Pointer optimizer;

#if 0
  std::cout << "Transform Number of Parameters = " << m_Transform->GetNumberOfParameters() << std::endl;
#endif
  if( m_Transform->GetNumberOfParameters() == 12 )     //  Affine -> estimate scales automatically
    {
    typename MultiMetricType::Pointer multiMetric =
                  dynamic_cast<MultiMetricType *>( this->m_CostMetricObject.GetPointer() );
    typename ImageMetricType::Pointer firstMetricComponent =
                  dynamic_cast<ImageMetricType *>( multiMetric->GetMetricQueue()[0].GetPointer() );

    using ScalesEstimatorType = itk::RegistrationParameterScalesFromPhysicalShift<ImageMetricType>;
    typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
    scalesEstimator->SetMetric( firstMetricComponent );
    scalesEstimator->SetTransformForward( true );

    using ConjugateGradientDescentOptimizerType = itk::ConjugateGradientLineSearchOptimizerv4Template<double>;
    typename ConjugateGradientDescentOptimizerType::Pointer affineOptimizer = ConjugateGradientDescentOptimizerType::New();
    // Set the parameters of ConjugateGradient optimizer
    affineOptimizer->SetLowerLimit( 0 );
    affineOptimizer->SetUpperLimit( 2 );
    affineOptimizer->SetEpsilon( 0.2 );
    affineOptimizer->SetLearningRate( m_MaximumStepLength );
    affineOptimizer->SetMaximumStepSizeInPhysicalUnits(m_MaximumStepLength);
    affineOptimizer->SetNumberOfIterations(m_NumberOfIterations);
    const double convergenceThreshold = 1e-6;
    constexpr int convergenceWindowSize = 10;
    affineOptimizer->SetMinimumConvergenceValue( convergenceThreshold );
    affineOptimizer->SetConvergenceWindowSize( convergenceWindowSize );
    affineOptimizer->SetDoEstimateLearningRateAtEachIteration( true );
    affineOptimizer->SetDoEstimateLearningRateOnce( false );
    affineOptimizer->SetScalesEstimator( scalesEstimator );

    optimizer = affineOptimizer;
    }
  else // versor transforms -> we get scales from input
    {
    //
    // For versor transforms we need to initialize the center of rotation. It can be done by setting
    // Fixed parameters that are derived from the center of mass of the fixed image.
    //
#if 0
    std::cout << "Versor transform fixed parameters are set from fixed image's center of mass ... " << std::endl;
#endif
    using FixedImageCalculatorType = typename itk::ImageMomentsCalculator< FixedImageType >;
    typename FixedImageCalculatorType::Pointer fixedCalculator = FixedImageCalculatorType::New();
    fixedCalculator->SetImage( m_FixedImage );
    fixedCalculator->Compute();
    typename FixedImageCalculatorType::VectorType fixedCenter = fixedCalculator->GetCenterOfGravity();

    const unsigned int numberOfFixedParameters = m_Transform->GetFixedParameters().Size(); // =3
    typename TransformType::ParametersType fixedParameters( numberOfFixedParameters );
    for (unsigned int i = 0; i < numberOfFixedParameters; ++i)
      {
      fixedParameters[i] = fixedCenter[i];
      }
    m_Transform->SetFixedParameters( fixedParameters );
#if 0
    std::cout << "Versor Transform Fixed Parameters: " << fixedParameters << "." << std::endl;
#endif

    const double translationScale  = 1.0 / m_TranslationScale;
    const double reproportionScale = 1.0 / m_ReproportionScale;
    const double skewScale         = 1.0 / m_SkewScale;

    OptimizerScalesType optimizerScales( m_Transform->GetNumberOfParameters() );

    if( m_Transform->GetNumberOfParameters() == 15 )     //  ScaleSkewVersor3D
      {
      for( unsigned int i = 0; i < m_Transform->GetNumberOfParameters(); ++i )
        {
        optimizerScales[i] = 1.0;
        }
      for( unsigned int i = 3; i < 6; ++i )
        {
        optimizerScales[i] = translationScale;
        }
      for( unsigned int i = 6; i < 9; ++i )
        {
        optimizerScales[i] = reproportionScale;
        }
      for( unsigned int i = 9; i < 15; ++i )
        {
        optimizerScales[i] = skewScale;
        }
      }
    else if( m_Transform->GetNumberOfParameters() == 9 )    // ScaleVersor3D
      {
      for( unsigned int i = 0; i < 3; ++i )
        {
        optimizerScales[i] = 1.0;
        }
      for( unsigned int i = 3; i < 6; ++i )
        {
        optimizerScales[i] = translationScale;
        }
      for( unsigned int i = 6; i < 9; ++i )
        {
        optimizerScales[i] = reproportionScale;
        }
      }
    else if( m_Transform->GetNumberOfParameters() == 6 )     //  VersorRigid3D
      {
      for( unsigned int i = 0; i < 3; ++i )
        {
        optimizerScales[i] = 1.0;
        }
      for( unsigned int i = 3; i < 6; ++i )
        {
        optimizerScales[i] = translationScale;
        }
      }
    else
      { // we only support Affine, VersorRigid3D, ScaleVersor3D or ScaleSkewVersor3D.
      itkGenericExceptionMacro(<< "ERROR: The optimization transform does not have sufficient number of parameters.");
      }
      // end of versor scaling
#if 0
    std::cout << "Initializer, optimizerScales: " << optimizerScales << "." << std::endl;
#endif

    using VersorOptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
    typename VersorOptimizerType::Pointer versorOptimizer = VersorOptimizerType::New();

    versorOptimizer->SetScales( optimizerScales );
    versorOptimizer->SetNumberOfIterations( m_NumberOfIterations );
    versorOptimizer->SetLearningRate( m_MaximumStepLength );
    versorOptimizer->SetRelaxationFactor( m_RelaxationFactor );
    versorOptimizer->SetMinimumStepLength( m_MinimumStepLength );
    versorOptimizer->SetGradientMagnitudeTolerance( 1e-4 );
    versorOptimizer->SetReturnBestParametersAndValue(true);

    optimizer = versorOptimizer;
    }

  std::vector<FixedImagePointer>  preprocessedFixedImagesList;
  std::vector<MovingImagePointer> preprocessedMovingImagesList;
  preprocessedFixedImagesList.push_back( m_FixedImage );
  preprocessedMovingImagesList.push_back( m_MovingImage );
  if( m_FixedImage2 && m_MovingImage2 )
    {
    preprocessedFixedImagesList.push_back( m_FixedImage2 );
    preprocessedMovingImagesList.push_back( m_MovingImage2 );
    }
  for (unsigned int n=0; n<preprocessedFixedImagesList.size(); n++)
    {
    m_Registration->SetFixedImage( n, preprocessedFixedImagesList[n] );
    m_Registration->SetMovingImage( n, preprocessedMovingImagesList[n] );
    }

  m_Registration->SetInitialTransform(m_Transform);
  m_Registration->SetMetric(this->m_CostMetricObject);
  m_Registration->SetOptimizer(optimizer);

////////////////////////HARD CODED PART//////////////////////
  constexpr unsigned int numberOfLevels = 1;

  m_Registration->SetNumberOfLevels( numberOfLevels );

  std::vector<unsigned int>  factors( numberOfLevels );
  factors[0] = 1;
  using ShrinkFactorsPerDimensionContainerType = typename RegistrationType::ShrinkFactorsPerDimensionContainerType;
  std::vector<ShrinkFactorsPerDimensionContainerType> shrinkFactorsPerDimensionForAllLevels;
  for( unsigned int n = 0; n < numberOfLevels; n++ )
    {
    ShrinkFactorsPerDimensionContainerType shrinkFactorsPerDimension(3);
    shrinkFactorsPerDimension.Fill(0);
    for (unsigned int d = 0; d <  3; ++d) // here we set all dimensions have the same shrink factor
      {
      shrinkFactorsPerDimension[d] = factors[n];
      }
    shrinkFactorsPerDimensionForAllLevels.push_back( shrinkFactorsPerDimension );
#if 0
    std::cout << "  Shrink factors (level " << n+1 << " out of " << numberOfLevels << "): " << shrinkFactorsPerDimension << std::endl;
#endif
    }

   // Get smoothing sigmas
  std::vector<float> sigmas( numberOfLevels );
  sigmas[0] = 0;
  typename RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize( sigmas.size() );
  for( unsigned int n = 0; n < smoothingSigmasPerLevel.Size(); n++ )
    {
    smoothingSigmasPerLevel[n] = sigmas[n];
    }

  // set smoothing sigma and shrink factors
  for( unsigned int level = 0; level < numberOfLevels; ++level )
    {
    m_Registration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
    }
  m_Registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
  m_Registration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( true );
////////////////////////////////////////////////////////////

  m_Registration->SetMetricSamplingStrategy(
                static_cast<typename RegistrationType::MetricSamplingStrategyType>( m_SamplingStrategy ));
  m_Registration->SetMetricSamplingPercentage(this->m_SamplingPercentage);
  m_Registration->MetricSamplingReinitializeSeed(121212);

  // Create the Command observer and register it with the optimizer.
  // TODO:  make this output optional.
  //
  if( this->m_ObserveIterations == true )
    {
    using CommandIterationUpdateType = BRAINSFit::CommandIterationUpdate<TOptimizer, TTransformType, TMovingImage>;
    typename CommandIterationUpdateType::Pointer observer =
      CommandIterationUpdateType::New();
    observer->SetDisplayDeformedImage(m_DisplayDeformedImage);
    observer->SetPromptUserAfterDisplay(m_PromptUserAfterDisplay);
    observer->SetPrintParameters(true);
    observer->SetMovingImage(m_MovingImage);
    observer->SetFixedImage(m_FixedImage);
    observer->SetTransform(m_Transform);

    using MattesMutualInformationMetricType = itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType>;
    typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
      dynamic_cast<MattesMutualInformationMetricType *>(this->m_CostMetricObject.GetPointer() );
    if( test_MMICostMetric.IsNotNull() )
      {
      observer->SetMetricObject(test_MMICostMetric);
      }
    optimizer->AddObserver(itk::IterationEvent(), observer);
    }

  std::cout << std::flush;
}

/*
 * Starts the Registration Process
 */
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::Update(void)
{
  if( this->m_Transform.IsNotNull() )
    {
    // Update only if the input image has been modified
    const ModifiedTimeType t1 = this->GetPrimaryOutput()->GetPipelineMTime();
    const ModifiedTimeType t2 = this->m_Transform->GetMTime();
    const ModifiedTimeType t = ( t1 > t2 ? t1 : t2 );
    if( t == this->m_InternalTransformTime )
      {
      return; // No need to update
      }
    this->m_InternalTransformTime = t;
    }

  try
    {
    m_Registration->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    // Attempt to auto-recover if too many samples were requested.
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    // Pass exception to caller
    throw err;
    }

    OptimizerPointer optimizer =
      dynamic_cast<OptimizerPointer>( m_Registration->GetOptimizer() );
    if( optimizer == nullptr )
      {
      itkExceptionMacro(<< "Failed to convert pointer to Optimizer type");
      }
    std::cout << "Stop condition from optimizer." << optimizer->GetStopConditionDescription() << std::endl;
    m_FinalMetricValue = optimizer->GetValue();
    m_ActualNumberOfIterations = optimizer->GetCurrentIteration();
    {
      this->m_InternalTransformTime = this->m_Transform->GetMTime();
    }

  // The initial transform has already affected the registration result,
  // so it should be removed from the output transform queue here.
  // Now, the final composite transform should only include the registration results.
  this->m_CompositeTransform->ClearTransformQueue();
  this->m_CompositeTransform->AddTransform( this->m_Transform );

#if 0
  if ( true ) // add DebugLevel here.
    {
    std::cout << "Write the output transform of the registration filter ..." << std::endl;
    itk::TransformFileWriter::Pointer dwriter2 = itk::TransformFileWriter::New();
    dwriter2->SetInput( this->m_Transform );
    dwriter2->SetFileName("DEBUGTransform_RegFilterOutput.mat");
#if ITK_VERSION_MAJOR >= 5
    dwriter2->SetUseCompression(true);
#endif
    try
      {
      dwriter2->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Exception Object caught: " << std::endl;
      std::cerr << err << std::endl;
      throw;
      }
    }
  std::cout << "\n---After registration:---";
#endif
#if 0
  typename TransformType::MatrixType matrix = m_Transform->GetMatrix();
  typename TransformType::OffsetType offset = m_Transform->GetOffset();
  std::cout << std::endl << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << offset << std::endl;
  std::cout << "--------------" << std::endl << std::endl;
#endif
}

/*
 *
 */
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
ModifiedTimeType
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::GetMTime() const
{
  unsigned long mtime = Superclass::GetMTime();
  unsigned long m;

  // Some of the following should be removed once ivars are put in the
  // input and output lists

  if( m_Transform )
    {
    m = m_Transform->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if( m_Registration )
    {
    m = m_Registration->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if( m_CompositeTransform )
    {
    m = m_CompositeTransform->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if( m_FixedImage )
    {
    m = m_FixedImage->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if( m_MovingImage )
    {
    m = m_MovingImage->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  return mtime;
}

/*
  * Generate Data
  */
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::GenerateData()
{
  this->Update();
}

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::SetFixedImage(FixedImagePointer fixedImage)
{
  itkDebugMacro("setting Fixed Image to " << fixedImage);

  if( this->m_FixedImage.GetPointer() != fixedImage )
    {
    this->m_FixedImage = fixedImage;
    this->ProcessObject::SetNthInput(0, this->m_FixedImage);
    this->Modified();
    }
}

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::SetMovingImage(MovingImagePointer movingImage)
{
  itkDebugMacro("setting Moving Image to " << movingImage);

  if( this->m_MovingImage.GetPointer() != movingImage )
    {
    this->m_MovingImage = movingImage;
    this->ProcessObject::SetNthInput(1, this->m_MovingImage);
    this->Modified();
    }
}

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::SetFixedImage2(FixedImagePointer fixedImage2)
{
itkDebugMacro("setting Fixed Image to " << fixedImage2);

if( this->m_FixedImage2.GetPointer() != fixedImage2 )
  {
  this->m_FixedImage2 = fixedImage2;
  this->ProcessObject::SetNthInput(2, this->m_FixedImage2);
  this->Modified();
  }
}

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::SetMovingImage2(MovingImagePointer movingImage2)
{
itkDebugMacro("setting Moving Image to " << movingImage2);

if( this->m_MovingImage2.GetPointer() != movingImage2 )
  {
  this->m_MovingImage2 = movingImage2;
  this->ProcessObject::SetNthInput(3, this->m_MovingImage2);
  this->Modified();
  }
}

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::SetInitialTransform(typename CompositeTransformType::Pointer initialTransform)
{
  itkDebugMacro("setting Initial Transform to " << initialTransform);
  if( initialTransform.IsNull() )
    {
    itkGenericExceptionMacro("initial transform is null");
    }
  typename CompositeTransformType::Pointer compToAdd;
  compToAdd = initialTransform->Clone();
  this->m_CompositeTransform = compToAdd;
}

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
typename MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer,
                                              TFixedImage, TMovingImage, MetricType>::CompositeTransformType::Pointer
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::GetTransform(void)
{
  return this->m_CompositeTransform;
}

/*
 *  Get Output
 */

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
const typename MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer,
                                                    TFixedImage, TMovingImage, MetricType>::TransformOutputType *
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                       TMovingImage, MetricType>
::GetOutput() const
  {
  return static_cast<const TransformOutputType *>(
    this->ProcessObject::GetOutput(0) );
  }

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
DataObject::Pointer
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::MakeOutput(unsigned int output)
{
  switch( output )
    {
    case 0:
      {
      return static_cast<DataObject *>( TransformOutputType::New().GetPointer() );
      }
      break;
    default:
      itkExceptionMacro(
        "MakeOutput request for an output number larger than the expected number of outputs");
      return nullptr;
    }
}

/*
 * PrintSelf
 */
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
void
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Fixed Image: " << m_FixedImage.GetPointer() << std::endl;
  os << indent << "Moving Image: " << m_MovingImage.GetPointer() << std::endl;
  if( m_FixedImage2 && m_MovingImage2 )
    {
    os << indent << "Fixed Image2: " << m_FixedImage2.GetPointer() << std::endl;
    os << indent << "Moving Image2: " << m_MovingImage2.GetPointer() << std::endl;
    }
  else
    {
    os << indent << "Fixed Image2: IS NULL" << std::endl;
    os << indent << "Moving Image2: IS NULL" << std::endl;
    }
}
} // end namespace itk

#endif
