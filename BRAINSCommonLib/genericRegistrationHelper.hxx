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

#include "itkVersor.h"
#include "itkMatrix.h"
#include "ConvertToRigidAffine.h"

#include "itkWindowedSincInterpolateImageFunction.h"

#include "itkStatisticsImageFilter.h"
#include "itkImageDuplicator.h"

#include "itkTransformFileWriter.h"

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
  m_FixedImage(0),                         // has to be provided by the user.
  m_MovingImage(0),                        // has to be provided by the user.
  m_CompositeTransform(NULL),              /* It is set by initial moving transform and
                                              integrates that with registration output transform.*/
  m_Transform(0),                          // has to be provided by
                                           // this->Initialize().
  m_Registration(0),                       // has to be provided by
                                           // this->Initialize().
  m_PermitParameterVariation(0),
  m_CostMetricObject(NULL),
  m_SamplingPercentage(1),
  m_NumberOfHistogramBins(200),
  m_NumberOfIterations(0),
  m_RelaxationFactor(0.5),
  m_MaximumStepLength(0.2000),
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
  m_ForceMINumberOfThreads(-1),
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

  //
  // Connect the transform to the Decorator.
  //
  //TransformOutputType *transformOutput = static_cast<TransformOutputType *>( this->ProcessObject::GetOutput(0) );

  //transformOutput->Set( m_Transform.GetPointer() );

  if( this->m_CompositeTransform.IsNull() )
    {
    this->m_CompositeTransform = CompositeTransformType::New();
    }

  typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
  typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
  scalesEstimator->SetMetric( this->m_CostMetricObject );
  scalesEstimator->SetTransformForward( true );
  //scalesEstimator->DebugOn();

  //HACK: optimizer is hardcoded.
  //typename OptimizerType::Pointer optimizer = OptimizerType::New();
  typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> ConjugateGradientDescentOptimizerType;
  ConjugateGradientDescentOptimizerType::Pointer optimizer = ConjugateGradientDescentOptimizerType::New();//////?????

  // Added for v4 optimizer
  optimizer->SetLowerLimit( 0 );
  optimizer->SetUpperLimit( 2 );
  optimizer->SetEpsilon( 0.2 );
  optimizer->SetLearningRate( m_MaximumStepLength );
  optimizer->SetMaximumStepSizeInPhysicalUnits(m_MaximumStepLength);
  optimizer->SetNumberOfIterations(m_NumberOfIterations);
  optimizer->SetScalesEstimator( scalesEstimator );

  const double convergenceThreshold = 1e-6;
  const int convergenceWindowSize = 10;
  optimizer->SetMinimumConvergenceValue( convergenceThreshold );
  optimizer->SetConvergenceWindowSize( convergenceWindowSize );

  optimizer->SetDoEstimateLearningRateAtEachIteration( true );
  optimizer->SetDoEstimateLearningRateOnce( false );

  //optimizer->DebugOn();

  m_Registration = RegistrationType::New();

  /*
  // writing fixed and moving images right before passing to the registration filter
  typedef typename  itk::ImageFileWriter< TFixedImage  > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("fixedImageRightBeforeRegistrationStarts_DEBUG.nii.gz");
  writer->SetInput(m_FixedImage);
  writer->Update();

  typename WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName("movingImageRightBeforeRegistrationStarts_DEBUG.nii.gz");
  writer2->SetInput(m_MovingImage);
  writer2->Update();
  *//////////////////////////////////////////////////////////////

  //m_CostMetricObject->DebugOn();

  m_Registration->SetFixedImage(0, m_FixedImage);
  m_Registration->SetMovingImage(0, m_MovingImage);
  m_Registration->SetMetric(this->m_CostMetricObject);

////////////////////////HARD CODED PART//////////////////////
  const unsigned int numberOfLevels = 2;

  m_Registration->SetNumberOfLevels( numberOfLevels );

  typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactorsList;
//  shrinkFactorsList[0]=1;
//  shrinkFactorsList[1]=1;
//  shrinkFactorsList[2]=1;
//  m_Registration->SetShrinkFactorsPerDimension( 0, shrinkFactorsList );

  shrinkFactorsList[0]=3;
  shrinkFactorsList[1]=1;

  for( unsigned int level = 0; level < numberOfLevels; ++level )
    {
    m_Registration->SetShrinkFactorsPerDimension( level, shrinkFactorsList[level] );
    }

  typename RegistrationType::SmoothingSigmasArrayType smoothingSigma;
  smoothingSigma.SetSize(2);
  smoothingSigma[0]=2;
  smoothingSigma[1]=0;

  m_Registration->SetSmoothingSigmasPerLevel(smoothingSigma);
  m_Registration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(false);
////////////////////////////////////////////////////////////

  m_Registration->SetMetricSamplingStrategy(
                static_cast<typename RegistrationType::MetricSamplingStrategyType>( m_SamplingStrategy ));
  m_Registration->SetMetricSamplingPercentage(this->m_SamplingPercentage);

  m_Registration->SetOptimizer(optimizer);

  if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
    {
    //this->m_CompositeTransform->SetOnlyMostRecentTransformToOptimizeOn();

    /*
    ////////
    if(this->m_CompositeTransform->GetNumberOfTransforms() == 1)
      {
      std::cout << "write the initial transform to the disk right before registration starts." << std::endl;
      this->m_CompositeTransform->Print(std::cout);
      itk::TransformFileWriter::Pointer dwriter1 = itk::TransformFileWriter::New();
      dwriter1->SetInput( this->m_CompositeTransform->GetFrontTransform() );
      dwriter1->SetFileName("initialTransformRightBeforeRegistrationStarts_DEBUG.mat");
      dwriter1->Update();
      }
    *//////

    m_Registration->SetMovingInitialTransform( this->m_CompositeTransform );
    }

  // Create the Command observer and register it with the optimizer.
  // TODO:  make this output optional.
  //

  if( this->m_ObserveIterations == true )
    {
    typedef BRAINSFit::CommandIterationUpdate<TOptimizer, TTransformType, TMovingImage>
      CommandIterationUpdateType;
    typename CommandIterationUpdateType::Pointer observer =
      CommandIterationUpdateType::New();
    observer->SetDisplayDeformedImage(m_DisplayDeformedImage);
    observer->SetPromptUserAfterDisplay(m_PromptUserAfterDisplay);
    observer->SetPrintParameters(true);
    observer->SetMovingImage(m_MovingImage);
    observer->SetFixedImage(m_FixedImage);
    observer->SetTransform(m_Transform);

    typedef COMMON_MMI_METRIC_TYPE<FixedImageType, MovingImageType> MattesMutualInformationMetricType;
    typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
      dynamic_cast<MattesMutualInformationMetricType *>(this->m_CostMetricObject.GetPointer() );
    if( test_MMICostMetric.IsNotNull() )
      {
      observer->SetMetricObject(test_MMICostMetric);
      }
    optimizer->AddObserver(itk::IterationEvent(), observer);
    // std::cout << "Observer Configured" << std::endl;
    }
  else
    {
    // std::cout << "Skipping Observer Configuration" << std::endl;
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

  //m_Registration->DebugOn();

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

  std::cout << "METRIC       THREADS USED: " << this->m_CostMetricObject->GetNumberOfThreadsUsed()
  << " of " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads() <<  std::endl;
  std::cout << "REGISTRATION THREADS USED: " << this->m_Registration->GetNumberOfThreads()
  << " of " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads() <<  std::endl;

    OptimizerPointer optimizer =
      dynamic_cast<OptimizerPointer>( m_Registration->GetOptimizer() );
    if( optimizer == NULL )
      {
      itkExceptionMacro(<< "Failed to convert pointer to Optimizer type");
      }
    std::cout << "Stop condition from optimizer." << optimizer->GetStopConditionDescription() << std::endl;
    m_FinalMetricValue = optimizer->GetValue();
    m_ActualNumberOfIterations = optimizer->GetCurrentIteration();

    m_Transform->SetFixedParameters( m_Registration->GetOutput()->Get()->GetFixedParameters() );
    m_Transform->SetParameters( m_Registration->GetOutput()->Get()->GetParameters() );
    /*
    /////////
    std::cout << "\n----------------\n";
    std::cout << "Registration output transform:" << std::endl;
    m_Transform->Print(std::cout);
    /////////
    */
      {
      this->m_InternalTransformTime = this->m_Transform->GetMTime();
      }

  this->m_CompositeTransform->AddTransform( this->m_Transform ); // The registration output transform is added to the
                                                                 // initial transform, and the result composite transform
                                                                 // is returned as the final transform of this class.

  std::cout << "\n---After registration:---";
  typename TransformType::MatrixType matrix = m_Transform->GetMatrix();
  typename TransformType::OffsetType offset = m_Transform->GetOffset();
  std::cout << std::endl << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << offset << std::endl;
  std::cout << "--------------" << std::endl << std::endl;

}

/*
 *
 */
template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
unsigned long
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
/*
  typename CompositeTransformType::Pointer compToAdd;
  typename CompositeTransformType::ConstPointer compXfrm =
    dynamic_cast<const CompositeTransformType *>( initialTransform.GetPointer() );
  if( compXfrm.IsNotNull() )
    {
    compToAdd = compXfrm->Clone();
    this->m_CompositeTransform = compToAdd;
    }
  else
    {
    compToAdd = CompositeTransformType::New();
    typename TransformType::Pointer xfrm = initialTransform->Clone();
    compToAdd->AddTransform( xfrm );
    this->m_CompositeTransform = compToAdd;
    }
*/
}

template <typename TTransformType, typename TOptimizer, typename TFixedImage,
          typename TMovingImage, typename MetricType>
typename MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer,
                                              TFixedImage, TMovingImage, MetricType>::CompositeTransformType::Pointer
MultiModal3DMutualRegistrationHelper<TTransformType, TOptimizer, TFixedImage,
                                     TMovingImage, MetricType>
::GetTransform(void)
{
  //this->Update();
  return this->m_CompositeTransform;
  //return this->m_Transform;
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
      return 0;
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
}
} // end namespace itk

#endif
