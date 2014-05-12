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
#include "itkIO.h"
#include "BRAINSCommonLib.h"

#ifdef USE_DebugImageViewer
#include "DebugImageViewerClient.h"
#endif

#include "BRAINSFitUtils.h"
#include "itkEuler3DTransform.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkOtsuHistogramMatchingImageFilter.h"
#include <fstream>
#include "BRAINSFitHelperTemplate.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkLBFGSBOptimizerv4.h"

#include "itkLabelObject.h"
#include "itkStatisticsLabelObject.h"
#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkMacro.h"

namespace itk
{
inline
void
ValidateTransformRankOrdering(const std::vector<std::string> & transformType)
{
  // Need to review transform types in the transform type vector to ensure that
  // they are ordered appropriately
  // for processing.  I.e. that they go from low dimensional to high
  // dimensional, and throw an error if
  // B-Spline is before Rigid, or any other non-sensical ordering of the
  // transform types.
  // Rigid=1, ScaleVersor3D=2, ScaleSkewVersor3D=3, Affine=4, and (BSpline or SyN)=5
#define VTROExceptionMacroMacro()                                       \
  itkGenericExceptionMacro(<< "Ordering of transforms does not proceed from\n" \
                           << "smallest to largest.  Please review settings for transformType.\n" \
                           << "Rigid < ScaleVersor3D < ScaleSkewVersor3D < Affine < (BSpline | SyN)")

  if( transformType.empty() )
    {
    itkGenericExceptionMacro(<< " ERROR:  No transform type specified.");
    }
  unsigned int CurrentTransformRank = 0;
  for( unsigned int l = 0; l < transformType.size(); ++l )
    {
    if( transformType[l] == "Rigid" )
      {
      if( CurrentTransformRank <= 1 )
        {
        CurrentTransformRank = 1;
        }
      else
        {
        VTROExceptionMacroMacro();
        }
      }
    else if( transformType[l] == "ScaleVersor3D" )
      {
      if( CurrentTransformRank <= 2 )
        {
        CurrentTransformRank = 2;
        }
      else
        {
        VTROExceptionMacroMacro();
        }
      }
    else if( transformType[l] == "ScaleSkewVersor3D" )
      {
      if( CurrentTransformRank <= 3 )
        {
        CurrentTransformRank = 3;
        }
      else
        {
        VTROExceptionMacroMacro();
        }
      }
    else if( transformType[l] == "Affine" )
      {
      if( CurrentTransformRank <= 4 )
        {
        CurrentTransformRank = 4;
        }
      else
        {
        VTROExceptionMacroMacro();
        }
      }
    else if( transformType[l] == "BSpline" )
      {
      if( CurrentTransformRank <= 5 )
        {
        CurrentTransformRank = 5;
        }
      else
        {
        VTROExceptionMacroMacro();
        }
      }
    else if( transformType[l] == "SyN" )
      {
      if( CurrentTransformRank <= 5 )
        {
        CurrentTransformRank = 5;
        }
      else
        {
        VTROExceptionMacroMacro();
        }
      }

    else
      {
      std::cerr << " ERROR:  Invalid transform type specified for element " << l << " of --transformType: "
                << transformType[l] << std::endl;
      VTROExceptionMacroMacro();
      }
    }
}

template <class FixedImageType, class MovingImageType, class TransformType,
          class SpecificInitializerType, typename MetricType>
typename TransformType::Pointer
DoCenteredInitialization( typename FixedImageType::Pointer & orientedFixedVolume,
                          typename MovingImageType::Pointer & orientedMovingVolume,
                          ImageMaskPointer & fixedMask,             // NOTE:  This is both input and output
                                                                    // variable,  the Mask is updated by
                                                                    // this function
                          ImageMaskPointer & movingMask,            // NOTE:  This is both input and output
                                                                    // variable,  the Mask is updated by
                                                                    // this function
                          std::string & initializeTransformMode,
                          typename MetricType::Pointer & CostMetricObject )
{
  typename TransformType::Pointer initialITKTransform = TransformType::New();
  initialITKTransform->SetIdentity();

  if( initializeTransformMode == "useGeometryAlign" )
    {
    // useGeometryAlign assumes objects are center in field of view, with
    // different
    typedef itk::CenteredTransformInitializer<TransformType, FixedImageType,
                                              MovingImageType> OrdinaryInitializerType;
    typename OrdinaryInitializerType::Pointer CenteredInitializer =
      OrdinaryInitializerType::New();

    CenteredInitializer->SetFixedImage(orientedFixedVolume);
    CenteredInitializer->SetMovingImage(orientedMovingVolume);
    CenteredInitializer->SetTransform(initialITKTransform);
    CenteredInitializer->GeometryOn();              // Use the image spce center
    CenteredInitializer->InitializeTransform();
    }
  else if( initializeTransformMode == "useCenterOfROIAlign" )
    {
    if( movingMask.IsNull() || fixedMask.IsNull() )
      {
      itkGenericExceptionMacro(<< "FAILURE:  Improper mode for initializeTransformMode: "
                               << initializeTransformMode);
      }

    typedef typename itk::ImageMaskSpatialObject<FixedImageType::ImageDimension> CROIImageMaskSpatialObjectType;
    typedef itk::Image<unsigned char, 3>                                         CROIMaskImageType;
    typename MovingImageType::PointType movingCenter;
    typename FixedImageType::PointType fixedCenter;

    typename CROIImageMaskSpatialObjectType::Pointer movingImageMask(
      dynamic_cast<CROIImageMaskSpatialObjectType *>( movingMask.GetPointer() ) );
    typename CROIMaskImageType::Pointer tempOutputMovingVolumeROI =
      const_cast<CROIMaskImageType *>( movingImageMask->GetImage() );

    typename CROIImageMaskSpatialObjectType::Pointer fixedImageMask(
      dynamic_cast<CROIImageMaskSpatialObjectType *>( fixedMask.GetPointer() ) );
    typename CROIMaskImageType::Pointer tempOutputFixedVolumeROI =
      const_cast<CROIMaskImageType *>( fixedImageMask->GetImage() );
    typedef itk::CastImageFilter<CROIMaskImageType, FixedImageType>  MaskToFixedCastType;
    typedef itk::CastImageFilter<CROIMaskImageType, MovingImageType> MaskToMovingCastType;

    typename MaskToFixedCastType::Pointer mask2fixedCast = MaskToFixedCastType::New();
    typename MaskToMovingCastType::Pointer mask2movingCast = MaskToMovingCastType::New();

    mask2fixedCast->SetInput(tempOutputFixedVolumeROI);
    mask2fixedCast->Update();

    mask2movingCast->SetInput(tempOutputMovingVolumeROI);
    mask2movingCast->Update();

    typename SpecificInitializerType::Pointer CenteredInitializer =
      SpecificInitializerType::New();

    CenteredInitializer->SetFixedImage(mask2fixedCast->GetOutput() );
    CenteredInitializer->SetMovingImage(mask2movingCast->GetOutput() );
    CenteredInitializer->SetTransform(initialITKTransform);
    CenteredInitializer->MomentsOn();              // Use intensity center of
                                                   // mass

    CenteredInitializer->InitializeTransform();
    }
  else if( initializeTransformMode == "useCenterOfHeadAlign" )
    {
    typedef itk::Image<unsigned char, 3> CHMMaskImageType;
    typename MovingImageType::PointType movingCenter;
    typename FixedImageType::PointType fixedCenter;

    typedef typename itk::FindCenterOfBrainFilter<MovingImageType>
      MovingFindCenterFilter;
    typename MovingFindCenterFilter::Pointer movingFindCenter =
      MovingFindCenterFilter::New();
    movingFindCenter->SetInput(orientedMovingVolume);
    if( movingMask.IsNotNull() )
      {
      typename ImageMaskSpatialObjectType::Pointer movingImageMask(
        dynamic_cast<ImageMaskSpatialObjectType *>( movingMask.GetPointer() ) );
      typename CHMMaskImageType::Pointer tempOutputMovingVolumeROI =
        const_cast<CHMMaskImageType *>( movingImageMask->GetImage() );
      movingFindCenter->SetImageMask(tempOutputMovingVolumeROI);
      }
    movingFindCenter->Update();
    movingCenter = movingFindCenter->GetCenterOfBrain();
      {
      // convert mask image to mask
      typename ImageMaskSpatialObjectType::Pointer mask = ImageMaskSpatialObjectType::New();
      mask->SetImage( movingFindCenter->GetClippedImageMask() );

      mask->ComputeObjectToWorldTransform();
      typename SpatialObjectType::Pointer p = dynamic_cast<SpatialObjectType *>( mask.GetPointer() );
      if( p.IsNull() )
        {
        itkGenericExceptionMacro(<< "Can't convert mask pointer to SpatialObject");
        }
      movingMask = p;
      }

    typedef typename itk::FindCenterOfBrainFilter<FixedImageType>
      FixedFindCenterFilter;
    typename FixedFindCenterFilter::Pointer fixedFindCenter =
      FixedFindCenterFilter::New();
    fixedFindCenter->SetInput(orientedFixedVolume);
    if( fixedMask.IsNotNull() )
      {
      typename ImageMaskSpatialObjectType::Pointer fixedImageMask(
        dynamic_cast<ImageMaskSpatialObjectType *>( fixedMask.GetPointer() ) );
      typename CHMMaskImageType::Pointer tempOutputFixedVolumeROI =
        const_cast<CHMMaskImageType *>( fixedImageMask->GetImage() );
      fixedFindCenter->SetImageMask(tempOutputFixedVolumeROI);
      }
    fixedFindCenter->Update();
    fixedCenter = fixedFindCenter->GetCenterOfBrain();

      {
      // convert mask image to mask
      typename ImageMaskSpatialObjectType::Pointer mask = ImageMaskSpatialObjectType::New();
      mask->SetImage( fixedFindCenter->GetClippedImageMask() );

      mask->ComputeObjectToWorldTransform();
      typename SpatialObjectType::Pointer p = dynamic_cast<SpatialObjectType *>( mask.GetPointer() );
      if( p.IsNull() )
        {
        itkGenericExceptionMacro(<< "Can't convert mask pointer to SpatialObject");
        }
      fixedMask = p;
      }

    const double movingHeadScaleGuessRatio = 1;
    /*
      *
      *fixedFindCenter->GetHeadSizeEstimate()/movingFindCenter->GetHeadSizeEstimate();
      */

    typename TransformType::InputPointType rotationCenter;
    typename TransformType::OutputVectorType translationVector;
    itk::Vector<double, 3> scaleValue;
    for( unsigned int i = 0; i < Dimension; i++ )
      {
      rotationCenter[i]    = fixedCenter[i];
      translationVector[i] = movingCenter[i] - fixedCenter[i];
      scaleValue[i] = movingHeadScaleGuessRatio;
      }
    typedef itk::Euler3DTransform<double> EulerAngle3DTransformType;
    typename EulerAngle3DTransformType::Pointer bestEulerAngles3D = EulerAngle3DTransformType::New();
    bestEulerAngles3D->SetCenter(rotationCenter);
    bestEulerAngles3D->SetTranslation(translationVector);

    typedef itk::Euler3DTransform<double> EulerAngle3DTransformType;
    typename EulerAngle3DTransformType::Pointer currentEulerAngles3D = EulerAngle3DTransformType::New();
    currentEulerAngles3D->SetCenter(rotationCenter);
    currentEulerAngles3D->SetTranslation(translationVector);

    CostMetricObject->SetMovingTransform(currentEulerAngles3D);
    CostMetricObject->Initialize();
      {
      currentEulerAngles3D->SetRotation(0, 0, 0);
      // Initialize with current guess;
      //double max_cc = CostMetricObject->GetValue( currentEulerAngles3D->GetParameters() );
      double max_cc = CostMetricObject->GetValue();

      // rough search in neighborhood.
      const double one_degree = 1.0F * vnl_math::pi / 180.0F;
      const double HAStepSize = 3.0 * one_degree;
      const double PAStepSize = 3.0 * one_degree;
      // Quick search just needs to get an approximate angle correct.
        {
        const double HARange = 12.0;
        for( double HA = -HARange * one_degree; HA <= HARange * one_degree; HA += HAStepSize )
          {
          const double PARange = 12.0;
          for( double PA = -PARange * one_degree; PA <= PARange * one_degree; PA += PAStepSize )
            {
            currentEulerAngles3D->SetRotation(PA, 0, HA);
            //const double current_cc = CostMetricObject->GetValue( currentEulerAngles3D->GetParameters() );
            const double current_cc = CostMetricObject->GetValue();
            if( current_cc < max_cc )
              {
              max_cc = current_cc;
              bestEulerAngles3D->SetFixedParameters( currentEulerAngles3D->GetFixedParameters() );
              bestEulerAngles3D->SetParameters( currentEulerAngles3D->GetParameters() );
              }
            // #define DEBUGGING_PRINT_IMAGES
#ifdef DEBUGGING_PRINT_IMAGES
              {
              std::cout << "quick search "
                        << " HA= " << ( currentEulerAngles3D->GetParameters()[2] ) * 180.0 / vnl_math::pi
                        << " PA= " << ( currentEulerAngles3D->GetParameters()[0] ) * 180.0 / vnl_math::pi
                        << " cc="  <<  current_cc
                        << std::endl;
              }
            if( 0 )
              {
              typedef itk::ResampleImageFilter<FixedImageType, MovingImageType, double> ResampleFilterType;
              typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

              resampler->SetTransform(currentEulerAngles3D);
              resampler->SetInput(orientedMovingVolume);
              // Remember:  the Data is Moving's, the shape is Fixed's.
              resampler->SetOutputParametersFromImage(orientedFixedVolume);
              resampler->Update();            //  Explicit Update() required
              // here.
              typename FixedImageType::Pointer ResampledImage = resampler->GetOutput();

              typedef itk::CheckerBoardImageFilter<FixedImageType> Checkerfilter;
              typename Checkerfilter::Pointer checker = Checkerfilter::New();
              unsigned int array[3] = { 36, 36, 36 };

              checker->SetInput1(orientedFixedVolume);
              checker->SetInput2(ResampledImage);
              checker->SetCheckerPattern(array);
              try
                {
                checker->Update();
                }
              catch( itk::ExceptionObject & err )
                {
                std::cout << "Caught an ITK exception: " << std::endl;
                std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
                throw;
                }
              char filename[300];
              sprintf(filename, "%05.2f_%05.2f_%05.2f.nii.gz",
                      ( currentEulerAngles3D->GetParameters()[2] ) * 180 / vnl_math::pi,
                      ( currentEulerAngles3D->GetParameters()[0] ) * 180 / vnl_math::pi, current_cc);

                {
                typedef typename itk::ImageFileWriter<FixedImageType> WriterType;
                typename WriterType::Pointer writer = WriterType::New();
                writer->UseCompressionOn();
                writer->SetFileName(filename);
                writer->SetInput( checker->GetOutput() );
                try
                  {
                  writer->Update();
                  }
                catch( itk::ExceptionObject & err )
                  {
                  std::cout << "Exception Object caught: " << std::endl;
                  std::cout << err << std::endl;
                  throw;
                  }
                }
              }
#endif
            }
          }
        }
      // DEBUGGING_PRINT_IMAGES INFORMATION
#ifdef DEBUGGING_PRINT_IMAGES
        {
        std::cout << "FINAL: quick search "
                  << " HA= " << ( bestEulerAngles3D->GetParameters()[2] ) * 180.0 / vnl_math::pi
                  << " PA= " << ( bestEulerAngles3D->GetParameters()[0] ) * 180.0 / vnl_math::pi
                  << " cc="  <<  max_cc
                  << std::endl;
        }
#endif
      }
    typename VersorRigid3DTransformType::Pointer quickSetVersor = VersorRigid3DTransformType::New();
    quickSetVersor->SetCenter( bestEulerAngles3D->GetCenter() );
    quickSetVersor->SetTranslation( bestEulerAngles3D->GetTranslation() );
      {
      itk::Versor<double> localRotation;
      localRotation.Set( bestEulerAngles3D->GetMatrix() );
      quickSetVersor->SetRotation(localRotation);
      }
#ifdef DEBUGGING_PRINT_IMAGES
      {
      typedef itk::ResampleImageFilter<FixedImageType, MovingImageType, double> ResampleFilterType;
      typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

      resampler->SetTransform(quickSetVersor);
      resampler->SetInput(orientedMovingVolume);
      // Remember:  the Data is Moving's, the shape is Fixed's.
      resampler->SetOutputParametersFromImage(orientedFixedVolume);
      resampler->Update();              //  Explicit Update() required here.
      typename FixedImageType::Pointer ResampledImage = resampler->GetOutput();

      typedef itk::CheckerBoardImageFilter<FixedImageType> Checkerfilter;
      typename Checkerfilter::Pointer checker = Checkerfilter::New();
      unsigned int array[3] = { 18, 18, 18 };

      checker->SetInput1(orientedFixedVolume);
      checker->SetInput2(ResampledImage);
      checker->SetCheckerPattern(array);
      try
        {
        checker->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "Caught an ITK exception: " << std::endl;
        std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
        throw;
        }
      char filename[300];
      sprintf(filename, "FINAL_%05.2f_%05.2f_%05.2f.nii.gz",
              ( bestEulerAngles3D->GetParameters()[2] ) * 180 / vnl_math::pi,
              ( bestEulerAngles3D->GetParameters()[0] ) * 180 / vnl_math::pi, max_cc);

        {
        typedef typename itk::ImageFileWriter<FixedImageType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        wirter->UseCompressionOn();
        writer->SetFileName(filename);
        writer->SetInput(ResampledImage);
        try
          {
          writer->Update();
          }
        catch( itk::ExceptionObject & err )
          {
          std::cout << "Exception Object caught: " << std::endl;
          std::cout << err << std::endl;
          throw;
          }
        }
      }
#endif
    AssignRigid::AssignConvertedTransform( initialITKTransform, quickSetVersor.GetPointer() );
    }
  else if( initializeTransformMode == "useMomentsAlign" )
    {
    // useMomentsAlign assumes that the structures being registered have same
    // amount of mass approximately uniformly distributed.
    typename SpecificInitializerType::Pointer CenteredInitializer =
      SpecificInitializerType::New();

    CenteredInitializer->SetFixedImage(orientedFixedVolume);
    CenteredInitializer->SetMovingImage(orientedMovingVolume);
    CenteredInitializer->SetTransform(initialITKTransform);
    CenteredInitializer->MomentsOn();                    // Use intensity center of mass

    CenteredInitializer->InitializeTransform();
    }
  else                     // can't happen unless an unimplemented CLP option was added:
    {
    itkGenericExceptionMacro(<< "FAILURE:  Improper mode for initializeTransformMode: "
                             << initializeTransformMode);
    }
  std::cout << "Initializing transform with "  << initializeTransformMode
            << " to " << std::endl;
  std::cout << initialITKTransform << std::endl;
  std::cout << "===============================================" << std::endl;
  return initialITKTransform;
}

template <class FixedImageType, class MovingImageType>
BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::BRAINSFitHelperTemplate() :
  m_FixedVolume(NULL),
  m_MovingVolume(NULL),
  m_FixedBinaryVolume(NULL),
  m_MovingBinaryVolume(NULL),
  m_OutputFixedVolumeROI(""),
  m_OutputMovingVolumeROI(""),
  m_SamplingPercentage(1),
  m_NumberOfHistogramBins(50),
  m_HistogramMatch(false),
  m_RemoveIntensityOutliers(0.00),
  m_NumberOfMatchPoints(10),
  m_NumberOfIterations(1, 1500),
  m_MaximumStepLength(0.2),
  m_RelaxationFactor(0.5),
  m_TranslationScale(1000.0),
  m_ReproportionScale(1.0),
  m_SkewScale(1.0),
  m_BackgroundFillValue(0.0),
  m_TransformType(1, "Rigid"),
  m_InitializeTransformMode("Off"),
  m_MaskInferiorCutOffFromCenter(1000),
  m_SplineGridSize(3, 10),
  m_CostFunctionConvergenceFactor(1e+9),
  m_ProjectedGradientTolerance(1e-5),
  m_MaxBSplineDisplacement(0.0),
  m_ActualNumberOfIterations(0),
  m_PermittedNumberOfIterations(0),
  m_DebugLevel(0),
  m_CurrentGenericTransform(NULL),
  m_DisplayDeformedImage(false),
  m_PromptUserAfterDisplay(false),
  m_FinalMetricValue(0.0),
  m_ObserveIterations(true),
  m_CostMetricObject(NULL),
  m_UseROIBSpline(0),
  m_PermitParameterVariation(0),
  m_SamplingStrategy(AffineRegistrationType::NONE),
  m_DoBSplineRegByWarpedMovingImage(true),
  m_ForceMINumberOfThreads(-1)
{
  m_SplineGridSize[0] = 14;
  m_SplineGridSize[1] = 10;
  m_SplineGridSize[2] = 12;
}

template <class FixedImageType, class MovingImageType>
template <class TransformType, class OptimizerType, class MetricType>
void
BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::FitCommonCode(
  int numberOfIterations,
  typename CompositeTransformType::Pointer & initialITKTransform)
{
  // FitCommonCode
  typedef typename itk::MultiModal3DMutualRegistrationHelper
    <TransformType,
     OptimizerType,
     FixedImageType,
     MovingImageType,
     MetricType> MultiModal3DMutualRegistrationHelperType;

  typename MultiModal3DMutualRegistrationHelperType::Pointer
  appMutualRegistration = MultiModal3DMutualRegistrationHelperType::New();

  appMutualRegistration->SetNumberOfHistogramBins(m_NumberOfHistogramBins);
  appMutualRegistration->SetNumberOfIterations( numberOfIterations);
  appMutualRegistration->SetPermitParameterVariation( m_PermitParameterVariation );
  appMutualRegistration->SetSamplingPercentage(m_SamplingPercentage);
  appMutualRegistration->SetRelaxationFactor( m_RelaxationFactor );
  appMutualRegistration->SetMaximumStepLength( m_MaximumStepLength );
  appMutualRegistration->SetTranslationScale( m_TranslationScale );
  appMutualRegistration->SetReproportionScale( m_ReproportionScale );
  appMutualRegistration->SetSkewScale( m_SkewScale );

  // NOTE: binary masks are set for the cost metric object!!!
  appMutualRegistration->SetFixedImage(    m_FixedVolume    );
  appMutualRegistration->SetMovingImage(   m_MovingVolume   );
  appMutualRegistration->SetCostMetricObject( this->m_CostMetricObject );
  appMutualRegistration->SetForceMINumberOfThreads( this->m_ForceMINumberOfThreads );

  appMutualRegistration->SetBackgroundFillValue(   m_BackgroundFillValue   );

  appMutualRegistration->SetInitialTransform( initialITKTransform.GetPointer() );
  appMutualRegistration->SetDisplayDeformedImage(m_DisplayDeformedImage);
  appMutualRegistration->SetPromptUserAfterDisplay(m_PromptUserAfterDisplay);
  appMutualRegistration->SetObserveIterations(m_ObserveIterations);
  appMutualRegistration->SetSamplingStrategy(m_SamplingStrategy);
  /*
   *  At this point appMutualRegistration should be all set to make
   *  an itk pipeline class templated in TransformType etc.
   *  with all its inputs in place;
   */
  // initialize the interconnects between components
  appMutualRegistration->Initialize();

  typename CompositeTransformType::Pointer finalTransform;
  try
    {
    appMutualRegistration->Update();
    finalTransform = appMutualRegistration->GetTransform(); // finalTransform is a composite transform

    // Find the metric value (It is needed when logFileReport flag is ON).
    //this->m_FinalMetricValue = appMutualRegistration->GetFinalMetricValue();

    this->m_ActualNumberOfIterations = appMutualRegistration->GetActualNumberOfIterations();
    this->m_PermittedNumberOfIterations = numberOfIterations;
    }
  catch( itk::ExceptionObject& err )
    {
    // pass exception to caller
    itkGenericExceptionMacro(<< "Exception caught during registration: " << err);
    }

  // Put the result transform in the CurrentGenericTransform.
  // This composite transform will be used for initailization of the next stage.
  //
  // Note1: The "finalTransform" is a composite transform that includes the result of the current registratin stage.
  //
  // Note2: Since the "finalTransform" has already been affected by the previous initial transform,
  //        we can remove the old initial transform from the queue of the CurrentGenericTransform.
  //
  this->m_CurrentGenericTransform->ClearTransformQueue();

  typename CompositeTransformType::Pointer compToAdd;
  typename CompositeTransformType::ConstPointer compXfrm =
                            dynamic_cast<const CompositeTransformType *>( finalTransform.GetPointer() );
  if( compXfrm.IsNotNull() )
    {
    compToAdd = compXfrm->Clone();
    this->m_CurrentGenericTransform = compToAdd;
    }
  else
    {
    itkExceptionMacro(<< "The registration output composite transform is a NULL transform.");
    }
}

template <class FixedImageType, class MovingImageType>
void
BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::Update(void)
{
  if( this->m_DebugLevel > 3 )
    {
    this->PrintSelf(std::cout, 3);
    }

  std::vector<int> localNumberOfIterations( m_TransformType.size() );
  if( m_NumberOfIterations.size() != m_TransformType.size() )
    {
    if( m_NumberOfIterations.size() != 1 )
      {
      itkGenericExceptionMacro(<< "ERROR:  Wrong number of parameters for NumberOfIterations."
                               << " It either needs to be 1 or the same size as TransformType.")
      }
    for( unsigned int q = 0; q < m_TransformType.size(); ++q )
      {
      localNumberOfIterations[q] = m_NumberOfIterations[0];
      }
    }
  else
    {
    localNumberOfIterations = m_NumberOfIterations;
    }
  std::string localInitializeTransformMode(this->m_InitializeTransformMode);
  for( unsigned int currentTransformIndex = 0;
       currentTransformIndex < m_TransformType.size();
       currentTransformIndex++ )
    {
    const std::string currentTransformType(m_TransformType[currentTransformIndex]);
    std::cout << "TranformTypes: "
              << currentTransformType << "(" << currentTransformIndex + 1 << " of " << m_TransformType.size() << ")."
              << std::endl;
    std::cout << std::flush << std::endl;
    }

  // Initialize Transforms
  // Note that we don't want to estimate initialization if initial moving transform is provided already.
  if( m_CurrentGenericTransform.IsNull() && localInitializeTransformMode != "Off" )
      // Use CenteredVersorTranformInitializer
  {
  typedef VersorRigid3DTransformType TransformType;
  std::cout << "Initializing transform with " << localInitializeTransformMode << std::endl;
  typedef itk::CenteredVersorTransformInitializer<FixedImageType,
  MovingImageType> InitializerType;

  TransformType::Pointer initialITKTransform =
  DoCenteredInitialization<FixedImageType, MovingImageType,
  TransformType, InitializerType, MetricType>(
                                              m_FixedVolume,
                                              m_MovingVolume,
                                              m_FixedBinaryVolume,
                                              m_MovingBinaryVolume,
                                              localInitializeTransformMode, this->m_CostMetricObject);

  // The currentGenericTransform will be initialized by estimated initial transform.
  this->m_CurrentGenericTransform = CompositeTransformType::New();
  this->m_CurrentGenericTransform->AddTransform( initialITKTransform.GetPointer() );

  localInitializeTransformMode = "Off";  // Now reset to Off once initialization is done.

  // Now if necessary clip the images based on m_MaskInferiorCutOffFromCenter
  DoCenteredTransformMaskClipping<TransformType, FixedImageType::ImageDimension>(
                                                                                 m_FixedBinaryVolume,
                                                                                 m_MovingBinaryVolume,
                                                                                 initialITKTransform,
                                                                                 m_MaskInferiorCutOffFromCenter);

    {   // Write out some debugging information if requested
    if( ( !this->m_FixedBinaryVolume.IsNull() ) && ( m_OutputFixedVolumeROI != "" ) )
      {
      const MaskImageType::ConstPointer tempOutputFixedVolumeROI =
      ExtractConstPointerToImageMaskFromImageSpatialObject(this->m_FixedBinaryVolume.GetPointer() );
      itkUtil::WriteConstImage<MaskImageType>(tempOutputFixedVolumeROI, m_OutputFixedVolumeROI);
      }
    if( ( !this->m_MovingBinaryVolume.IsNull() ) && ( m_OutputMovingVolumeROI != "" ) )
      {
      const MaskImageType::ConstPointer tempOutputMovingVolumeROI =
      ExtractConstPointerToImageMaskFromImageSpatialObject(this->m_MovingBinaryVolume.GetPointer() );
      itkUtil::WriteConstImage<MaskImageType>(tempOutputMovingVolumeROI, m_OutputMovingVolumeROI);
      }
    }
  }

  for( unsigned int currentTransformIndex = 0;
       currentTransformIndex < m_TransformType.size();
       currentTransformIndex++ )
    {
    const std::string currentTransformType(m_TransformType[currentTransformIndex]);
    std::cout << "\n\n\n=============================== "
              << "ITKv4 Registration: Starting Transform Estimations for "
              << currentTransformType << "(" << currentTransformIndex + 1
              << " of " << m_TransformType.size() << ")."
              << "==============================="
              << std::endl;
    std::cout << std::flush << std::endl;
    //
    // Break into cases on TransformType:
    //
    if( currentTransformType == "Rigid" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef VersorRigid3DTransformType                          TransformType;
      typedef itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
      //
      // Process the initialITKTransform as VersorRigid3DTransform:
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();

      if( m_CurrentGenericTransform.IsNotNull() ) //When null, m_CurrentGenericTransform will be initialized by identity.
        {
        /* NOTE1: m_CurrentGenericTransform is a composite transform. When not null, it is initialized by:
                 {Intial Moving Transform
                        OR
                  Estimated initial transform indicated by initialTransformMode
                        OR
                  Output of previous level}

         * NOTE2: If the transform Type of the current level is linear,
                  the m_CurrentGenericTransform should contain only one transform as we collapse all linear transforms together.
         */
        if( m_CurrentGenericTransform->GetNumberOfTransforms() != 1 )
          {
          itkGenericExceptionMacro("Linear initial composite transform should have only one component \
                                   as all linaear transforms are collapsed together.");
          }
        const GenericTransformType::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
            dynamic_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( ( transformFileType == "ScaleVersor3DTransform" )
                  || ( transformFileType == "ScaleSkewVersor3DTransform" )
                  || ( transformFileType == "AffineTransform" ) )
            {
            // CONVERTING TO RIGID TRANSFORM TYPE from other type:
            std::cout << "WARNING:  Extracting Rigid component type from transform." << std::endl;
            VersorRigid3DTransformType::Pointer tempInitializerITKTransform =
                                                ComputeRigidTransformFromGeneric(currInitTransformFormGenericComposite.GetPointer() );

            AssignRigid::AssignConvertedTransform( initialITKTransform, tempInitializerITKTransform.GetPointer() );
            }
          else
            {
            std::cout
            <<
            "Unsupported initial transform file -- TransformBase first transform typestring, "
            << transformFileType
            << " not equal to required type VersorRigid3DTransform"
            << std::endl;
            return;
            }
          }
        catch( itk::ExceptionObject & excp )
          {
          std::cout << "[FAILED]" << std::endl;
          std::cerr
          << "Error while reading the m_CurrentGenericTransform" << std::endl;
          std::cerr << excp << std::endl;
          throw;
          }
        }
      else
        {
        m_CurrentGenericTransform = CompositeTransformType::New();
        }
      // replace the original initial transform with the extracted version.
      m_CurrentGenericTransform->ClearTransformQueue();
      m_CurrentGenericTransform->AddTransform( initialITKTransform );

      this->FitCommonCode<TransformType, OptimizerType, MetricType>
        (localNumberOfIterations[currentTransformIndex],
        this->m_CurrentGenericTransform);
      // NOW, after running the above function, the m_CurrentGenericTransform contains the integration of initial transform and rigid registration results.
      ///////////////////////
      }
    else if( currentTransformType == "ScaleVersor3D" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef ScaleVersor3DTransformType                          TransformType; // NumberOfEstimatedParameter = 9;
      typedef itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
      //
      // Process the initialITKTransform as ScaleVersor3DTransform:
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        /* NOTE1: m_CurrentGenericTransform is a composite transform. When not null, it is initialized by:
         {Intial Moving Transform
         OR
         Estimated initial transform indicated by initialTransformMode
         OR
         Output of previous level}

         * NOTE2: If the transform Type of the current level is linear,
         the m_CurrentGenericTransform should contain only one transform as we collapse all linear transforms together.
         */
        if( m_CurrentGenericTransform->GetNumberOfTransforms() != 1 )
          {
          itkGenericExceptionMacro("Linear initial composite transform should have only one component \
                                   as all linaear transforms are collapsed together.");
          }
        const GenericTransformType::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleVersor3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( ( transformFileType == "ScaleSkewVersor3DTransform" )
                   || ( transformFileType == "AffineTransform" ) )
            {
            // CONVERTING TO RIGID TRANSFORM TYPE from other type:
            // TODO: we should preserve the Scale components
            std::cout << "WARNING:  Extracting Rigid component type from transform." << std::endl;
            VersorRigid3DTransformType::Pointer tempInitializerITKTransform = ComputeRigidTransformFromGeneric(
                currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform( initialITKTransform, tempInitializerITKTransform.GetPointer() );
            }
          else
            {
            std::cout
              <<
              "Unsupported initial transform file -- TransformBase first transform typestring, "
              << transformFileType
              <<
              " not equal to required type VersorRigid3DTransform OR ScaleVersor3DTransform"
              << std::endl;
            return;
            }
          }
        catch( itk::ExceptionObject & excp )
          {
          std::cout << "[FAILED]" << std::endl;
          std::cerr
            << "Error while reading the m_CurrentGenericTransform"
            << std::endl;
          std::cerr << excp << std::endl;
          throw;
          }
        }
      else
        {
        m_CurrentGenericTransform = CompositeTransformType::New();
        }

      // replace the original initial transform with the above converted version.
      m_CurrentGenericTransform->ClearTransformQueue();
      m_CurrentGenericTransform->AddTransform( initialITKTransform );

      this->FitCommonCode<TransformType, OptimizerType, MetricType>
      (localNumberOfIterations[currentTransformIndex],
       this->m_CurrentGenericTransform);
      // NOW, after running the above function, the m_CurrentGenericTransform contains the integration of initial transform and ScaleVersor registration results.
      /////////////////////
      }
    else if( currentTransformType == "ScaleSkewVersor3D" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef ScaleSkewVersor3DTransformType                       TransformType;  // NumberOfEstimatedParameter = 15;
      typedef itk::RegularStepGradientDescentOptimizerv4<double>   OptimizerType;
      //
      // Process the initialITKTransform as ScaleSkewVersor3D:
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        /* NOTE1: m_CurrentGenericTransform is a composite transform. When not null, it is initialized by:
         {Intial Moving Transform
         OR
         Estimated initial transform indicated by initialTransformMode
         OR
         Output of previous level}

         * NOTE2: If the transform Type of the current level is linear,
         the m_CurrentGenericTransform should contain only one transform as we collapse all linear transforms together.
         */
        if( m_CurrentGenericTransform->GetNumberOfTransforms() != 1 )
          {
          itkGenericExceptionMacro("Linear initial composite transform should have only one component \
                                   as all linaear transforms are collapsed together.");
          }
        const GenericTransformType::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleVersor3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleSkewVersor3DTransform" )
            {
            const ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleSkewVersor3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( ( transformFileType == "AffineTransform" ) )
            {
            // CONVERTING TO RIGID TRANSFORM TYPE from other type:
            // TODO:  We should really preserve the Scale and Skew components
            std::cout << "WARNING:  Extracting Rigid component type from transform." << std::endl;
            VersorRigid3DTransformType::Pointer tempInitializerITKTransform = ComputeRigidTransformFromGeneric(
                currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform( initialITKTransform, tempInitializerITKTransform.GetPointer() );
            }
          else
            {
            std::cout << "Unsupported initial transform file -- TransformBase first transform typestring, "
                      << transformFileType
                      << " not equal to required type VersorRigid3DTransform "
                      << "OR ScaleVersor3DTransform OR ScaleSkewVersor3DTransform"
                      << std::endl;
            return;
            }
          }
        catch( itk::ExceptionObject & excp )
          {
          std::cout << "[FAILED]" << std::endl;
          std::cerr << "Error while reading the m_CurrentGenericTransform"
                    << std::endl
                    << excp << std::endl;
          throw;
          }
        }
      else
        {
        m_CurrentGenericTransform = CompositeTransformType::New();
        }
      // replace the original initial transform with the above converted version.
      m_CurrentGenericTransform->ClearTransformQueue();
      m_CurrentGenericTransform->AddTransform( initialITKTransform );

      this->FitCommonCode<TransformType, OptimizerType, MetricType>
      (localNumberOfIterations[currentTransformIndex],
       this->m_CurrentGenericTransform);
      // NOW, after running the above function, the m_CurrentGenericTransform contains the integration of initial transform and ScaleSkew registration results that is an "Affine" transform.
      /////////////////////
      }
    else if( currentTransformType == "Affine" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef itk::AffineTransform<double, Dimension>                      TransformType;
      typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double>  OptimizerType;
      //
      // Process the initialITKTransform
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        /* NOTE1: m_CurrentGenericTransform is a composite transform. When not null, it is initialized by:
         {Intial Moving Transform
         OR
         Estimated initial transform indicated by initialTransformMode
         OR
         Output of previous level}

         * NOTE2: If the transform Type of the current level is linear,
         the m_CurrentGenericTransform should contain only one transform as we collapse all linear transforms together.
         */
        if( m_CurrentGenericTransform->GetNumberOfTransforms() != 1 )
          {
          itkGenericExceptionMacro("Linear initial composite transform should have only one component \
                                   as all linaear transforms are collapsed together.");
          }
        const GenericTransformType::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleVersor3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleSkewVersor3DTransform" )
            {
            const ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleSkewVersor3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "AffineTransform" )
            {
            const typename AffineTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<AffineTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else              //  NO SUCH CASE!!
            {
            std::cout
              << "Unsupported initial transform file -- TransformBase first transform typestring, "
              << transformFileType
              << " not equal to any recognized type VersorRigid3DTransform OR "
              << "ScaleVersor3DTransform OR ScaleSkewVersor3DTransform OR AffineTransform"
              << std::endl;
            return;
            }
          }
        catch( itk::ExceptionObject & excp )
          {
          std::cout << "[FAILED]" << std::endl;
          std::cerr << "Error while reading the m_CurrentGenericTransform"
                    << std::endl << excp << std::endl;
          throw;
          }
        }
      else
        {
        m_CurrentGenericTransform = CompositeTransformType::New();
        }

      // replace the original initial transform with the above converted version.
      m_CurrentGenericTransform->ClearTransformQueue();
      m_CurrentGenericTransform->AddTransform( initialITKTransform );

      this->FitCommonCode<TransformType, OptimizerType, MetricType>
      (localNumberOfIterations[currentTransformIndex],
       this->m_CurrentGenericTransform);
      // NOW, after running the above function, the m_CurrentGenericTransform contains the integration of initial transform and Affine registration results.
      /////////////////////
      }
    else if( currentTransformType == "BSpline" )
      {
      const unsigned int SpaceDimension = 3;
      const unsigned int SplineOrder = 3;
      typedef itk::BSplineTransform<double, SpaceDimension, SplineOrder> BSplineTransformType;

      typedef itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, BSplineTransformType> BSplineRegistrationType;
      typename BSplineRegistrationType::Pointer bsplineRegistration = BSplineRegistrationType::New();

      typename BSplineTransformType::Pointer initialBSplineTransform = BSplineTransformType::New();
      // Initialize the BSpline transform
      BSplineTransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
      BSplineTransformType::MeshSizeType             meshSize;
      BSplineTransformType::OriginType               fixedOrigin;

      for( unsigned int i=0; i< SpaceDimension; i++ )
        {
        fixedOrigin[i] = m_FixedVolume->GetOrigin()[i];
        fixedPhysicalDimensions[i] = m_FixedVolume->GetSpacing()[i]
        * static_cast<double>( m_FixedVolume->GetLargestPossibleRegion().GetSize()[i] - 1 );
        meshSize[i] = m_SplineGridSize[i];
        }

      initialBSplineTransform->SetTransformDomainOrigin( fixedOrigin );
      initialBSplineTransform->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
      initialBSplineTransform->SetTransformDomainMeshSize( meshSize );
      initialBSplineTransform->SetTransformDomainDirection( m_FixedVolume->GetDirection() );

      typedef BSplineTransformType::ParametersType     ParametersType;
      const unsigned int numberOfParameters =  initialBSplineTransform->GetNumberOfParameters();
      ParametersType parameters( numberOfParameters );
      parameters.Fill( 0.0 );

      initialBSplineTransform->SetParameters( parameters );
      // std::cout << "Intial Parameters = " << std::endl
      //           << initialBSplineTransform->GetParameters() << std::endl;

      bsplineRegistration->InitializeOutputTransformFromReference( initialBSplineTransform );

      // TODO:  Expose these to the command line for consistancy.
      const int m_MaximumNumberOfIterations = 1500;
      const int m_MaximumNumberOfEvaluations = 900;
      const int m_MaximumNumberOfCorrections = 12;

      typedef typename itk::LBFGSBOptimizerv4                  LBFGSBOptimizerType;
      typedef typename LBFGSBOptimizerType::Pointer            LBFGSBOptimizerTypePointer;
      typedef typename LBFGSBOptimizerType::ParametersType     OptimizerParameterType;
      typedef typename LBFGSBOptimizerType::ScalesType         OptimizerScalesType;
      typedef typename LBFGSBOptimizerType::BoundSelectionType OptimizerBoundSelectionType;
      typedef typename LBFGSBOptimizerType::BoundValueType     OptimizerBoundValueType;

      LBFGSBOptimizerTypePointer      LBFGSBoptimizer     = LBFGSBOptimizerType::New();

      /**
       *
       * Set the boundary condition for each variable, where
       * select[i] = 0 if x[i] is unbounded,
       *           = 1 if x[i] has only a lower bound,
       *           = 2 if x[i] has both lower and upper bounds, and
       *           = 3 if x[1] has only an upper bound
       */
        // TODO:  For control points outside the fixed image mask, it might be good to
        // constrian
        // the parameters to something different than those control points inside the
        // fixed image mask.
      OptimizerBoundSelectionType boundSelect( initialBSplineTransform->GetNumberOfParameters() );
      if( vcl_abs(m_MaxBSplineDisplacement) < 1e-12 )
        {
        boundSelect.Fill(0);
        }
      else
        {
        boundSelect.Fill(2);
        }
      OptimizerBoundValueType upperBound( initialBSplineTransform->GetNumberOfParameters() );
      upperBound.Fill(m_MaxBSplineDisplacement);
      OptimizerBoundValueType lowerBound( initialBSplineTransform->GetNumberOfParameters() );
      lowerBound.Fill(-m_MaxBSplineDisplacement);

      LBFGSBoptimizer->SetBoundSelection(boundSelect);
      std::cout << "PRE : " << LBFGSBoptimizer->GetUpperBound().size() << " " << upperBound.size() << std::endl;
      LBFGSBoptimizer->SetUpperBound(upperBound);
      std::cout << "POST: " << LBFGSBoptimizer->GetUpperBound().size() << " " << upperBound.size() << std::endl;

      std::cout << "PRE : " << LBFGSBoptimizer->GetLowerBound().size() << " " << lowerBound.size() << std::endl;
      LBFGSBoptimizer->SetLowerBound(lowerBound);
      std::cout << "POST: " << LBFGSBoptimizer->GetLowerBound().size() << " " << lowerBound.size() << std::endl;

      LBFGSBoptimizer->SetCostFunctionConvergenceFactor(m_CostFunctionConvergenceFactor);
      LBFGSBoptimizer->SetGradientConvergenceTolerance(m_ProjectedGradientTolerance);
      LBFGSBoptimizer->SetMaximumNumberOfIterations(m_MaximumNumberOfIterations);
      LBFGSBoptimizer->SetMaximumNumberOfFunctionEvaluations(m_MaximumNumberOfEvaluations);
      LBFGSBoptimizer->SetMaximumNumberOfCorrections(m_MaximumNumberOfCorrections);

      bsplineRegistration->SetFixedImage( 0, m_FixedVolume );

      if( this->m_DoBSplineRegByWarpedMovingImage )
        {
        if( this->m_CurrentGenericTransform.IsNotNull() )
          {
          std::cout << "Moving image is warped by initial transform, before it is passed to the BSpline registration." << std::endl;
          typedef float                                                                     VectorComponentType;
          typedef itk::Vector<VectorComponentType, GenericTransformImageNS::SpaceDimension> VectorPixelType;
          typedef itk::Image<VectorPixelType,  GenericTransformImageNS::SpaceDimension>     DisplacementFieldType;
          typename MovingImageType::Pointer warpedMoving = GenericTransformImage<
                                                                  MovingImageType,
                                                                  MovingImageType,
                                                                  DisplacementFieldType>( this->m_MovingVolume,
                                                                                         this->m_FixedVolume,
                                                                                         this->m_CurrentGenericTransform.GetPointer(),
                                                                                         m_BackgroundFillValue,
                                                                                         "Linear",
                                                                                         false );

          bsplineRegistration->SetMovingImage( 0, warpedMoving );

          // Then warp the moving image mask as well if it is not null.
          typedef typename MetricType::MovingImageMaskType     MovingImageMaskType;
          typename MovingImageMaskType::ConstPointer movingMask =
          dynamic_cast<MovingImageMaskType const *>( this->m_CostMetricObject->GetMovingImageMask() );
          if( movingMask.IsNotNull() )
            {
            typename MaskImageType::Pointer warpedMovingMaskImage = GenericTransformImage<
                                                                      MaskImageType,
                                                                      MaskImageType,
                                                                      DisplacementFieldType>(
                                                                                             ExtractConstPointerToImageMaskFromImageSpatialObject(movingMask).GetPointer(),
                                                                                             this->m_FixedVolume,
                                                                                             this->m_CurrentGenericTransform.GetPointer(),
                                                                                             m_BackgroundFillValue,
                                                                                             "Linear",
                                                                                             false );


            typename MovingImageMaskType::ConstPointer warpedMask = ConvertMaskImageToSpatialMask( warpedMovingMaskImage.GetPointer() );
            m_CostMetricObject->SetMovingImageMask( warpedMask );
            }
          }
        else
          {
          bsplineRegistration->SetMovingImage( 0, m_MovingVolume );
          m_CurrentGenericTransform = CompositeTransformType::New();
          m_CurrentGenericTransform->ClearTransformQueue();
          if( this->m_DebugLevel > 4 )
            {
            this->m_CurrentGenericTransform->Print(std::cout);
            }
          }
        }
      else
        {
        bsplineRegistration->SetMovingImage( 0, m_MovingVolume );

        if( this->m_CurrentGenericTransform.IsNull() )
          {
          m_CurrentGenericTransform = CompositeTransformType::New();
          m_CurrentGenericTransform->ClearTransformQueue();

          if( this->m_DebugLevel > 4 )
            {
            this->m_CurrentGenericTransform->Print(std::cout);
            }
          }
        else
          {
          bsplineRegistration->SetMovingInitialTransform( this->m_CurrentGenericTransform );

          if( this->m_DebugLevel > 4 )
            {
            std::cout << "write the initial transform to the disk right before registration starts." << std::endl;
            this->m_CurrentGenericTransform->Print(std::cout);
            itk::TransformFileWriter::Pointer dwriter1 = itk::TransformFileWriter::New();
            dwriter1->SetInput( this->m_CurrentGenericTransform->GetNthTransform(0) );
            dwriter1->SetFileName("DEBUG_initial_transform_for_bspline.mat");
            dwriter1->Update();
            }
          }
        }

      // we have a 1 level BSpline registration.
      const unsigned int numberOfLevels = 1;

      typename BSplineRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
      shrinkFactorsPerLevel.SetSize( 1 );
      shrinkFactorsPerLevel[0] = 1;

      typename BSplineRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
      smoothingSigmasPerLevel.SetSize( 1 );
      smoothingSigmasPerLevel[0] = 0;

      bsplineRegistration->SetNumberOfLevels( numberOfLevels );
      bsplineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      bsplineRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      bsplineRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( true );
      bsplineRegistration->SetMetricSamplingStrategy(
                                                     static_cast<typename BSplineRegistrationType::MetricSamplingStrategyType>( m_SamplingStrategy ) );
      bsplineRegistration->SetMetricSamplingPercentage( m_SamplingPercentage );
      bsplineRegistration->SetMetric( this->m_CostMetricObject );
      bsplineRegistration->SetOptimizer( LBFGSBoptimizer );

      try
        {
        std::cout << "*** Running bspline registration (meshSizeAtBaseLevel = " << meshSize << ") ***"
        << std::endl << std::endl;
        bsplineRegistration->Update();
        }
      catch( itk::ExceptionObject & e )
        {
        itkGenericExceptionMacro( << "Exception caught: " << e << std::endl );
        }

      typename BSplineTransformType::Pointer outputBSplineTransform = BSplineTransformType::New();
      outputBSplineTransform->SetFixedParameters( bsplineRegistration->GetOutput()->Get()->GetFixedParameters() );
      outputBSplineTransform->SetParameters( bsplineRegistration->GetOutput()->Get()->GetParameters() );

      if( outputBSplineTransform.IsNotNull() )
        {
        this->m_CurrentGenericTransform->AddTransform( outputBSplineTransform );
        }
      else
        {
        itkGenericExceptionMacro( << "******* Error: the BSpline output transform is null." << std::endl );
        }
      const LBFGSBOptimizerType::ConstPointer  outputOptimizer =
                                                  dynamic_cast<const LBFGSBOptimizerType *>( bsplineRegistration->GetOptimizer() );
      std::cout << "Stop condition from LBFGSBoptimizer." << outputOptimizer->GetStopConditionDescription() << std::endl;
      std::cout << "Final Parameters = " << std::endl;
      std::cout << outputBSplineTransform->GetParameters() << std::endl;
      }
    else if( currentTransformType == "SyN" )
      {
#ifdef USE_ANTS
      //
      // Process SyN transform initializer
      //
      // current m_CurrentGenericTransform will be used as an initializer for SyN registration.
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        // Notice that we use only one SyN stage, and the outputs of all the previous linear levels
        // are composed in "one" transform.
        //
        if( m_CurrentGenericTransform->GetNumberOfTransforms() != 1 )
          {
          itkGenericExceptionMacro("Linear initial composite transform should have only one component \
                                   as all linaear transforms are collapsed together.");
          }
        const GenericTransformType::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();

        // Bspline transform cannot be used as an initializer for SyN registration.
        if( transformFileType == "BSplineTransform" )
          {
          itkGenericExceptionMacro( << "ERROR: Improper transform initializer for SyN registration: "
                                    << "BSpline Transform cannot be used as a transform initializer for SyN registration"
                                    << std::endl);
          }
        }
      else
        {
        // Initialize the registeration process with an Identity transform
        typename AffineTransformType::Pointer initialITKTransform = AffineTransformType::New();
        initialITKTransform->SetIdentity();

        m_CurrentGenericTransform = CompositeTransformType::New();
        m_CurrentGenericTransform->AddTransform( initialITKTransform );
        }

        typename CompositeTransformType::Pointer outputSyNTransform =
          simpleSynReg<FixedImageType, MovingImageType>( m_FixedVolume,
                                                        m_MovingVolume,
                                                        m_CurrentGenericTransform );

      if( outputSyNTransform.IsNull() )
        {
        std::cout << "\n*******Error: the SyN registration has failed.********\n" << std::endl;
        itkGenericExceptionMacro( << "******* Error: the SyN registration has failed." << std::endl );
        }
      else
        {
        // Update m_CurrentGenericTransform after SyN registration.
        m_CurrentGenericTransform = outputSyNTransform.GetPointer();
        }
#else
      std::cout << "******* Error: BRAINSFit cannot do the SyN registration ***"
                << "\n******* To use SyN option, you should also build ANTS ***" << std::endl;
      itkGenericExceptionMacro( << "******* Error: To use SyN registration, build ANTS too." << std::endl );
#endif
      }
    else
      {
      itkGenericExceptionMacro(
        << "Error choosing what kind of transform to fit \""
        << currentTransformType << "(" << currentTransformIndex + 1 << " of " << m_TransformType.size() << "). ");
      }
    }
  return;
}

template <class FixedImageType, class MovingImageType>
void
BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::PrintSelf(std::ostream & os, Indent indent) const
{
  // Superclass::PrintSelf(os,indent);
  os << indent << "FixedVolume:\n"  <<   this->m_FixedVolume << std::endl;
  os << indent << "MovingVolume:\n" <<   this->m_MovingVolume << std::endl;
  if( this->m_FixedBinaryVolume.IsNotNull() )
    {
    os << indent << "FixedBinaryVolume:\n" << this->m_FixedBinaryVolume << std::endl;
    }
  else
    {
    os << indent << "FixedBinaryVolume: IS NULL" << std::endl;
    }
  if( this->m_MovingBinaryVolume.IsNotNull() )
    {
    os << indent << "MovingBinaryVolume:\n" << this->m_MovingBinaryVolume << std::endl;
    }
  else
    {
    os << indent << "MovingBinaryVolume: IS NULL" << std::endl;
    }
  os << indent << "SamplingPercentage:      " << this->m_SamplingPercentage << std::endl;

  os << indent << "NumberOfIterations:    [";
  for( unsigned int q = 0; q < this->m_NumberOfIterations.size(); ++q )
    {
    os << this->m_NumberOfIterations[q] << " ";
    }
  os << "]" << std::endl;
  os << indent << "NumberOfHistogramBins:" << this->m_NumberOfHistogramBins << std::endl;
  os << indent << "MaximumStepLength:    " << this->m_MaximumStepLength << std::endl;
  os << indent << "TransformType:     [";
  for( unsigned int q = 0; q < this->m_TransformType.size(); ++q )
    {
    os << this->m_TransformType[q] << " ";
    }
  os << "]" << std::endl;

  os << indent << "RelaxationFactor:    " << this->m_RelaxationFactor << std::endl;
  os << indent << "TranslationScale:    " << this->m_TranslationScale << std::endl;
  os << indent << "ReproportionScale:   " << this->m_ReproportionScale << std::endl;
  os << indent << "SkewScale:           " << this->m_SkewScale << std::endl;
  os << indent << "BackgroundFillValue:            " << this->m_BackgroundFillValue << std::endl;
  os << indent << "InitializeTransformMode:        " << this->m_InitializeTransformMode << std::endl;
  os << indent << "MaskInferiorCutOffFromCenter:   " << this->m_MaskInferiorCutOffFromCenter << std::endl;
  os << indent << "ActualNumberOfIterations:       " << this->m_ActualNumberOfIterations << std::endl;
  os << indent << "PermittedNumberOfIterations:       " << this->m_PermittedNumberOfIterations << std::endl;

  os << indent << "SplineGridSize:     [";
  for( unsigned int q = 0; q < this->m_SplineGridSize.size(); ++q )
    {
    os << this->m_SplineGridSize[q] << " ";
    }
  os << "]" << std::endl;

  os << indent << "PermitParameterVariation:     [";
  for( unsigned int q = 0; q < this->m_PermitParameterVariation.size(); ++q )
    {
    os << this->m_PermitParameterVariation[q] << " ";
    }
  os << "]" << std::endl;

  if( m_CurrentGenericTransform.IsNotNull() )
    {
    os << indent << "CurrentGenericTransform:\n" << this->m_CurrentGenericTransform << std::endl;
    }
  else
    {
    os << indent << "CurrentGenericTransform: IS NULL" << std::endl;
    }
}

template <class FixedImageType, class MovingImageType>
void
BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::GenerateData()
{
  this->Update();
}
} // end namespace itk
