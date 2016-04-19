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
          class SpecificInitializerType, typename DoCenteredInitializationMetricType>
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
                          typename DoCenteredInitializationMetricType::Pointer & CostMetricObject )
{
  typedef itk::Image<unsigned char, 3>                               MaskImageType;
  typedef itk::ImageMaskSpatialObject<MaskImageType::ImageDimension> ImageMaskSpatialObjectType;

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
  else if( initializeTransformMode == "useCenterOfROIAlign"
           || initializeTransformMode == "useCenterOfHeadAlign")
    {

    typedef typename itk::ImageMaskSpatialObject<FixedImageType::ImageDimension> CROIImageMaskSpatialObjectType;
    typedef itk::Image<unsigned char, 3>                                         CROIMaskImageType;
    typename MovingImageType::PointType movingCenter;
    typename FixedImageType::PointType fixedCenter;

    if(initializeTransformMode == "useCenterOfROIAlign")
      {

      if( movingMask.IsNull() || fixedMask.IsNull() )
        {
        itkGenericExceptionMacro(<< "FAILURE:  Improper mode for initializeTransformMode: "
                                 << initializeTransformMode);
        }

      typedef itk::StatisticsLabelObject< unsigned char, 3 > LabelObjectType;
      typedef itk::LabelImageToStatisticsLabelMapFilter< MaskImageType, MaskImageType >
        LabelStatisticsFilterType;

      typename CROIImageMaskSpatialObjectType::Pointer movingImageMask(
        dynamic_cast<CROIImageMaskSpatialObjectType *>( movingMask.GetPointer() ) );
      typename CROIMaskImageType::Pointer tempOutputMovingVolumeROI =
        const_cast<CROIMaskImageType *>( movingImageMask->GetImage() );

      typename CROIImageMaskSpatialObjectType::Pointer fixedImageMask(
        dynamic_cast<CROIImageMaskSpatialObjectType *>( fixedMask.GetPointer() ) );
      typename CROIMaskImageType::Pointer tempOutputFixedVolumeROI =
        const_cast<CROIMaskImageType *>( fixedImageMask->GetImage() );

      LabelStatisticsFilterType::Pointer movingImageToLabel = LabelStatisticsFilterType::New();
      movingImageToLabel->SetInput( movingImageMask->GetImage() );
      movingImageToLabel->SetFeatureImage( movingImageMask->GetImage() );
      movingImageToLabel->SetComputePerimeter(false);
      movingImageToLabel->Update();

      LabelStatisticsFilterType::Pointer fixedImageToLabel = LabelStatisticsFilterType::New();
      fixedImageToLabel->SetInput( fixedImageMask->GetImage() );
      fixedImageToLabel->SetFeatureImage( fixedImageMask->GetImage() );
      fixedImageToLabel->SetComputePerimeter(false);
      fixedImageToLabel->Update();

      LabelObjectType *movingLabel = movingImageToLabel->GetOutput()->GetNthLabelObject(0);
      LabelObjectType *fixedLabel = fixedImageToLabel->GetOutput()->GetNthLabelObject(0);

      LabelObjectType::CentroidType movingCentroid = movingLabel->GetCentroid();
      LabelObjectType::CentroidType fixedCentroid = fixedLabel->GetCentroid();

      movingCenter[0] = movingCentroid[0];
      movingCenter[1] = movingCentroid[1];
      movingCenter[2] = movingCentroid[2];

      fixedCenter[0] = fixedCentroid[0];
      fixedCenter[1] = fixedCentroid[1];
      fixedCenter[2] = fixedCentroid[2];

      }
    else if( initializeTransformMode == "useCenterOfHeadAlign" )
      {
      typedef itk::Image<unsigned char, 3> CHMMaskImageType;

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
    typedef itk::VersorRigid3DTransform<double>              VersorRigid3DTransformType;
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
  m_FixedVolume(ITK_NULLPTR),
  m_FixedVolume2(ITK_NULLPTR),
  m_MovingVolume(ITK_NULLPTR),
  m_MovingVolume2(ITK_NULLPTR),
  m_FixedBinaryVolume(ITK_NULLPTR),
  m_MovingBinaryVolume(ITK_NULLPTR),
  m_OutputFixedVolumeROI(""),
  m_OutputMovingVolumeROI(""),
  m_SamplingPercentage(1),
  m_NumberOfHistogramBins(50),
  m_HistogramMatch(false),
  m_RemoveIntensityOutliers(0.00),
  m_NumberOfMatchPoints(10),
  m_NumberOfIterations(1, 1500),
  m_MaximumStepLength(0.2),
  m_MinimumStepLength(1, 0.005),
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
  m_CurrentGenericTransform(ITK_NULLPTR),
  m_RestoreState(ITK_NULLPTR),
  m_DisplayDeformedImage(false),
  m_PromptUserAfterDisplay(false),
  m_FinalMetricValue(0.0),
  m_ObserveIterations(true),
  m_CostMetricObject(ITK_NULLPTR),
  m_UseROIBSpline(0),
  m_SamplingStrategy(AffineRegistrationType::NONE),
  m_InitializeRegistrationByCurrentGenericTransform(true),
  m_MaximumNumberOfEvaluations(900),
  m_MaximumNumberOfCorrections(12),
  m_SyNMetricType(""),
  m_SaveState(""),
  m_SyNFull(true),
  m_ForceMINumberOfThreads(-1)
{
  m_SplineGridSize[0] = 14;
  m_SplineGridSize[1] = 10;
  m_SplineGridSize[2] = 12;
}

template <class FixedImageType, class MovingImageType>
template <class TransformType, class OptimizerType, class FitCommonCodeMetricType>
void
BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::FitCommonCode(
  int numberOfIterations,
  double minimumStepLength,
  typename CompositeTransformType::Pointer & initialITKTransform)
{
  // FitCommonCode
  typedef typename itk::MultiModal3DMutualRegistrationHelper
    <TransformType,
     OptimizerType,
     FixedImageType,
     MovingImageType,
     FitCommonCodeMetricType> MultiModal3DMutualRegistrationHelperType;

  typename MultiModal3DMutualRegistrationHelperType::Pointer
  appMutualRegistration = MultiModal3DMutualRegistrationHelperType::New();

  appMutualRegistration->SetNumberOfHistogramBins(m_NumberOfHistogramBins);
  appMutualRegistration->SetNumberOfIterations( numberOfIterations);
  appMutualRegistration->SetSamplingStrategy(m_SamplingStrategy);
  appMutualRegistration->SetSamplingPercentage(m_SamplingPercentage);

  appMutualRegistration->SetRelaxationFactor( m_RelaxationFactor );
  appMutualRegistration->SetMaximumStepLength( m_MaximumStepLength );
  appMutualRegistration->SetMinimumStepLength( minimumStepLength );
  appMutualRegistration->SetTranslationScale( m_TranslationScale );
  appMutualRegistration->SetReproportionScale( m_ReproportionScale );
  appMutualRegistration->SetSkewScale( m_SkewScale );

  // NOTE: binary masks are set for the cost metric object!!!
  appMutualRegistration->SetFixedImage(    m_FixedVolume    );
  appMutualRegistration->SetMovingImage(   m_MovingVolume   );
  if( m_FixedVolume2.IsNotNull() && m_MovingVolume2.IsNotNull() )
    {
    appMutualRegistration->SetFixedImage2(    m_FixedVolume2    );
    appMutualRegistration->SetMovingImage2(   m_MovingVolume2   );
    }
  appMutualRegistration->SetCostMetricObject( this->m_CostMetricObject );

  appMutualRegistration->SetBackgroundFillValue(   m_BackgroundFillValue   );

  appMutualRegistration->SetInitialTransform( initialITKTransform.GetPointer() );
  appMutualRegistration->SetDisplayDeformedImage(m_DisplayDeformedImage);
  appMutualRegistration->SetPromptUserAfterDisplay(m_PromptUserAfterDisplay);
  appMutualRegistration->SetObserveIterations(m_ObserveIterations);
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
  std::string transformFileType;
  if ( finalTransform.IsNotNull() )
    {
    transformFileType = finalTransform->GetNameOfClass();
    }
  if( transformFileType == "CompositeTransform" )
    {
    typename CompositeTransformType::ConstPointer compXfrm =
                              static_cast<const CompositeTransformType *>( finalTransform.GetPointer() );
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
  typedef itk::Image<unsigned char, 3>                               MaskImageType;
  typedef itk::ImageMaskSpatialObject<MaskImageType::ImageDimension> ImageMaskSpatialObjectType;
  typedef itk::VersorRigid3DTransform<double>                        VersorRigid3DTransformType;

  // The m_CostMetricObject should be a multi metric type.
  typename MultiMetricType::Pointer multiMetric =
                                      dynamic_cast<MultiMetricType *>( this->m_CostMetricObject.GetPointer() );
  if( multiMetric.IsNull() )
    {
    itkGenericExceptionMacro("Error in metric type conversion");
    }

  std::vector<FixedImagePointer> preprocessedFixedImagesList;
  std::vector<MovingImagePointer> preprocessedMovingImagesList;
  preprocessedFixedImagesList.push_back( m_FixedVolume );
  preprocessedMovingImagesList.push_back( m_MovingVolume );
  if( m_FixedVolume2.IsNotNull() && m_MovingVolume2.IsNotNull() )
    {
    preprocessedFixedImagesList.push_back( m_FixedVolume2 );
    preprocessedMovingImagesList.push_back( m_MovingVolume2 );
    }

  if( this->m_DebugLevel > 3 )
    {
    this->PrintSelf(std::cout, 3);
    }
  std::vector<double> localMinimumStepLength( m_TransformType.size() );
  if( m_MinimumStepLength.size() != m_TransformType.size() )
    {
    if( m_MinimumStepLength.size() != 1 )
      {
      itkGenericExceptionMacro(<< "ERROR:  Wrong number of parameters for MinimumStepLength."
                               << "It either needs to be 1 or the same size as TransformType.");
      }
    for( unsigned int q = 0; q < m_TransformType.size(); ++q )
      {
      localMinimumStepLength[q] = m_MinimumStepLength[0];
      }
    }
  else
    {
    localMinimumStepLength = m_MinimumStepLength;
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
  typedef itk::VersorRigid3DTransform<double> TransformType;
  std::cout << "Initializing transform with " << localInitializeTransformMode << std::endl;
  typedef itk::CenteredVersorTransformInitializer<FixedImageType,
  MovingImageType> InitializerType;

  TransformType::Pointer initialITKTransform =
  DoCenteredInitialization<FixedImageType,
                           MovingImageType,
                           TransformType,
                           InitializerType,
                           MultiMetricType>(
                                            m_FixedVolume,
                                            m_MovingVolume,
                                            m_FixedBinaryVolume,
                                            m_MovingBinaryVolume,
                                            localInitializeTransformMode,
                                            multiMetric );

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
      typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
      //
      // Process the initialITKTransform as VersorRigid3DTransform:
      //
      typedef VersorRigid3DTransformType TransformType;
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
        const itk::Transform<double, 3, 3>::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const itk::VersorRigid3DTransform<double>::ConstPointer tempInitializerITKTransform =
            static_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
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
         localMinimumStepLength[currentTransformIndex],
        this->m_CurrentGenericTransform);
      // NOW, after running the above function, the m_CurrentGenericTransform contains the integration of initial transform and rigid registration results.
      ///////////////////////
      }
    else if( currentTransformType == "ScaleVersor3D" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef itk::ScaleVersor3DTransform<double>                          TransformType; // NumberOfEstimatedParameter = 9;
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
        const itk::Transform<double, 3, 3>::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              static_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const itk::ScaleVersor3DTransform<double>::ConstPointer tempInitializerITKTransform =
              static_cast<itk::ScaleVersor3DTransform<double> const *>( currInitTransformFormGenericComposite.GetPointer() );
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
       localMinimumStepLength[currentTransformIndex],
       this->m_CurrentGenericTransform);
      // NOW, after running the above function, the m_CurrentGenericTransform contains the integration of initial transform and ScaleVersor registration results.
      /////////////////////
      }
    else if( currentTransformType == "ScaleSkewVersor3D" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef itk::ScaleSkewVersor3DTransform<double>                       TransformType;  // NumberOfEstimatedParameter = 15;
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
        const itk::Transform<double, 3, 3>::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              static_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const itk::ScaleVersor3DTransform<double>::ConstPointer tempInitializerITKTransform =
              static_cast<itk::ScaleVersor3DTransform<double> const *>( currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleSkewVersor3DTransform" )
            {
            const itk::ScaleSkewVersor3DTransform<double>::ConstPointer tempInitializerITKTransform =
              static_cast<itk::ScaleSkewVersor3DTransform<double> const *>( currInitTransformFormGenericComposite.GetPointer() );
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
       localMinimumStepLength[currentTransformIndex],
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
        const itk::Transform<double, 3, 3>::ConstPointer currInitTransformFormGenericComposite =
                                                  m_CurrentGenericTransform->GetFrontTransform();
        try
          {
          const std::string transformFileType = currInitTransformFormGenericComposite->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              static_cast<VersorRigid3DTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const itk::ScaleVersor3DTransform<double>::ConstPointer tempInitializerITKTransform =
              static_cast<itk::ScaleVersor3DTransform<double> const *>( currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "ScaleSkewVersor3DTransform" )
            {
            const itk::ScaleSkewVersor3DTransform<double>::ConstPointer tempInitializerITKTransform =
              static_cast<itk::ScaleSkewVersor3DTransform<double> const *>( currInitTransformFormGenericComposite.GetPointer() );
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "AffineTransform" )
            {
            const typename AffineTransformType::ConstPointer tempInitializerITKTransform =
              static_cast<AffineTransformType const *>( currInitTransformFormGenericComposite.GetPointer() );
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
       localMinimumStepLength[currentTransformIndex],
       this->m_CurrentGenericTransform);
      // NOW, after running the above function, the m_CurrentGenericTransform contains the integration of initial transform and Affine registration results.
      /////////////////////
      }
    else if( currentTransformType == "BSpline" )
      {
      const unsigned int SpaceDimension = 3;
      const unsigned int SplineOrder = 3;
      typedef itk::BSplineTransform<double, SpaceDimension, SplineOrder> BSplineTransformType;

      typedef itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType> BSplineRegistrationType;
      typename BSplineRegistrationType::Pointer bsplineRegistration = BSplineRegistrationType::New();


      typename FixedImageType::Pointer initializationImage = FixedImageType::New();
      if ( m_UseROIBSpline  )
        {
        ImageMaskSpatialObjectType::Pointer roiMask = ImageMaskSpatialObjectType::New();
        if( m_MovingBinaryVolume.GetPointer() != ITK_NULLPTR )
          {
          ImageMaskSpatialObjectType::Pointer movingImageMask =
          dynamic_cast<ImageMaskSpatialObjectType *>(m_MovingBinaryVolume.GetPointer() );

          typedef itk::ResampleImageFilter<MaskImageType, MaskImageType, double> ResampleFilterType;
          ResampleFilterType::Pointer resampler = ResampleFilterType::New();

          if( m_CurrentGenericTransform.IsNotNull() )
            {
              // resample the moving mask, if available
            resampler->SetTransform(m_CurrentGenericTransform);
            resampler->SetInput( movingImageMask->GetImage() );
            resampler->SetOutputParametersFromImage( m_FixedVolume );
            resampler->Update();
            }
          if( m_FixedBinaryVolume.GetPointer() != ITK_NULLPTR )
            {
            typedef itk::AddImageFilter<MaskImageType, MaskImageType> AddFilterType;
            ImageMaskSpatialObjectType::Pointer fixedImageMask =
            dynamic_cast<ImageMaskSpatialObjectType *>(m_FixedBinaryVolume.GetPointer() );
            AddFilterType::Pointer adder = AddFilterType::New();
            adder->SetInput1(fixedImageMask->GetImage() );
            adder->SetInput2(resampler->GetOutput() );
            adder->Update();
            roiMask->SetImage(adder->GetOutput() );
            }
          else
            {
            roiMask->SetImage(resampler->GetOutput() );
            }
          }
        else if( m_FixedBinaryVolume.GetPointer() != ITK_NULLPTR )
          {
          ImageMaskSpatialObjectType::Pointer fixedImageMask =
          dynamic_cast<ImageMaskSpatialObjectType *>(m_FixedBinaryVolume.GetPointer() );
          roiMask->SetImage(fixedImageMask->GetImage() );
          }

        else
          {
          itkGenericExceptionMacro( << "ERROR: ROIBSpline mode can only be used with ROI(s) specified!");
          return;
          }

        typename FixedImageType::PointType roiOriginPt;
        typename FixedImageType::IndexType roiOriginIdx;
        typename FixedImageType::Pointer    roiImage = FixedImageType::New();
        typename FixedImageType::RegionType roiRegion =
        roiMask->GetAxisAlignedBoundingBoxRegion();
        typename FixedImageType::SpacingType roiSpacing =
        m_FixedVolume->GetSpacing();

        roiOriginIdx.Fill(0);
        m_FixedVolume->TransformIndexToPhysicalPoint(roiRegion.GetIndex(), roiOriginPt);
        roiRegion.SetIndex(roiOriginIdx);
        roiImage->SetRegions(roiRegion);
        roiImage->Allocate();
        roiImage->FillBuffer(1.);
        roiImage->SetSpacing(roiSpacing);
        roiImage->SetOrigin(roiOriginPt);
        roiImage->SetDirection( m_FixedVolume->GetDirection() );

        initializationImage = roiImage;
        }
      else
        {
        initializationImage = this->m_FixedVolume;
        }

      typename BSplineTransformType::Pointer bsplineTx =
                                                BSplineTransformType::New();
      // Initialize the BSpline transform
      // Using BSplineTransformInitializer
      //
      BSplineTransformType::MeshSizeType    meshSize;
      for( unsigned int i=0; i< SpaceDimension; i++ )
        {
        meshSize[i] = m_SplineGridSize[i];
        }

      typedef itk::BSplineTransformInitializer< BSplineTransformType,
                                                FixedImageType>         InitializerType;
      typename InitializerType::Pointer transformInitializer = InitializerType::New();

      transformInitializer->SetTransform( bsplineTx );
      transformInitializer->SetImage( initializationImage );
      transformInitializer->SetTransformDomainMeshSize( meshSize );
      transformInitializer->InitializeTransform();

      bsplineTx->SetIdentity();

      const int psize = bsplineTx->GetNumberOfParameters();
      std::cout << "Initialized BSpline transform is set to be an identity transform." << std::endl;
      std::cout << "  - Number of parameters = " << psize << std::endl;
      if( psize/15 > 1 )
        {
        std::cout << "-- WARNING: Only one in every " << (psize/15) << " parameters is printed on screen.\n" << std::endl;
        }
      //std::cout << "Intial Parameters = " << std::endl
      //          << bsplineTx->GetParameters() << std::endl;

      bsplineRegistration->SetInitialTransform( bsplineTx );
      bsplineRegistration->InPlaceOn(); // So bsplineTx is also the output transform
                                        //  of the registration filter.

      typedef typename itk::LBFGSBOptimizerv4                  LBFGSBOptimizerType;
      typedef typename LBFGSBOptimizerType::Pointer            LBFGSBOptimizerTypePointer;
      typedef typename LBFGSBOptimizerType::BoundSelectionType OptimizerBoundSelectionType;
      typedef typename LBFGSBOptimizerType::BoundValueType     OptimizerBoundValueType;

      LBFGSBOptimizerTypePointer      LBFGSBoptimizer     = LBFGSBOptimizerType::New();

      // TODO:  For control points outside the fixed image mask, it might be good to
      // constrian
      // the parameters to something different than those control points inside the
      // fixed image mask.

      OptimizerBoundSelectionType boundSelect( bsplineTx->GetNumberOfParameters() );
      if( std::abs(m_MaxBSplineDisplacement) < 1e-12 )
        {
        boundSelect.Fill( LBFGSBOptimizerType::UNBOUNDED );
        }
      else
        {
        boundSelect.Fill( LBFGSBOptimizerType::BOTHBOUNDED );
        }
      OptimizerBoundValueType upperBound( bsplineTx->GetNumberOfParameters() );
      upperBound.Fill(m_MaxBSplineDisplacement);
      OptimizerBoundValueType lowerBound( bsplineTx->GetNumberOfParameters() );
      lowerBound.Fill(-m_MaxBSplineDisplacement);

      LBFGSBoptimizer->SetBoundSelection(boundSelect);
      //std::cout << "PRE : " << LBFGSBoptimizer->GetUpperBound().size() << " " << upperBound.size() << std::endl;
      LBFGSBoptimizer->SetUpperBound(upperBound);
      //std::cout << "POST: " << LBFGSBoptimizer->GetUpperBound().size() << " " << upperBound.size() << std::endl;

      //std::cout << "PRE : " << LBFGSBoptimizer->GetLowerBound().size() << " " << lowerBound.size() << std::endl;
      LBFGSBoptimizer->SetLowerBound(lowerBound);
      //std::cout << "POST: " << LBFGSBoptimizer->GetLowerBound().size() << " " << lowerBound.size() << std::endl;

      LBFGSBoptimizer->SetCostFunctionConvergenceFactor( m_CostFunctionConvergenceFactor );
      LBFGSBoptimizer->SetGradientConvergenceTolerance( m_ProjectedGradientTolerance );
      LBFGSBoptimizer->SetNumberOfIterations( localNumberOfIterations[currentTransformIndex] );
      LBFGSBoptimizer->SetMaximumNumberOfFunctionEvaluations( m_MaximumNumberOfEvaluations );
      LBFGSBoptimizer->SetMaximumNumberOfCorrections( m_MaximumNumberOfCorrections );

      std::cout << "LBFGSB optimizer is used for BSpline registration using following parameters set:" << std::endl;
      std::cout << "-----------------------------------------------------------------------------------" << std::endl;
      std::cout << "NOTICE: You can use commandline options to adjust these parameters to find a\n "
                << "        probable better compromise between running time and registration precision." << std::endl;
      std::cout << "-----------------------------------------------------------------------------------" << std::endl;
      std::cout << " - Cost Function Convergence Factor : " << m_CostFunctionConvergenceFactor << std::endl;
      std::cout << " - Projected Gradient Tolerance     : " << m_ProjectedGradientTolerance << std::endl;
      std::cout << " - Maximum Number of Corrections    : " << m_MaximumNumberOfCorrections << std::endl;
      std::cout << " - Maximum Number of Evaluations    : " << m_MaximumNumberOfEvaluations << std::endl;
      std::cout << " - Maximum Number of Iterations     : " << localNumberOfIterations[currentTransformIndex] << std::endl << std::endl;

      typedef itk::Image<float, 3> RegisterImageType;

      // TODO: pass this option by commandline
      const bool ObserveIterations = true;
      if( ObserveIterations == true )
        {
        typedef BRAINSFit::CommandIterationUpdate<LBFGSBOptimizerType, BSplineTransformType, RegisterImageType>
        CommandIterationUpdateType;
        typename CommandIterationUpdateType::Pointer observer =
        CommandIterationUpdateType::New();
        observer->SetDisplayDeformedImage(m_DisplayDeformedImage);
        observer->SetPromptUserAfterDisplay(m_PromptUserAfterDisplay);
        observer->SetPrintParameters(true);
        observer->SetMovingImage(m_MovingVolume);
        observer->SetFixedImage(m_FixedVolume);
        observer->SetTransform(bsplineTx);
        LBFGSBoptimizer->AddObserver(itk::IterationEvent(), observer);
        }

      for (unsigned int n=0; n<preprocessedFixedImagesList.size(); n++)
         {
         bsplineRegistration->SetFixedImage( n, preprocessedFixedImagesList[n] );
         }

      if( !this->m_InitializeRegistrationByCurrentGenericTransform )
        {
        if( this->m_CurrentGenericTransform.IsNotNull() )
          {
          std::cout << "\nMoving image is warped by initial transform, "
                    << "before it is passed to the BSpline registration.\n" << std::endl;
          typedef float                                                                     VectorComponentType;
          typedef itk::Vector<VectorComponentType, 3> VectorPixelType;
          typedef itk::Image<VectorPixelType,  3>     DisplacementFieldType;
          typename MovingImageType::Pointer warpedMoving1 =
                                              GenericTransformImage<
                                                     MovingImageType,
                                                     MovingImageType,
                                                     DisplacementFieldType>( this->m_MovingVolume,
                                                                             this->m_FixedVolume,
                                                                             this->m_CurrentGenericTransform.GetPointer(),
                                                                             m_BackgroundFillValue,
                                                                             "Linear",
                                                                             false );
          preprocessedMovingImagesList.clear();
          preprocessedMovingImagesList.push_back(warpedMoving1);
          if( m_MovingVolume2.IsNotNull() )
            {
            typename MovingImageType::Pointer warpedMoving2 =
                                                GenericTransformImage<
                                                MovingImageType,
                                                MovingImageType,
                                                DisplacementFieldType>( this->m_MovingVolume2,
                                                                        this->m_FixedVolume2,
                                                                        this->m_CurrentGenericTransform.GetPointer(),
                                                                        m_BackgroundFillValue,
                                                                        "Linear",
                                                                        false );
            preprocessedMovingImagesList.push_back(warpedMoving2);
            }
          for (unsigned int n=0; n<preprocessedMovingImagesList.size(); n++)
            {
            bsplineRegistration->SetMovingImage( n, preprocessedMovingImagesList[n] );
            }

          // Then warp the moving image mask as well if it is not null.
          // Note: m_CostMetricObject should be multi metric object, and masks are get from the first metric object.
          //
          typename ImageMetricType::Pointer firstMetricComponent =
                                                  dynamic_cast<ImageMetricType *>( multiMetric->GetMetricQueue()[0].GetPointer() );
          if( firstMetricComponent.IsNull() )
            {
            std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
            }

          // Moving mask is the same in both metrics.
          typename SpatialObjectType::ConstPointer movingMask = dynamic_cast<SpatialObjectType const *>(
                                                                              firstMetricComponent->GetMovingImageMask() );
          if( movingMask.IsNotNull() )
            {
            typename MaskImageType::Pointer
            warpedMovingMaskImage =
                  GenericTransformImage<
                           MaskImageType,
                           MaskImageType,
                           DisplacementFieldType>(
                                                   ExtractConstPointerToImageMaskFromImageSpatialObject(movingMask).GetPointer(),
                                                   this->m_FixedVolume,
                                                   this->m_CurrentGenericTransform.GetPointer(),
                                                   m_BackgroundFillValue,
                                                   "Linear",
                                                   false );


            const SpatialObjectType *warpedMask =
              ConvertMaskImageToSpatialMask( warpedMovingMaskImage.GetPointer() ).GetPointer();
            firstMetricComponent->SetMovingImageMask( warpedMask );
            // The same for the second metric component in the case of multi modality
            if( preprocessedFixedImagesList.size() > 1)
              {
              typename ImageMetricType::Pointer secondMetricComponent =
                                                    dynamic_cast<ImageMetricType *>( multiMetric->GetMetricQueue()[0].GetPointer() );
              if( secondMetricComponent.IsNull() )
                {
                std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
                }
              secondMetricComponent->SetMovingImageMask( warpedMask );
              }
            }
          }
        else
          {
          for (unsigned int n=0; n<preprocessedMovingImagesList.size(); n++)
            {
            bsplineRegistration->SetMovingImage( n, preprocessedMovingImagesList[n] );
            }
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
        for (unsigned int n=0; n<preprocessedMovingImagesList.size(); n++)
          {
          bsplineRegistration->SetMovingImage( n, preprocessedMovingImagesList[n] );
          }

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

        std::cout << "Stop condition from LBFGSBoptimizer."
                  << bsplineRegistration->GetOptimizer()->GetStopConditionDescription() << std::endl;
        }
      catch( itk::ExceptionObject & e )
        {
        itkGenericExceptionMacro( << "Exception caught: " << e << std::endl );
        }

      // Add the optimized bspline tranform to the current generic transform
      //
      this->m_CurrentGenericTransform->AddTransform( bsplineTx );

      if( this->m_DebugLevel > 4 )
        {
        std::cout << "Final Parameters = " << std::endl;
        std::cout << bsplineTx->GetParameters() << std::endl;
        }
      }
    else if( currentTransformType == "SyN" )
      {
#ifdef USE_ANTS
      //
      // SyN registration metric
      //
      std::string whichmetric = "cc"; // default value
      if( this->m_SyNMetricType == "MMI" )
        {
        whichmetric = "mattes";
        }
      else if( this->m_SyNMetricType == "MSE" )
        {
        whichmetric = "meansquares";
        }
      else if( this->m_SyNMetricType == "NC" )
        {
        whichmetric = "cc";
        }
      else if( this->m_SyNMetricType == "MIH" )
        {
        whichmetric = "mi";
        }
      // Either current m_CurrentGenericTransform or m_RestoreState
      // are used to initialize the SyN registration.
      typename CompositeTransformType::Pointer internalSyNSavedState = CompositeTransformType::New();

      typename CompositeTransformType::Pointer outputSyNTransform =
        simpleSynReg<FixedImageType, MovingImageType>( m_FixedVolume,
                                                       m_MovingVolume,
                                                       m_CurrentGenericTransform,
                                                       internalSyNSavedState,
                                                       m_FixedVolume2,
                                                       m_MovingVolume2,
                                                       m_SamplingPercentage,
                                                       whichmetric,
                                                       m_SyNFull,
                                                       m_RestoreState);

      if( outputSyNTransform.IsNull() )
        {
        std::cout << "\n*******Error: the SyN registration has failed.********\n" << std::endl;
        itkGenericExceptionMacro( << "******* Error: the SyN registration has failed." << std::endl );
        }
      else
        {
        if( this->m_SaveState != "" )
          {
          // Write the state to the disk
          if( internalSyNSavedState.IsNotNull() )
            {
            unsigned int numTransforms = internalSyNSavedState->GetNumberOfTransforms();
            // If the last two transforms are displacement field transforms,
            // we add their inverse displacement field to the saved state composite.

            typename CompositeTransformType::TransformType::Pointer oneToEndTransformGeneric =
                internalSyNSavedState->GetNthTransform( numTransforms-2 ).GetPointer();
            std::string oneToEndTransformFileType;
            if ( oneToEndTransformGeneric.IsNotNull() )
              {
              oneToEndTransformFileType = oneToEndTransformGeneric->GetNameOfClass();
              }

            typename CompositeTransformType::TransformType::Pointer endTransformGeneric =
                internalSyNSavedState->GetNthTransform( numTransforms-1 ).GetPointer();
            std::string endTransformFileType;
            if ( endTransformGeneric.IsNotNull() )
              {
              endTransformFileType = endTransformGeneric->GetNameOfClass();
              }

            if( oneToEndTransformFileType == "DisplacementFieldTransform"
               && endTransformFileType == "DisplacementFieldTransform")
              {
              typedef itk::DisplacementFieldTransform<double, 3>                  DisplacementFieldTransformType;
              DisplacementFieldTransformType::Pointer oneToEndTransform =
                static_cast<DisplacementFieldTransformType *>( oneToEndTransformGeneric.GetPointer() );
              DisplacementFieldTransformType::Pointer endTransform =
                static_cast<DisplacementFieldTransformType *>( endTransformGeneric.GetPointer() );
              if ( oneToEndTransform->GetInverseDisplacementField()
                  && endTransform->GetInverseDisplacementField() )
                {
                internalSyNSavedState->RemoveTransform();
                internalSyNSavedState->AddTransform( oneToEndTransform->GetInverseTransform() );
                internalSyNSavedState->AddTransform( endTransform );
                internalSyNSavedState->AddTransform( endTransform->GetInverseTransform() );
                }
              }
            std::cout << "Writing the registration state: " << this->m_SaveState << std::endl;
            typedef itk::TransformFileWriterTemplate<double>                TransformWriterType;
            typename TransformWriterType::Pointer transformWriter =  TransformWriterType::New();
            transformWriter->SetFileName( this->m_SaveState );
            transformWriter->AddTransform( internalSyNSavedState.GetPointer() );
            try
              {
              transformWriter->Update();
              }
            catch( itk::ExceptionObject & excp )
              {
              itkGenericExceptionMacro( << "Exception caught: Cannot write state file: " << excp << std::endl );
              }
            }
          else
            {
            itkGenericExceptionMacro( << "******* Error: Could save the registration state to the disk." << std::endl );
            }
          }

        if( m_CurrentGenericTransform.IsNull() )
          {
          m_CurrentGenericTransform = CompositeTransformType::New();
          }
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
  if( this->m_FixedVolume2.IsNotNull() )
    {
    os << indent << "FixedVolume2:\n" << this->m_FixedVolume2 << std::endl;
    }
  else
    {
    os << indent << "FixedVolume2: IS NULL" << std::endl;
    }
  os << indent << "MovingVolume:\n" <<   this->m_MovingVolume << std::endl;
  if( this->m_MovingVolume2.IsNotNull() )
    {
    os << indent << "MovingVolume2:\n" << this->m_MovingVolume2 << std::endl;
    }
  else
    {
    os << indent << "MovingVolume2: IS NULL" << std::endl;
    }
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
  os << indent << "MinimumStepLength:     [";
  for( unsigned int q = 0; q < this->m_MinimumStepLength.size(); ++q )
    {
    os << this->m_MinimumStepLength[q] << " ";
    }
  os << "]" << std::endl;
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
