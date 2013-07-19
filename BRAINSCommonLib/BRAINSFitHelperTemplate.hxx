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

#include "itkLabelObject.h"
#include "itkStatisticsLabelObject.h"
#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkMacro.h"
#include "itkBSplineDeformableTransformInitializer.h"

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

    //     typename MovingImageType::PointType movingCenter =
    // GetCenterOfBrain<MovingImageType>(orientedMovingVolume);
    //     typename FixedImageType::PointType fixedCenter
    // GetCenterOfBrain<FixedImageType>(orientedFixedVolume);
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

      // typename CHMMaskImageType::ConstPointer ClippedMask = movingFindCenter->GetClippedImageMask();
      // itkUtil::WriteImage<CHMMaskImageType>( ClippedMask , std::string("MOVING_MASK.nii.gz"));

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

      // typename CHMMaskImageType::ConstPointer ClippedMask = fixedFindCenter->GetClippedImageMask();

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

    CostMetricObject->SetTransform(currentEulerAngles3D);
    CostMetricObject->Initialize();
    // void QuickSampleParameterSpace(void)
      {
      currentEulerAngles3D->SetRotation(0, 0, 0);
      // Initialize with current guess;
      double max_cc = CostMetricObject->GetValue( currentEulerAngles3D->GetParameters() );
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
            const double current_cc = CostMetricObject->GetValue( currentEulerAngles3D->GetParameters() );
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
        // writer->SetInput(checker->GetOutput());
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
    // amount
    // of mass approximately uniformly distributed.
    typename SpecificInitializerType::Pointer CenteredInitializer =
      SpecificInitializerType::New();

    CenteredInitializer->SetFixedImage(orientedFixedVolume);
    CenteredInitializer->SetMovingImage(orientedMovingVolume);
    CenteredInitializer->SetTransform(initialITKTransform);
    CenteredInitializer->MomentsOn();                    // Use intensity center
                                                         // of
    // mass

    CenteredInitializer->InitializeTransform();
    }
  else                     // can't happen unless an unimplemented CLP option
                           // was
  // added:
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
  m_NumberOfSamples(500000),
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
  m_UseExplicitPDFDerivativesMode("Off"),
  m_MaskInferiorCutOffFromCenter(1000),
  m_SplineGridSize(3, 10),
  m_CostFunctionConvergenceFactor(1e+9),
  m_ProjectedGradientTolerance(1e-5),
  m_MaxBSplineDisplacement(0.0),
  m_ActualNumberOfIterations(0),
  m_PermittedNumberOfIterations(0),
  // m_AccumulatedNumberOfIterationsForAllLevels(0),
  m_DebugLevel(0),
  m_CurrentGenericTransform(NULL),
  m_GenericTransformList(0),
  m_DisplayDeformedImage(false),
  m_PromptUserAfterDisplay(false),
  m_FinalMetricValue(0.0),
  m_ObserveIterations(true),
  m_CostMetricObject(NULL),
  m_UseROIBSpline(0),
  m_PermitParameterVariation(0),
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
  double minimumStepLength,
  typename TransformType::Pointer &
  initialITKTransform)
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

  // TODO:  What do the following two lines really accomplish
  // debug parameter, suppressed from command line
  const bool initialTransformPassThru(false);
  appMutualRegistration->SetInitialTransformPassThruFlag( initialTransformPassThru );
  appMutualRegistration->SetPermitParameterVariation( m_PermitParameterVariation );
  appMutualRegistration->SetNumberOfSamples( m_NumberOfSamples );
  appMutualRegistration->SetRelaxationFactor( m_RelaxationFactor );
  appMutualRegistration->SetMaximumStepLength( m_MaximumStepLength );
  appMutualRegistration->SetMinimumStepLength( minimumStepLength );
  appMutualRegistration->SetTranslationScale( m_TranslationScale );
  appMutualRegistration->SetReproportionScale( m_ReproportionScale );
  appMutualRegistration->SetSkewScale( m_SkewScale );

  // NOTE: binary masks are set for the cost metric object!!!
  appMutualRegistration->SetFixedImage(    m_FixedVolume    );
  appMutualRegistration->SetMovingImage(   m_MovingVolume   );
  appMutualRegistration->SetCostMetricObject( this->m_CostMetricObject );
  appMutualRegistration->SetForceMINumberOfThreads( this->m_ForceMINumberOfThreads );

  appMutualRegistration->SetBackgroundFillValue(   m_BackgroundFillValue   );

  appMutualRegistration->SetInitialTransform( initialITKTransform );
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

  typename TransformType::Pointer finalTransform;
  try
    {
    appMutualRegistration->Update();
    finalTransform = appMutualRegistration->GetTransform();

    // Find the metric value (It is needed when logFileReport flag is ON).
    //this->m_FinalMetricValue = appMutualRegistration->GetFinalMetricValue();

    this->m_ActualNumberOfIterations = appMutualRegistration->GetActualNumberOfIterations();
    this->m_PermittedNumberOfIterations = numberOfIterations;
    // this->m_AccumulatedNumberOfIterationsForAllLevels +=
    // appMutualRegistration->GetActualNumberOfIterations();
    }
  catch( itk::ExceptionObject& err )
    {
    // pass exception to caller
    itkGenericExceptionMacro(<< "Exception caught during registration: " << err);
    }

  // Put the transform on the CurrentTransformList
  // Initialize next level of transformations with previous transform
  // result
  this->m_CurrentGenericTransform = finalTransform;
}

template <class FixedImageType, class MovingImageType>
void
BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::Update(void)
{
  typedef COMMON_MMI_METRIC_TYPE<FixedImageType, MovingImageType> MattesMutualInformationMetricType;
  unsigned currentTransformId = 0;

  if( std::string(this->m_InitializeTransformMode) != "Off" )
    {
    m_GenericTransformList.resize(m_TransformType.size() + 1);
    }
  else
    {
    m_GenericTransformList.resize( m_TransformType.size() );
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
  if( localInitializeTransformMode != "Off" )
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
    m_CurrentGenericTransform = initialITKTransform.GetPointer();
    localInitializeTransformMode = "Off";        // Now reset to Off once
    // initialization is done.

    // Now if necessary clip the images based on m_MaskInferiorCutOffFromCenter
    DoCenteredTransformMaskClipping<TransformType,
                                    FixedImageType::ImageDimension>(
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

    m_GenericTransformList[currentTransformId++] = initialITKTransform;
    }
  for( unsigned int currentTransformIndex = 0;
       currentTransformIndex < m_TransformType.size();
       currentTransformIndex++ )
    {
    // m_AccumulatedNumberOfIterationsForAllLevels +=
    // localNumberOfIterations[currentTransformIndex];
    const std::string currentTransformType(m_TransformType[currentTransformIndex]);
    std::cout << "\n\n\n=============================== "
              << "Starting Transform Estimations for "
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
      typedef VersorRigid3DTransformType           TransformType;
      typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
      // const int NumberOfEstimatedParameter = 6;

      //
      // Process the initialITKTransform as VersorRigid3DTransform:
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        try
          {
          const std::string transformFileType = m_CurrentGenericTransform->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
            VersorRigid3DTransformType::Pointer tempInitializerITKTransform = ComputeRigidTransformFromGeneric(
                m_CurrentGenericTransform.GetPointer() );
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
        {
        // Special optimizations only for the MMI metric
        // that need adjusting based on both the type of metric, and
        // the "dimensionality" of the transform being adjusted.
        typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
          dynamic_cast<MattesMutualInformationMetricType *>(this->m_CostMetricObject.GetPointer() );
        if( test_MMICostMetric.IsNotNull() )
          {
          const bool UseExplicitPDFDerivatives =
            ( m_UseExplicitPDFDerivativesMode == "ON" || m_UseExplicitPDFDerivativesMode == "AUTO" ) ? true : false;
          test_MMICostMetric->SetUseExplicitPDFDerivatives(UseExplicitPDFDerivatives);
          }
        }

      this->FitCommonCode<TransformType, OptimizerType, MetricType>
        (localNumberOfIterations[currentTransformIndex],
        localMinimumStepLength[currentTransformIndex],
        initialITKTransform);
      localInitializeTransformMode = "Off";   // Now turn of the initiallize
                                              // code to off
      }
    else if( currentTransformType == "ScaleVersor3D" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef ScaleVersor3DTransformType    TransformType;
      typedef itk::VersorTransformOptimizer OptimizerType;
      // const int NumberOfEstimatedParameter = 9;

      //
      // Process the initialITKTransform as ScaleVersor3DTransform:
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        try
          {
          const std::string transformFileType = m_CurrentGenericTransform->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
              dynamic_cast<ScaleVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
                m_CurrentGenericTransform.GetPointer() );
            AssignRigid::AssignConvertedTransform( initialITKTransform, tempInitializerITKTransform.GetPointer() );
            }
          else              // || transformFileType ==
                            // "ScaleSkewVersor3DTransform"
          // ||
          // transformFileType == "AffineTransform"
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
        {
        // Special optimizations only for the MMI metric
        // that need adjusting based on both the type of metric, and
        // the "dimensionality" of the transform being adjusted.
        typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
          dynamic_cast<MattesMutualInformationMetricType *>(this->m_CostMetricObject.GetPointer() );
        if( test_MMICostMetric.IsNotNull() )
          {
          const bool UseExplicitPDFDerivatives =
            ( m_UseExplicitPDFDerivativesMode == "ON" || m_UseExplicitPDFDerivativesMode == "AUTO" ) ? true : false;
          test_MMICostMetric->SetUseExplicitPDFDerivatives(UseExplicitPDFDerivatives);
          }
        }
      // #include "FitCommonCode.tmpl"
      this->FitCommonCode<TransformType, OptimizerType, MetricType>
        (localNumberOfIterations[currentTransformIndex],
        localMinimumStepLength[currentTransformIndex],
        initialITKTransform);
      localInitializeTransformMode = "Off";   // Now turn of the initiallize
                                              // code to off
      }
    else if( currentTransformType == "ScaleSkewVersor3D" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef ScaleSkewVersor3DTransformType TransformType;
      typedef itk::VersorTransformOptimizer  OptimizerType;
      // const int NumberOfEstimatedParameter = 15;

      //
      // Process the initialITKTransform as ScaleSkewVersor3D:
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        try
          {
          const std::string transformFileType = m_CurrentGenericTransform->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
              dynamic_cast<ScaleVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
              dynamic_cast<ScaleSkewVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
                m_CurrentGenericTransform.GetPointer() );
            AssignRigid::AssignConvertedTransform( initialITKTransform, tempInitializerITKTransform.GetPointer() );
            }
          else              // || transformFileType == "AffineTransform" ||
          // transformFileType
          // == "ScaleVersor3DTransform"
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
        {
        // Special optimizations only for the MMI metric
        // that need adjusting based on both the type of metric, and
        // the "dimensionality" of the transform being adjusted.
        typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
          dynamic_cast<MattesMutualInformationMetricType *>(this->m_CostMetricObject.GetPointer() );
        if( test_MMICostMetric.IsNotNull() )
          {
          const bool UseExplicitPDFDerivatives =
            ( m_UseExplicitPDFDerivativesMode == "ON" || m_UseExplicitPDFDerivativesMode == "AUTO" ) ? true : false;
          test_MMICostMetric->SetUseExplicitPDFDerivatives(UseExplicitPDFDerivatives);
          }
        }
      // #include "FitCommonCode.tmpl"
      this->FitCommonCode<TransformType, OptimizerType, MetricType>
        (localNumberOfIterations[currentTransformIndex],
        localMinimumStepLength[currentTransformIndex],
        initialITKTransform);
      localInitializeTransformMode = "Off";   // Now turn of the initiallize
                                              // code to off
      }
    else if( currentTransformType == "Affine" )
      {
      //  Choose TransformType for the itk registration class template:
      typedef itk::AffineTransform<double, Dimension>  TransformType;
      typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
      // const int NumberOfEstimatedParameter = 12;

      //
      // Process the initialITKTransform
      //
      TransformType::Pointer initialITKTransform = TransformType::New();
      initialITKTransform->SetIdentity();
      if( m_CurrentGenericTransform.IsNotNull() )
        {
        try
          {
          const std::string transformFileType = m_CurrentGenericTransform->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
              dynamic_cast<ScaleVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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
              dynamic_cast<ScaleSkewVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(initialITKTransform,
                                                  tempInitializerITKTransform);
            }
          else if( transformFileType == "AffineTransform" )
            {
            const AffineTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<AffineTransformType const *>( m_CurrentGenericTransform.GetPointer() );
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

        {
        // Special optimizations only for the MMI metric
        // that need adjusting based on both the type of metric, and
        // the "dimensionality" of the transform being adjusted.
        typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
          dynamic_cast<MattesMutualInformationMetricType *>(this->m_CostMetricObject.GetPointer() );
        if( test_MMICostMetric.IsNotNull() )
          {
          const bool UseExplicitPDFDerivatives =
            ( m_UseExplicitPDFDerivativesMode == "ON" || m_UseExplicitPDFDerivativesMode == "AUTO" ) ? true : false;
          test_MMICostMetric->SetUseExplicitPDFDerivatives(UseExplicitPDFDerivatives);
          }
        }
      // #include "FitCommonCode.tmpl"
      this->FitCommonCode<TransformType, OptimizerType, MetricType>
        (localNumberOfIterations[currentTransformIndex],
        localMinimumStepLength[currentTransformIndex],
        initialITKTransform);
      localInitializeTransformMode = "Off";   // Now turn of the initiallize
                                              // code to off
      }
    else if( currentTransformType == "BSpline" )
      {
      //
      // Process the bulkAffineTransform for BSpline's BULK
      //
      AffineTransformType::Pointer bulkAffineTransform =
        AffineTransformType::New();
      bulkAffineTransform->SetIdentity();

      typedef itk::Image<float, 3> RegisterImageType;

      BSplineTransformType::Pointer outputBSplineTransform =
        BSplineTransformType::New();
      outputBSplineTransform->SetIdentity();

      BSplineTransformType::Pointer initialBSplineTransform =
        BSplineTransformType::New();
      initialBSplineTransform->SetIdentity();

        {
        typedef BSplineTransformType::RegionType
          TransformRegionType;
        typedef TransformRegionType::SizeType
          TransformSizeType;
        typedef itk::BSplineDeformableTransformInitializer<BSplineTransformType,
                                                                 RegisterImageType> InitializerType;
        InitializerType::Pointer transformInitializer = InitializerType::New();
        transformInitializer->SetTransform(initialBSplineTransform);

        if( m_UseROIBSpline )
          {
          ImageMaskSpatialObjectType::Pointer roiMask = ImageMaskSpatialObjectType::New();
          if( m_MovingBinaryVolume.GetPointer() != NULL )
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
            if( m_FixedBinaryVolume.GetPointer() != NULL )
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
          else if( m_FixedBinaryVolume.GetPointer() != NULL )
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

          transformInitializer->SetImage(roiImage);
          }
        else
          {
          transformInitializer->SetImage(m_FixedVolume);
          }

        TransformSizeType tempGridSize;
        tempGridSize[0] = m_SplineGridSize[0];
        tempGridSize[1] = m_SplineGridSize[1];
        tempGridSize[2] = m_SplineGridSize[2];
        transformInitializer->SetGridSizeInsideTheImage(tempGridSize);
        transformInitializer->InitializeTransform();

        std::cout << "BSpline initialized: " << initialBSplineTransform << std::endl;
        }

      if( m_CurrentGenericTransform.IsNotNull() )
        {
        try
          {
          const std::string transformFileType = m_CurrentGenericTransform->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform,
                                                  tempInitializerITKTransform);
            initialBSplineTransform->SetBulkTransform(bulkAffineTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform,
                                                  tempInitializerITKTransform);
            initialBSplineTransform->SetBulkTransform(bulkAffineTransform);
            }
          else if( transformFileType == "ScaleSkewVersor3DTransform" )
            {
            const ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleSkewVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform,
                                                  tempInitializerITKTransform);
            initialBSplineTransform->SetBulkTransform(bulkAffineTransform);
            }
          else if( transformFileType == "AffineTransform" )
            {
            const AffineTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<AffineTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform,
                                                  tempInitializerITKTransform);
            initialBSplineTransform->SetBulkTransform(bulkAffineTransform);
            }
          else if( transformFileType == "BSplineDeformableTransform" )
            {
            const BSplineTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<BSplineTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }

            initialBSplineTransform->SetBulkTransform(
              tempInitializerITKTransform->GetBulkTransform() );
            BSplineTransformType::ParametersType tempFixedInitialParameters =
              tempInitializerITKTransform->GetFixedParameters();
            BSplineTransformType::ParametersType initialFixedParameters =
              initialBSplineTransform->GetFixedParameters();

            bool checkMatch = true;         // Assume true;
            if( initialFixedParameters.GetSize() != tempFixedInitialParameters.GetSize() )
              {
              checkMatch = false;
              std::cerr << "ERROR INITILIZATION FIXED PARAMETERS DO NOT MATCH: " << initialFixedParameters.GetSize()
                        << " != " << tempFixedInitialParameters.GetSize() << std::endl;
              }
            if( checkMatch ) //  This ramus covers the hypothesis that the FixedParameters
                             //  represent the grid locations of the spline nodes.
              {
              for( unsigned int i = 0; i < initialFixedParameters.GetSize(); ++i )
                {
                if( initialFixedParameters.GetElement(i) != tempFixedInitialParameters.GetElement(i) )
                  {
                  checkMatch = false;
                  std::cerr << "ERROR FIXED PARAMETERS DO NOT MATCH: " << initialFixedParameters.GetElement(i)
                            << " != " << tempFixedInitialParameters.GetElement(i) << std::endl;
                  }
                }
              BSplineTransformType::ParametersType tempInitialParameters =
                tempInitializerITKTransform->GetParameters();
              if( initialBSplineTransform->GetNumberOfParameters() ==
                  tempInitialParameters.Size() )
                {
                initialBSplineTransform->SetFixedParameters(
                  tempFixedInitialParameters);
                initialBSplineTransform->SetParametersByValue(tempInitialParameters);
                }
              else
                {
                // Error, initializing from wrong size transform parameters;
                //  Use its bulk transform only?
                itkGenericExceptionMacro(
                  << "Trouble using the m_CurrentGenericTransform"
                  << "for initializing a BSPlineDeformableTransform:"
                  << std::endl
                  << "The initializing BSplineDeformableTransform has a different"
                  << " number of Parameters, than what is required for the requested grid."
                  << std::endl
                  << "BRAINSFit was only able to use the bulk transform that was before it.");
                }
              }
            else
              {
              itkGenericExceptionMacro(
                << "ERROR:  initialization BSpline transform does not have the same "
                << "parameter dimensions as the one currently specified.");
              }
            }
          else if( transformFileType == "CompositeTransform" )
            {
            itkGenericExceptionMacro( << "Composite transform initializer type found:  "
                                      << transformFileType )
            }
          else
            {
            itkGenericExceptionMacro( << "ERROR:  Invalid transform initializer type found:  "
                                      << transformFileType )
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

      // Special optimizations only for the MMI metric
      // that need adjusting based on both the type of metric, and
      // the "dimensionality" of the transform being adjusted.
      typename MattesMutualInformationMetricType::Pointer test_MMICostMetric =
        dynamic_cast<MattesMutualInformationMetricType *>(this->m_CostMetricObject.GetPointer() );
      if( test_MMICostMetric.IsNotNull() )
        {
        // As recommended in documentation in
        // itkMattesMutualInformationImageToImageMetric
        // "UseExplicitPDFDerivatives = False ... This method is well suited
        // for Transforms with a large number of parameters, such as,
        // BSplineDeformableTransforms."
        const bool UseExplicitPDFDerivatives =
          ( m_UseExplicitPDFDerivativesMode != "ON" || m_UseExplicitPDFDerivativesMode == "AUTO" ) ? false : true;
        test_MMICostMetric->SetUseExplicitPDFDerivatives(UseExplicitPDFDerivatives);
        }

      outputBSplineTransform =
        DoBSpline<RegisterImageType, SpatialObjectType,
                  BSplineTransformType>(
          initialBSplineTransform,
          m_FixedVolume, m_MovingVolume,
          this->m_CostMetricObject.GetPointer(),
          this->m_MaxBSplineDisplacement,
          this->m_CostFunctionConvergenceFactor,
          this->m_ProjectedGradientTolerance,
          this->m_DisplayDeformedImage,
          this->m_PromptUserAfterDisplay);
      if( outputBSplineTransform.IsNull() )
        {
        std::cout
          << "Error -- the BSpline fit has failed." << std::endl;
        std::cout
          << "Error -- the BSpline fit has failed." << std::endl;

        m_ActualNumberOfIterations = 1;
        m_PermittedNumberOfIterations = 1;
        }
      else
        {
        // Initialize next level of transformations with previous transform
        // result
        // TransformList.clear();
        // TransformList.push_back(finalTransform);
        m_CurrentGenericTransform = outputBSplineTransform;
        // Now turn of the initiallize code to off
        localInitializeTransformMode = "Off";
          {
          // HACK:  The BSpline optimizer does not return the correct iteration
          // values.
          m_ActualNumberOfIterations = 1;
          m_PermittedNumberOfIterations = 3;
          }
        }
      }
    else if( currentTransformType == "Composite3D" )
      {
      itkGenericExceptionMacro(<< "Composite Transform is not yet Implemented");
      }
    else if( currentTransformType == "SyN" )
      {
#ifdef USE_ANTS
      //
      // Process the bulkAffineTransform for SyN's transform initializer
      //
      AffineTransformType::Pointer bulkAffineTransform = AffineTransformType::New();
      bulkAffineTransform->SetIdentity();

      CompositeTransformType::Pointer initialSyNTransform = CompositeTransformType::New();

      if( m_CurrentGenericTransform.IsNotNull() )
        {
        try
          {
          const std::string transformFileType = m_CurrentGenericTransform->GetNameOfClass();
          if( transformFileType == "VersorRigid3DTransform" )
            {
            const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<VersorRigid3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform, tempInitializerITKTransform);
            initialSyNTransform->AddTransform(bulkAffineTransform);
            }
          else if( transformFileType == "ScaleVersor3DTransform" )
            {
            const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform, tempInitializerITKTransform);
            initialSyNTransform->AddTransform(bulkAffineTransform);
            }
          else if( transformFileType == "ScaleSkewVersor3DTransform" )
            {
            const ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<ScaleSkewVersor3DTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform, tempInitializerITKTransform);
            initialSyNTransform->AddTransform(bulkAffineTransform);
            }
          else if( transformFileType == "AffineTransform" )
            {
            const AffineTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<AffineTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            AssignRigid::AssignConvertedTransform(bulkAffineTransform, tempInitializerITKTransform);
            initialSyNTransform->AddTransform(bulkAffineTransform);
            }
          else if( transformFileType == "BSplineDeformableTransform" )
            {
            itkGenericExceptionMacro( << "ERROR: Improper transform initializer for SyN registration: "
                                      << "BSpline Transform cannot be used as a transform initializer for SyN registration"
                                      << std::endl);
            }
          else if( transformFileType == "CompositeTransform" )
            {
            itkGenericExceptionMacro( << "ERROR:  Can not initialize SyN with CompositeTranform yet."
                                      << transformFileType );
            const CompositeTransformType::ConstPointer tempInitializerITKTransform =
              dynamic_cast<CompositeTransformType const *>( m_CurrentGenericTransform.GetPointer() );
            if( tempInitializerITKTransform.IsNull() )
              {
              std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
              }
            // AssignRigid::AssignConvertedTransform(bulkAffineTransform, tempInitializerITKTransform);
            initialSyNTransform->AddTransform(bulkAffineTransform);
            }
          else
            {
            itkGenericExceptionMacro( << "ERROR:  Invalid transform initializer type found:  "
                                      << transformFileType );
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

      if( initialSyNTransform.IsNull() )
        {
        std::cout << "\n**********" << std::endl;
        std::cout << "ERORR: Undefined intial transform for SyN registration:" << std::endl;
        std::cout << "SyN registration process cannot be done!" << std::endl;
        std::cout << "************" << std::endl;
        itkGenericExceptionMacro( << "******* Error: Undefined intial transform for SyN registration." << std::endl );
        }
      else
        {
        CompositeTransformType::Pointer outputSyNTransform =
          simpleSynReg<FixedImageType, MovingImageType>( m_FixedVolume,
            m_MovingVolume,
            initialSyNTransform );

        if( outputSyNTransform.IsNull() )
          {
          std::cout << "\n*******Error: the SyN registration has failed.********\n" << std::endl;
          itkGenericExceptionMacro( << "******* Error: the SyN registration has failed." << std::endl );
          }
        else
          {
          // CompositeTransformType has derived from itk::Transform, so we can directly assigne that to the
          // m_CurrentGenericTransform that is a GenericTransformType.
          m_CurrentGenericTransform = outputSyNTransform.GetPointer();
          // Now turn of the initiallize code to off
          localInitializeTransformMode = "Off";
          }
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

    if( currentTransformId > m_GenericTransformList.size() - 1 )
      {
      itkGenericExceptionMacro(
        << "Out of bounds access for transform vector!" << std::endl);
      }

    m_GenericTransformList[currentTransformId] = m_CurrentGenericTransform;
    currentTransformId++;
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
  os << indent << "NumberOfSamples:      " << this->m_NumberOfSamples << std::endl;

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
  // os << indent << "AccumulatedNumberOfIterationsForAllLevels: " <<
  // this->m_AccumulatedNumberOfIterationsForAllLevels << std::endl;

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
