#ifndef __BRAINSFitBSpline_h
#define __BRAINSFitBSpline_h

#include <itkBSplineDeformableTransform.h>
#include <itkBSplineDeformableTransformInitializer.h>
#include <itkLBFGSBOptimizer.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkImageRegistrationMethod.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>

#if (ITK_VERSION_MAJOR < 4)
#include "itkOptImageToImageMetric.h"
#else
#include "itkImageToImageMetric.h"
#endif
#include "genericRegistrationHelper.h"

/**
  * This class is the BSpline component of the BRAINSFit program developed at
  *the University of Iowa.
  * The master version of this file is always located from the nitric site.
  * http://www.nitrc.org/projects/multimodereg/
  */

template <class RegisterImageType, class ImageMaskSpatialObjectType, class BSplineTransformType>
typename BSplineTransformType::Pointer
DoBSpline(typename BSplineTransformType::Pointer InitializerBsplineTransform,
          typename RegisterImageType::Pointer m_FixedVolume,
          typename RegisterImageType::Pointer m_MovingVolume,
          typename itk::ImageToImageMetric<
            RegisterImageType, RegisterImageType>::Pointer CostMetricObject,
          const double m_MaxBSplineDisplacement,
          const float m_CostFunctionConvergenceFactor,
          const float m_ProjectedGradientTolerance,
          const bool m_DisplayDeformedImage,
          const bool m_PromptUserAfterDisplay)
{
  /*
    *  Begin straightline BSpline optimization, after
    *GTRACT/Common/itkAnatomicalBSplineFilter.
    */

  typedef typename RegisterImageType::Pointer       RegisterImagePointer;
  typedef typename RegisterImageType::ConstPointer  RegisterImageConstPointer;
  typedef typename RegisterImageType::RegionType    RegisterImageRegionType;
  typedef typename RegisterImageType::SizeType      RegisterImageSizeType;
  typedef typename RegisterImageType::SpacingType   RegisterImageSpacingType;
  typedef typename RegisterImageType::PointType     RegisterImagePointType;
  typedef typename RegisterImageType::PixelType     RegisterImagePixelType;
  typedef typename RegisterImageType::DirectionType RegisterImageDirectionType;
  typedef typename RegisterImageType::IndexType     RegisterImageIndexType;

  typedef typename BSplineTransformType::RegionType     TransformRegionType;
  typedef typename TransformRegionType::SizeType        TransformSizeType;
  typedef typename BSplineTransformType::SpacingType    TransformSpacingType;
  typedef typename BSplineTransformType::OriginType     TransformOriginType;
  typedef typename BSplineTransformType::DirectionType  TransformDirectionType;
  typedef typename BSplineTransformType::ParametersType TransformParametersType;


  typedef typename itk::LinearInterpolateImageFunction<
      RegisterImageType,
      double>          InterpolatorType;

  typedef typename itk::ImageRegistrationMethod<
      RegisterImageType,
      RegisterImageType>          RegistrationType;

  typedef typename BSplineTransformType::Pointer           TransformTypePointer;
  typedef typename itk::LBFGSBOptimizer                    LBFGSBOptimizerType;
  typedef typename LBFGSBOptimizerType::Pointer            LBFGSBOptimizerTypePointer;
  typedef typename LBFGSBOptimizerType::ParametersType     OptimizerParameterType;
  typedef typename LBFGSBOptimizerType::ScalesType         OptimizerScalesType;
  typedef typename LBFGSBOptimizerType::BoundSelectionType OptimizerBoundSelectionType;
  typedef typename LBFGSBOptimizerType::BoundValueType     OptimizerBoundValueType;

  typedef typename InterpolatorType::Pointer InterpolatorTypePointer;
  typedef typename RegistrationType::Pointer RegistrationTypePointer;

  typedef typename itk::ResampleImageFilter< RegisterImageType,
      RegisterImageType>     ResampleFilterType;

  // TODO:  Expose these to the command line for consistancy.
  const int m_MaximumNumberOfIterations = 1500;

  const int m_MaximumNumberOfEvaluations = 900;
  const int m_MaximumNumberOfCorrections = 12;


  typename BSplineTransformType::Pointer m_OutputBSplineTransform = BSplineTransformType::New();
  m_OutputBSplineTransform->SetIdentity();
  m_OutputBSplineTransform->SetBulkTransform( InitializerBsplineTransform->GetBulkTransform() );
  m_OutputBSplineTransform->SetFixedParameters( InitializerBsplineTransform->GetFixedParameters() );
  m_OutputBSplineTransform->SetParametersByValue( InitializerBsplineTransform->GetParameters() );

  /** Set up the Registration */
  RegistrationTypePointer registration = RegistrationType::New();
  registration->SetMetric(CostMetricObject);
  LBFGSBOptimizerTypePointer LBFGSBoptimizer = LBFGSBOptimizerType::New();
  registration->SetOptimizer(LBFGSBoptimizer);
  InterpolatorTypePointer interpolator = InterpolatorType::New();
  registration->SetInterpolator(interpolator);
  registration->SetTransform(m_OutputBSplineTransform);

  /** Setup the Registration */
  registration->SetFixedImage(m_FixedVolume);
  registration->SetMovingImage(m_MovingVolume);

  RegisterImageRegionType fixedImageRegion = m_FixedVolume->GetBufferedRegion();

  registration->SetFixedImageRegion(fixedImageRegion);

  registration->SetInitialTransformParameters( m_OutputBSplineTransform->GetParameters() );

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

  OptimizerBoundSelectionType boundSelect( m_OutputBSplineTransform->GetNumberOfParameters() );
  if( vcl_abs(m_MaxBSplineDisplacement) < 1e-12 )
    {
    boundSelect.Fill(0);
    }
  else
    {
    boundSelect.Fill(2);
    }
  OptimizerBoundValueType upperBound( m_OutputBSplineTransform->GetNumberOfParameters() );
  upperBound.Fill(m_MaxBSplineDisplacement);
  OptimizerBoundValueType lowerBound( m_OutputBSplineTransform->GetNumberOfParameters() );
  lowerBound.Fill(-m_MaxBSplineDisplacement);

  LBFGSBoptimizer->SetBoundSelection(boundSelect);
  std::cout << "PRE : " << LBFGSBoptimizer->GetUpperBound().size() << " " << upperBound.size() << std::endl;
  LBFGSBoptimizer->SetUpperBound(upperBound);
  std::cout << "POST: " << LBFGSBoptimizer->GetUpperBound().size() << " " << upperBound.size() << std::endl;

  std::cout << "PRE : " << LBFGSBoptimizer->GetLowerBound().size() << " " << lowerBound.size() << std::endl;
  LBFGSBoptimizer->SetLowerBound(lowerBound);
  std::cout << "POST: " << LBFGSBoptimizer->GetLowerBound().size() << " " << lowerBound.size() << std::endl;

  LBFGSBoptimizer->SetCostFunctionConvergenceFactor(m_CostFunctionConvergenceFactor);
  LBFGSBoptimizer->SetProjectedGradientTolerance(m_ProjectedGradientTolerance);
  LBFGSBoptimizer->SetMaximumNumberOfIterations(m_MaximumNumberOfIterations);
  LBFGSBoptimizer->SetMaximumNumberOfEvaluations(m_MaximumNumberOfEvaluations);
  LBFGSBoptimizer->SetMaximumNumberOfCorrections(m_MaximumNumberOfCorrections);

  // Create the Command observer and register it with the LBFGSBoptimizer.
  // TODO:  make this output optional.

  const bool ObserveIterations = true;
  if( ObserveIterations == true )
    {
    typedef BRAINSFit::CommandIterationUpdate<LBFGSBOptimizerType, BSplineTransformType, RegisterImageType>
      CommandIterationUpdateType;
    typename CommandIterationUpdateType::Pointer observer =
      CommandIterationUpdateType::New();
    observer->SetDisplayDeformedImage(m_DisplayDeformedImage);
    observer->SetPromptUserAfterDisplay(m_PromptUserAfterDisplay);
    observer->SetPrintParameters(false);
    observer->SetMovingImage(m_MovingVolume);
    observer->SetFixedImage(m_FixedVolume);
    observer->SetTransform(m_OutputBSplineTransform);
    LBFGSBoptimizer->AddObserver(itk::IterationEvent(), observer);
    }

  /* Now start the execute function */

  // Add a time probe
  itk::TimeProbesCollectorBase collector;

  std::cout << std::endl << "Starting Registration" << std::endl;

  try
    {
    collector.Start("Registration");
    registration->Update();
    collector.Stop("Registration");
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return NULL;
    }

  collector.Report();
  std::cout << "Stop condition from LBFGSBoptimizer." << LBFGSBoptimizer->GetStopConditionDescription() << std::endl;

  /* This call is required to copy the parameters */
  const LBFGSBOptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
  m_OutputBSplineTransform->SetParametersByValue(finalParameters);
  return m_OutputBSplineTransform;
}

#endif // __BRAINSFitBSpline_H_
