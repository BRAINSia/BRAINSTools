#ifndef __BRAINSFitSyN_h
#define __BRAINSFitSyN_h

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"

typedef  ants::RegistrationHelper<3>                       SyNRegistrationHelperType;
typedef  SyNRegistrationHelperType::ImageType              ImageType;
typedef  SyNRegistrationHelperType::CompositeTransformType CompositeTransformType;

template <class FixedImageType, class MovingimageType>
typename CompositeTransformType::Pointer
simpleSynReg( typename FixedImageType::Pointer & infixedImage,
              typename MovingimageType::Pointer & inmovingImage,
              typename CompositeTransformType::Pointer compositeInitialTransform )
{
  typename SyNRegistrationHelperType::Pointer regHelper = SyNRegistrationHelperType::New();
    {
    const float lowerQuantile = 0.025;
    const float upperQuantile = 0.975;
    const bool  doWinsorize(false);
    regHelper->SetWinsorizeImageIntensities(doWinsorize, lowerQuantile, upperQuantile);
    }
    {
    const bool doHistogramMatch(true);
    regHelper->SetUseHistogramMatching(doHistogramMatch);
    }
    {
    const bool doEstimateLearningRateAtEachIteration = true;
    regHelper->SetDoEstimateLearningRateAtEachIteration( doEstimateLearningRateAtEachIteration );
    }
    {
    std::vector<std::vector<unsigned int> > iterationList;
    std::vector<unsigned int>               iterations(4);
    iterations[0] = 10000;
    iterations[1] = 500;
    iterations[2] = 500;
    iterations[3] = 100;   // NOTE:  but it gives a reasonable result, 70 converges.
    iterationList.push_back(iterations);
    regHelper->SetIterations( iterationList );
    }
    {
    std::vector<double> convergenceThresholdList;
    const double        convergenceThreshold = 5e-7;
    convergenceThresholdList.push_back(convergenceThreshold);
    regHelper->SetConvergenceThresholds( convergenceThresholdList );
    }
    {
    std::vector<unsigned int> convergenceWindowSizeList;
    const unsigned int        convergenceWindowSize = 25;
    convergenceWindowSizeList.push_back(convergenceWindowSize);
    regHelper->SetConvergenceWindowSizes( convergenceWindowSizeList );
    }

    {
    // --shrink-factors 3x2x1
    std::vector<std::vector<unsigned int> > shrinkFactorsList;
    std::vector<unsigned int>               factors(4);
    factors[0] = 5;
    factors[1] = 4;
    factors[2] = 2;
    factors[3] = 1;
    shrinkFactorsList.push_back(factors);
    regHelper->SetShrinkFactors( shrinkFactorsList );
    }
    {
    // --smoothing-sigmas 3x2x0
    std::vector<std::vector<float> > smoothingSigmasList;
    std::vector<float>               sigmas(4);
    sigmas[0] = 5;
    sigmas[1] = 4;
    sigmas[2] = 2;
    sigmas[3] = 0;
    smoothingSigmasList.push_back(sigmas);
    regHelper->SetSmoothingSigmas( smoothingSigmasList );
    }
    {
    // Force all units to be in physcial space
    std::vector<bool> smoothingSigmasAreInPhysicalUnitsList;
    smoothingSigmasAreInPhysicalUnitsList.push_back(true);
    regHelper->SetSmoothingSigmasAreInPhysicalUnits( smoothingSigmasAreInPhysicalUnitsList );
    }

    {
    const std::string whichMetric = "cc";
    typename SyNRegistrationHelperType::MetricEnumeration curMetric = regHelper->StringToMetricType(whichMetric);
    const double weighting = 1.0;
    typename SyNRegistrationHelperType::SamplingStrategy samplingStrategy = SyNRegistrationHelperType::none;
    const int          bins = 32;
    const unsigned int radius = 4;
    const double       samplingPercentage = 1.0;

    const unsigned int stageID = 0;
    typename itk::CastImageFilter<FixedImageType,
                                  ImageType>::Pointer fixedCaster =
      itk::CastImageFilter<FixedImageType, ImageType>::New();
    fixedCaster->SetInput( infixedImage );
    fixedCaster->Update();
    typename ImageType::Pointer dblFixedImage = fixedCaster->GetOutput();
    typename itk::CastImageFilter<MovingimageType,
                                  ImageType>::Pointer movingCaster =
      itk::CastImageFilter<FixedImageType, ImageType>::New();
    movingCaster->SetInput( inmovingImage );
    movingCaster->Update();
    typename ImageType::Pointer dblMovingImage = movingCaster->GetOutput();
    regHelper->AddMetric(curMetric, dblFixedImage, dblMovingImage, stageID, weighting, samplingStrategy, bins, radius,
                         samplingPercentage);
    }
    {
    // --transform "SyN[0.33,3.0,0.0]"
    const float learningRate = 0.15;
    const float varianceForUpdateField = 3.0;
    const float varianceForTotalField = 0.0;
    regHelper->AddSyNTransform(learningRate, varianceForUpdateField, varianceForTotalField);
    }
  regHelper->SetMovingInitialTransform( compositeInitialTransform );
  regHelper->SetLogStream(std::cout);
  if( regHelper->DoRegistration() != EXIT_SUCCESS )
    {
    ants::antscout << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    ants::antscout << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    ants::antscout << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    ants::antscout << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    }
  else
    {
    ants::antscout << "Finshed SyN stage" << std::endl;
    }
  // Get the output transform
  typename CompositeTransformType::Pointer outputCompositeTransform = regHelper->GetModifiableCompositeTransform();
  // return composite result Transform;
  return outputCompositeTransform;
}

#endif
