#ifndef __BRAINSFitSyN_h
#define __BRAINSFitSyN_h

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"

typedef  ants::RegistrationHelper<3>                    RegistrationHelperType;
typedef  RegistrationHelperType::ImageType              ImageType;
typedef  RegistrationHelperType::CompositeTransformType CompositeTransformType;

template <class FixedImageType, class MovingimageType>
typename CompositeTransformType::Pointer
simpleSynReg( typename FixedImageType::Pointer & fixedImage,
              typename MovingimageType::Pointer & movingImage,
              typename CompositeTransformType::Pointer compositeInitialTransform )
{
  typename RegistrationHelperType::Pointer regHelper = RegistrationHelperType::New();

  std::string whichMetric = "mattes";
  typename RegistrationHelperType::MetricEnumeration curMetric = regHelper->StringToMetricType(whichMetric);
  float lowerQuantile = 0.0;
  float upperQuantile = 1.0;
  bool  doWinsorize(false);
  regHelper->SetWinsorizeImageIntensities(doWinsorize, lowerQuantile, upperQuantile);

  bool doHistogramMatch(true);
  regHelper->SetUseHistogramMatching(doHistogramMatch);

  bool doEstimateLearningRateAtEachIteration = true;
  regHelper->SetDoEstimateLearningRateAtEachIteration( doEstimateLearningRateAtEachIteration );

  std::vector<std::vector<unsigned int> > iterationList;
  std::vector<double>                     convergenceThresholdList;
  std::vector<unsigned int>               convergenceWindowSizeList;
  std::vector<std::vector<unsigned int> > shrinkFactorsList;
  std::vector<std::vector<float> >        smoothingSigmasList;

  // --convergence [150x100x70,1e-6,10]
  std::vector<unsigned int> iterations(3);
  iterations[0] = 150; iterations[1] = 100; iterations[2] = 75;
  double       convergenceThreshold = 1e-6;
  unsigned int convergenceWindowSize = 10;    // ?
  ants::antscout << "  number of levels = 3 " << std::endl;
  iterationList.push_back(iterations);
  convergenceThresholdList.push_back(convergenceThreshold);
  convergenceWindowSizeList.push_back(convergenceWindowSize);

  // --shrink-factors 3x2x1
  std::vector<unsigned int> factors(3);
  factors[0] = 3; factors[1] = 2; factors[2] = 1;
  shrinkFactorsList.push_back(factors);

  // --smoothing-sigmas 3x2x0
  std::vector<float> sigmas(3);
  sigmas[0] = 3; sigmas[1] = 2; sigmas[2] = 0;
  smoothingSigmasList.push_back(sigmas);

  float samplingPercentage = 1.0;
  typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::regular; // Regular
  unsigned int binOption = 32;                                                                          // bins
  regHelper->AddMetric(curMetric, fixedImage, movingImage, 1.0, samplingStrategy, binOption, 1, samplingPercentage);

  // --transform "SyN[0.25,3.0,0.0]"
  float       learningRate = 0.25;
  const float varianceForUpdateField = 3.0;
  const float varianceForTotalField = 0.0;

  regHelper->AddSyNTransform(learningRate, varianceForUpdateField, varianceForTotalField);
  regHelper->SetMovingInitialTransform( compositeInitialTransform );
  regHelper->SetIterations( iterationList );
  regHelper->SetConvergenceWindowSizes( convergenceWindowSizeList );
  regHelper->SetConvergenceThresholds( convergenceThresholdList );
  regHelper->SetSmoothingSigmas( smoothingSigmasList );
  regHelper->SetShrinkFactors( shrinkFactorsList );
  if( regHelper->DoRegistration() == EXIT_FAILURE )
    {
    ants::antscout << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    }
  // Get the output transform
  typename CompositeTransformType::Pointer outputCompositeTransform = regHelper->GetCompositeTransform();

  // return composite result Transform;
  return outputCompositeTransform;
}

#endif
