#ifndef __BRAINSFitSyN_h
#define __BRAINSFitSyN_h

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"

typedef  ants::RegistrationHelper<3>                    RegistrationHelperType;
typedef  RegistrationHelperType::ImageType              ImageType;
typedef  RegistrationHelperType::CompositeTransformType CompositeTransformType;

template <class FixedImageType, class MovingimageType>
typename CompositeTransformType::TransformTypePointer
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

  std::vector<unsigned int> iterations(3);
  iterations[0] = 100; iterations[1] = 70; iterations[2] = 20;
  double       convergenceThreshold = 1e-6;
  unsigned int convergenceWindowSize = 10;
  ants::antscout << "  number of levels = 3 " << std::endl;
  iterationList.push_back(iterations);
  convergenceThresholdList.push_back(convergenceThreshold);
  convergenceWindowSizeList.push_back(convergenceWindowSize);

  std::vector<unsigned int> factors(3);
  factors[0] = 3; factors[1] = 2; factors[2] = 1;
  shrinkFactorsList.push_back(factors);

  std::vector<float> sigmas(3);
  sigmas[0] = 2; sigmas[1] = 1; sigmas[2] = 0;
  smoothingSigmasList.push_back(sigmas);

  float samplingPercentage = 1.0;
  typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::none;
  unsigned int binOption = 200;
  regHelper->AddMetric(curMetric, fixedImage, movingImage, 1.0, samplingStrategy, binOption, 1, samplingPercentage);

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
  // write out transform actually computed, so skip the initial transform
  typename CompositeTransformType::TransformTypePointer resultTransform =
    outputCompositeTransform->GetNthTransform( 1 );

  return resultTransform;
}

#endif
