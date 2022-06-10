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
#ifndef __BRAINSFitSyN_h
#define __BRAINSFitSyN_h

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"
#include "GenericTransformImage.h"

namespace SyN
{
using RealType = double;
using RegistrationHelperType = ants::RegistrationHelper<SyN::RealType, 3>;
using ImageType = RegistrationHelperType::ImageType;
using CompositeTransformType = RegistrationHelperType::CompositeTransformType;
} // namespace SyN

template <typename FixedImageType, typename MovingimageType>
typename SyN::CompositeTransformType::Pointer
simpleSynReg(typename FixedImageType::Pointer &                    infixedImage,
             typename MovingimageType::Pointer &                   inmovingImage,
             const typename SyN::CompositeTransformType::Pointer & compositeInitialTransform,
             typename SyN::CompositeTransformType::Pointer &       internalSavedState,
             typename FixedImageType::Pointer &                    infixedImage2 = NULL,
             typename MovingimageType::Pointer &                   inmovingImage2 = NULL,
             SyN::RealType                                         samplingPercentage = 1.0,
             const std::string &                                   whichMetric = "cc",
             const bool                                            synFull = true,
             const typename SyN::CompositeTransformType::Pointer & restoreState = nullptr)
{
  typename SyN::RegistrationHelperType::Pointer regHelper = SyN::RegistrationHelperType::New();
  {
    /*
    constexpr float lowerQuantile = 0.025;
    constexpr float upperQuantile = 0.975;
    const bool  doWinsorize(false);
    regHelper->SetWinsorizeImageIntensities(doWinsorize, lowerQuantile, upperQuantile);
    */
  } {
    const bool doHistogramMatch(true);
    regHelper->SetUseHistogramMatching(doHistogramMatch);
  }

  {
    /*
    const bool doEstimateLearningRateAtEachIteration = true;
    regHelper->SetDoEstimateLearningRateAtEachIteration( doEstimateLearningRateAtEachIteration );
    */
  }

  {
    std::vector<std::vector<unsigned int>> iterationList;
    if (synFull == true)
    {
      std::vector<unsigned int> iterations(4);
      iterations[0] = 1000;
      iterations[1] = 500;
      iterations[2] = 250;
      iterations[3] = 100;
      iterationList.push_back(iterations);
    }
    else
    {
      std::vector<unsigned int> iterations(1);
      iterations[0] = 100;
      iterationList.push_back(iterations);
    }
    regHelper->SetIterations(iterationList);
  }

  {
    std::vector<SyN::RealType> convergenceThresholdList;
    const SyN::RealType        convergenceThreshold = 1e-6;
    convergenceThresholdList.push_back(convergenceThreshold);
    regHelper->SetConvergenceThresholds(convergenceThresholdList);
  }

  {
    std::vector<unsigned int> convergenceWindowSizeList;
    constexpr unsigned int    convergenceWindowSize = 10;
    convergenceWindowSizeList.push_back(convergenceWindowSize);
    regHelper->SetConvergenceWindowSizes(convergenceWindowSizeList);
  }

  {
    std::vector<std::vector<unsigned int>> shrinkFactorsList;
    if (synFull == true)
    {
      // --shrink-factors 3x2x1
      std::vector<unsigned int> factors(4);
      factors[0] = 8;
      factors[1] = 4;
      factors[2] = 2;
      factors[3] = 1;
      shrinkFactorsList.push_back(factors);
    }
    else
    {
      // --shrink-factors 1
      std::vector<unsigned int> factors(1);
      factors[0] = 1;
      shrinkFactorsList.push_back(factors);
    }
    regHelper->SetShrinkFactors(shrinkFactorsList);
  }

  {
    std::vector<std::vector<float>> smoothingSigmasList;
    if (synFull == true)
    {
      // --smoothing-sigmas 3x2x0
      std::vector<float> sigmas(4);
      sigmas[0] = 3;
      sigmas[1] = 2;
      sigmas[2] = 1;
      sigmas[3] = 0;
      smoothingSigmasList.push_back(sigmas);
    }
    else
    {
      // --smoothing-sigmas 0
      std::vector<float> sigmas(1);
      sigmas[0] = 0;
      smoothingSigmasList.push_back(sigmas);
    }
    regHelper->SetSmoothingSigmas(smoothingSigmasList);
  }

  {
    // Force all units to be in physcial space
    std::vector<bool> smoothingSigmasAreInPhysicalUnitsList;
    smoothingSigmasAreInPhysicalUnitsList.push_back(true);
    regHelper->SetSmoothingSigmasAreInPhysicalUnits(smoothingSigmasAreInPhysicalUnitsList);
  }
  //
  // Add metric to regHelper
  //
  {
    // Note that here we run only one stage of registration, so stageID does not change
    constexpr unsigned int stageID = 0;
    // However, this stage can have one or two metrics (for multi-modality registration).
    // All the parameters of the second metric is the same as the first metric except for fixed and moving volumes.
    // Common parameters are:
    // - Metric type (MMI->mattes, MSE->meansquares, NC->cc, MIH->mi)
    typename SyN::RegistrationHelperType::MetricEnumeration curMetric = regHelper->StringToMetricType(whichMetric);
    // - Metric weight
    constexpr SyN::RealType weighting = 1.0;
    // - Sampling strategy (alway random)
    typename SyN::RegistrationHelperType::SamplingStrategy samplingStrategy = SyN::RegistrationHelperType::random;
    // - Sampling percentage (defined by input)
    // - Number of bins
    constexpr int bins = 32;
    // - radius
    constexpr unsigned int radius = 4;
    //
    // Add the first metric with the first mandatory fixed and moving volumes
    //
    using FixedCasterType = itk::CastImageFilter<FixedImageType, SyN::ImageType>;
    typename FixedCasterType::Pointer fixedCaster = FixedCasterType::New();
    fixedCaster->SetInput(infixedImage);
    fixedCaster->Update();
    typename SyN::ImageType::Pointer dblFixedImage = fixedCaster->GetOutput();

    using MovingCasterType = itk::CastImageFilter<MovingimageType, SyN::ImageType>;
    typename MovingCasterType::Pointer movingCaster = MovingCasterType::New();
    movingCaster->SetInput(inmovingImage);
    movingCaster->Update();
    typename SyN::ImageType::Pointer dblMovingImage = movingCaster->GetOutput();

    {
      SyN::RegistrationHelperType::LabeledPointSetType *   fixedLabeledPointSet = nullptr;
      SyN::RegistrationHelperType::LabeledPointSetType *   movingLabeledPointSet = nullptr;
      SyN::RegistrationHelperType::IntensityPointSetType * fixedIntensityPointSet = nullptr;
      SyN::RegistrationHelperType::IntensityPointSetType * movingIntensityPointSet = nullptr;
      bool                                                 useGradientFilter = false;
      bool                                                 useBoundaryPointsOnly = false;
      SyN::RegistrationHelperType::RealType                pointSetSigma = 0.0;
      unsigned int                                         evaluationKNeighborhood = 0;
      SyN::RegistrationHelperType::RealType                alpha = 0.0;
      bool                                                 useAnisotropicCovariances = false;
      SyN::RegistrationHelperType::RealType                intensityDistanceSigma = 0.0;
      SyN::RegistrationHelperType::RealType                euclideanDistanceSigma = 0.0;


      regHelper->AddMetric(
        /*metricType*/ curMetric,
        /*fixedImage*/ dblFixedImage,
        /*movingImage*/ dblMovingImage,
        fixedLabeledPointSet,
        movingLabeledPointSet,
        fixedIntensityPointSet,
        movingIntensityPointSet,
        stageID,
        weighting,
        samplingStrategy,
        bins,
        radius,
        useGradientFilter,
        useBoundaryPointsOnly,
        pointSetSigma,
        evaluationKNeighborhood,
        alpha,
        useAnisotropicCovariances,
        samplingPercentage,
        intensityDistanceSigma,
        euclideanDistanceSigma);

      // Now if fixedVolume2 and movingVolume2 are not NULL, we add the second metric
      // to our only stage of registration. Second metric has the same parameters but different input images.
      if (infixedImage2.IsNotNull() && inmovingImage2.IsNotNull())
      {
        std::cout << "Do Multimodal Registration..." << std::endl;
        typename FixedCasterType::Pointer fixedCaster2 = FixedCasterType::New();
        fixedCaster2->SetInput(infixedImage2);
        fixedCaster2->Update();
        typename SyN::ImageType::Pointer dblFixedImage2 = fixedCaster2->GetOutput();

        typename MovingCasterType::Pointer movingCaster2 = MovingCasterType::New();
        movingCaster2->SetInput(inmovingImage2);
        movingCaster2->Update();
        typename SyN::ImageType::Pointer dblMovingImage2 = movingCaster2->GetOutput();

        regHelper->AddMetric(curMetric,
                             dblFixedImage2,
                             dblMovingImage2,
                             fixedLabeledPointSet,
                             movingLabeledPointSet,
                             fixedIntensityPointSet,
                             movingIntensityPointSet,
                             stageID,
                             weighting,
                             samplingStrategy,
                             bins,
                             radius,
                             useGradientFilter,
                             useBoundaryPointsOnly,
                             pointSetSigma,
                             evaluationKNeighborhood,
                             alpha,
                             useAnisotropicCovariances,
                             samplingPercentage,
                             intensityDistanceSigma,
                             euclideanDistanceSigma);
      }
    }
  }

  {
    // --transform "SyN[0.33,3.0,0.0]"
    constexpr float learningRate = 0.15;
    constexpr float varianceForUpdateField = 3.0;
    constexpr float varianceForTotalField = 0.0;
    regHelper->AddSyNTransform(learningRate, varianceForUpdateField, varianceForTotalField);
  }

  if (restoreState.IsNotNull())
  {
    regHelper->SetRestoreStateTransform(restoreState);
  }
  else
  {
    if (compositeInitialTransform.IsNotNull())
    {
      regHelper->SetMovingInitialTransform(compositeInitialTransform);
    }
  }

  regHelper->SetInitializeTransformsPerStage(true);

  regHelper->SetLogStream(std::cout);

  if (regHelper->DoRegistration() != EXIT_SUCCESS)
  {
    std::cerr << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    std::cerr << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    std::cerr << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
    std::cerr << "FATAL ERROR: REGISTRATION PROCESS WAS UNSUCCESSFUL" << std::endl;
  }
  else
  {
    std::cerr << "Finshed SyN stage" << std::endl;
  }
  // Get the output transform
  typename SyN::CompositeTransformType::Pointer outputCompositeTransform = regHelper->GetModifiableCompositeTransform();
  // Get the registration state file
  internalSavedState = dynamic_cast<SyN::CompositeTransformType *>(regHelper->GetModifiableRegistrationState());
  // return composite result Transform;
  return outputCompositeTransform;
}

#endif
