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

namespace // put in anon namespace to suppress shadow declaration warnings.
{
typedef  ants::RegistrationHelper<double,3>                SyNRegistrationHelperType;
typedef  SyNRegistrationHelperType::ImageType              ImageType;
}

template <class FixedImageType, class MovingimageType>
typename itk::CompositeTransform<double,3>::Pointer
simpleSynReg( typename FixedImageType::Pointer & infixedImage,
              typename MovingimageType::Pointer & inmovingImage,
              typename itk::CompositeTransform<double,3>::Pointer compositeInitialTransform,
              typename FixedImageType::Pointer & infixedImage2 = NULL,
              typename MovingimageType::Pointer & inmovingImage2 = NULL,
              double samplingPercentage = 1.0,
              std::string whichMetric = "cc",
              const bool synFull = true,
              typename itk::CompositeTransform<double,3>::Pointer restoreState = NULL)
{
  typename SyNRegistrationHelperType::Pointer regHelper = SyNRegistrationHelperType::New();
    {
    /*
    const float lowerQuantile = 0.025;
    const float upperQuantile = 0.975;
    const bool  doWinsorize(false);
    regHelper->SetWinsorizeImageIntensities(doWinsorize, lowerQuantile, upperQuantile);
    */
    }
    {
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
    std::vector<std::vector<unsigned int> > iterationList;
    if( synFull == true )
      {
      std::vector<unsigned int>  iterations(3);
      iterations[0] = 100;
      iterations[1] = 100;
      iterations[2] = 100;
      iterationList.push_back(iterations);
      }
    else
      {
      std::vector<unsigned int>  iterations(1);
      iterations[0] = 100;
      iterationList.push_back(iterations);
      }
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
    const unsigned int        convergenceWindowSize = 10;
    convergenceWindowSizeList.push_back(convergenceWindowSize);
    regHelper->SetConvergenceWindowSizes( convergenceWindowSizeList );
    }

    {
    std::vector<std::vector<unsigned int> > shrinkFactorsList;
    if( synFull == true )
      {
      // --shrink-factors 3x2x1
      std::vector<unsigned int>   factors(3);
      factors[0] = 3;
      factors[1] = 2;
      factors[2] = 1;
      shrinkFactorsList.push_back(factors);
      }
    else
      {
      // --shrink-factors 1
      std::vector<unsigned int>   factors(1);
      factors[0] = 1;
      shrinkFactorsList.push_back(factors);
      }
    regHelper->SetShrinkFactors( shrinkFactorsList );
    }

    {
    std::vector<std::vector<float> > smoothingSigmasList;
    if( synFull == true )
      {
      // --smoothing-sigmas 3x2x0
      std::vector<float>    sigmas(3);
      sigmas[0] = 2;
      sigmas[1] = 1;
      sigmas[2] = 0;
      smoothingSigmasList.push_back(sigmas);
      }
    else
      {
      // --smoothing-sigmas 0
      std::vector<float>    sigmas(1);
      sigmas[0] = 0;
      smoothingSigmasList.push_back(sigmas);
      }
    regHelper->SetSmoothingSigmas( smoothingSigmasList );
    }

    {
    // Force all units to be in physcial space
    std::vector<bool> smoothingSigmasAreInPhysicalUnitsList;
    smoothingSigmasAreInPhysicalUnitsList.push_back(true);
    regHelper->SetSmoothingSigmasAreInPhysicalUnits( smoothingSigmasAreInPhysicalUnitsList );
    }
    //
    // Add metric to regHelper
    //
    {
    // Note that here we run only one stage of registration, so stageID does not change
    const unsigned int stageID = 0;
    // However, this stage can have one or two metrics (for multi-modality registration).
    // All the parameters of the second metric is the same as the first metric except for fixed and moving volumes.
    // Common parameters are:
    // - Metric type (MMI->mattes, MSE->meansquares, NC->cc, MIH->mi)
    typename SyNRegistrationHelperType::MetricEnumeration curMetric = regHelper->StringToMetricType(whichMetric);
    // - Metric weight
    const double weighting = 1.0;
    // - Sampling strategy (alway random)
    typename SyNRegistrationHelperType::SamplingStrategy samplingStrategy = SyNRegistrationHelperType::random;
    // - Sampling percentage (defined by input)
    // - Number of bins
    const int          bins = 32;
    // - radius
    const unsigned int radius = 4;
    //
    // Add the first metric with the first mandatory fixed and moving volumes
    //
    typedef itk::CastImageFilter<FixedImageType,ImageType> FixedCasterType;
    typename FixedCasterType::Pointer fixedCaster = FixedCasterType::New();
    fixedCaster->SetInput( infixedImage );
    fixedCaster->Update();
    typename ImageType::Pointer dblFixedImage = fixedCaster->GetOutput();

    typedef itk::CastImageFilter<MovingimageType,ImageType> MovingCasterType;
    typename MovingCasterType::Pointer movingCaster = MovingCasterType::New();
    movingCaster->SetInput( inmovingImage );
    movingCaster->Update();
    typename ImageType::Pointer dblMovingImage = movingCaster->GetOutput();

    regHelper->AddMetric(curMetric,
                         dblFixedImage,
                         dblMovingImage,
                         stageID,
                         weighting,
                         samplingStrategy,
                         bins,
                         radius,
                         samplingPercentage);

    // Now if fixedVolume2 and movingVolume2 are not NULL, we add the second metric
    // to our only stage of registration. Second metric has the same parameters but different input images.
    if( infixedImage2.IsNotNull() && inmovingImage2.IsNotNull() )
      {
      std::cout<<"Do Multimodal Registration..."<<std::endl;
      typename FixedCasterType::Pointer fixedCaster2 = FixedCasterType::New();
      fixedCaster2->SetInput( infixedImage2 );
      fixedCaster2->Update();
      typename ImageType::Pointer dblFixedImage2 = fixedCaster2->GetOutput();

      typename MovingCasterType::Pointer movingCaster2 = MovingCasterType::New();
      movingCaster2->SetInput( inmovingImage2 );
      movingCaster2->Update();
      typename ImageType::Pointer dblMovingImage2 = movingCaster2->GetOutput();

      regHelper->AddMetric(curMetric,
                           dblFixedImage2,
                           dblMovingImage2,
                           stageID,
                           weighting,
                           samplingStrategy,
                           bins,
                           radius,
                           samplingPercentage);
      }
    }

    {
    // --transform "SyN[0.33,3.0,0.0]"
    const float learningRate = 0.15;
    const float varianceForUpdateField = 3.0;
    const float varianceForTotalField = 0.0;
    regHelper->AddSyNTransform(learningRate, varianceForUpdateField, varianceForTotalField);
    }

  if( restoreState.IsNotNull() )
    {
    regHelper->SetRestoreStateTransform( restoreState );
    }
  else
    {
    if( compositeInitialTransform.IsNotNull() )
      {
      regHelper->SetMovingInitialTransform( compositeInitialTransform );
      }
    }

  regHelper->SetInitializeTransformsPerStage( true );

  regHelper->SetLogStream(std::cout);

  if( regHelper->DoRegistration() != EXIT_SUCCESS )
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
  typename itk::CompositeTransform<double,3>::Pointer outputCompositeTransform =
    regHelper->GetModifiableCompositeTransform();
  // return composite result Transform;
  return outputCompositeTransform;
}

#endif
