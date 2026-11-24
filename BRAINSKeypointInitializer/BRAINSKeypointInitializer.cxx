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
#include "BRAINSKeypointMatcher.h"
#include "BRAINSSIFT3DFeatureDetector.h"
#include "BRAINSThreadControl.h"
#include "GenericTransformImage.h"
#include "itkAffineTransform.h"
#include "itkIO.h"
#include "itkSimilarity3DTransform.h"
#include "itkVersorRigid3DTransform.h"

#include "BRAINSKeypointInitializerCLP.h"

#include <chrono>
#include <fstream>
#include <iostream>

namespace
{
using ImageType = brains::SIFT3DFeatureDetector::ImageType;

/** Transform with its center at c: x' = A (x - c) + c + translation. */
itk::Transform<double, 3, 3>::Pointer
MakeTransform(brains::TransformModel model, const brains::RobustFitResult & fit, const itk::Point<double, 3> & c)
{
  itk::Vector<double, 3>       translation = fit.translation;
  const itk::Vector<double, 3> cv = c.GetVectorFromOrigin();
  translation += fit.matrix * cv - cv;
  if (model == brains::TransformModel::Affine)
  {
    auto tx = itk::AffineTransform<double, 3>::New();
    tx->SetCenter(c);
    tx->SetMatrix(fit.matrix);
    tx->SetTranslation(translation);
    return tx.GetPointer();
  }
  if (model == brains::TransformModel::Similarity)
  {
    auto tx = itk::Similarity3DTransform<double>::New();
    tx->SetCenter(c);
    tx->SetMatrix(fit.matrix, 1e-6);
    tx->SetTranslation(translation);
    return tx.GetPointer();
  }
  auto tx = itk::VersorRigid3DTransform<double>::New();
  tx->SetCenter(c);
  tx->SetMatrix(fit.matrix, 1e-6);
  tx->SetTranslation(translation);
  return tx.GetPointer();
}

double
Seconds(std::chrono::steady_clock::time_point t0)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
} // namespace

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads threadsHolder(numberOfThreads);

  if (fixedVolume.empty() || movingVolume.empty() || outputTransform.empty())
  {
    std::cerr << "--fixedVolume, --movingVolume, and --outputTransform are required" << std::endl;
    return EXIT_FAILURE;
  }
  brains::RobustFitOptions options;
  if (!brains::ParseTransformModel(transformType, options.model) ||
      !brains::ParseRobustEstimator(robustEstimator, options.estimator))
  {
    std::cerr << "unknown --transformType or --robustEstimator" << std::endl;
    return EXIT_FAILURE;
  }
  options.inlierTolerance = inlierTolerance;
  options.iterations = static_cast<unsigned int>(std::max(numberOfIterations, 1));
  options.seed = static_cast<unsigned int>(randomSeed);
  options.localOptimization = localOptimization;

  try
  {
    auto detector = brains::SIFT3DFeatureDetector::New();
    detector->SetSamplingSpacing(samplingSpacing);
    detector->SetPeakThreshold(peakThreshold);
    detector->SetCornerThreshold(cornerThreshold);
    detector->SetMaximumNumberOfKeypoints(static_cast<unsigned int>(std::max(maximumNumberOfKeypoints, 0)));

    auto                     t0 = std::chrono::steady_clock::now();
    const ImageType::Pointer fixed = itkUtil::ReadImage<ImageType>(fixedVolume);
    const ImageType::Pointer moving = itkUtil::ReadImage<ImageType>(movingVolume);
    const double             readSeconds = Seconds(t0);
    t0 = std::chrono::steady_clock::now();
    const brains::SIFT3DFeatures fixedFeatures = detector->Compute(fixed);
    detector->SetInvertIntensity(invertMovingIntensity);
    const brains::SIFT3DFeatures movingFeatures = detector->Compute(moving);
    const double                 detectSeconds = Seconds(t0);
    t0 = std::chrono::steady_clock::now();
    const brains::KeypointMatches matches = brains::MatchFeatures(fixedFeatures, movingFeatures, maximumRatio);
    const brains::RobustFitResult fit = brains::FitTransform(matches, options);
    const double                  fitSeconds = Seconds(t0);

    std::cout << "keypoints fixed " << fixedFeatures.size() << " moving " << movingFeatures.size() << "; matches "
              << matches.size() << "; inliers " << fit.numberOfInliers << " (mean residual " << fit.meanInlierResidual
              << " mm)\n"
              << "seconds: read " << readSeconds << " detect " << detectSeconds << " match+fit " << fitSeconds
              << std::endl;

    if (!outputMatches.empty())
    {
      std::ofstream csv(outputMatches);
      csv << "fixed_x,fixed_y,fixed_z,moving_x,moving_y,moving_z,ratio,inlier\n";
      for (size_t i = 0; i < matches.size(); ++i)
      {
        csv << matches.fixed[i][0] << ',' << matches.fixed[i][1] << ',' << matches.fixed[i][2] << ','
            << matches.moving[i][0] << ',' << matches.moving[i][1] << ',' << matches.moving[i][2] << ','
            << matches.ratio[i] << ',' << (fit.inliers[i] ? 1 : 0) << '\n';
      }
    }

    if (fit.numberOfInliers < static_cast<unsigned int>(std::max(minimumInliers, 0)))
    {
      std::cerr << "only " << fit.numberOfInliers << " inliers; need " << minimumInliers << std::endl;
      return EXIT_FAILURE;
    }

    itk::Point<double, 3> center;
    center.Fill(0.0);
    for (size_t i = 0; i < matches.size(); ++i)
    {
      if (fit.inliers[i])
      {
        for (unsigned int d = 0; d < 3; ++d)
        {
          center[d] += matches.fixed[i][d] / fit.numberOfInliers;
        }
      }
    }
    itk::WriteTransformToDisk<double>(MakeTransform(options.model, fit, center), outputTransform);
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << "BRAINSKeypointInitializer failed: " << err.GetDescription() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
