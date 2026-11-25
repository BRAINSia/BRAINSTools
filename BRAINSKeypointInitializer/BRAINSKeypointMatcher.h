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
#ifndef BRAINSKeypointMatcher_h
#define BRAINSKeypointMatcher_h

#include "BRAINSSIFT3DFeatureDetector.h"
#include "itkMatrix.h"
#include "itkVector.h"

#include <string>
#include <vector>

namespace brains
{

/** Putative correspondences: fixed[i] is matched to moving[i]. */
struct KeypointMatches
{
  std::vector<SIFT3DFeatures::PointType> fixed;
  std::vector<SIFT3DFeatures::PointType> moving;
  std::vector<double>                    ratio; // nearest / second-nearest descriptor distance

  size_t
  size() const
  {
    return fixed.size();
  }
};

/** Nearest-neighbor descriptor matching with the ratio test. */
KeypointMatches
MatchFeatures(const SIFT3DFeatures & fixed, const SIFT3DFeatures & moving, double maximumRatio);

enum class TransformModel
{
  Rigid,
  Similarity,
  Affine
};

/** How hypotheses are sampled and scored.
 *  RANSAC: uniform samples, scored by inlier count.
 *  MSAC:   uniform samples, scored by truncated squared residual.
 *  PROSAC: samples drawn first from the best-ratio matches, scored as MSAC. */
enum class RobustEstimator
{
  RANSAC,
  MSAC,
  PROSAC
};

bool
ParseTransformModel(const std::string & name, TransformModel & model);
bool
ParseRobustEstimator(const std::string & name, RobustEstimator & estimator);

struct RobustFitOptions
{
  TransformModel  model{ TransformModel::Affine };
  RobustEstimator estimator{ RobustEstimator::RANSAC };
  double          inlierTolerance{ 8.0 }; // mm
  unsigned int    iterations{ 4000 };
  unsigned int    seed{ 0 };
  double          minimumScale{ 0.6 };
  double          maximumScale{ 1.6 };
  /** Refit each new best hypothesis on its inliers (LO-RANSAC). */
  bool localOptimization{ false };
};

/** x_moving = matrix * x_fixed + translation. */
struct RobustFitResult
{
  itk::Matrix<double, 3, 3> matrix;
  itk::Vector<double, 3>    translation;
  std::vector<bool>         inliers;
  unsigned int              numberOfInliers{ 0 };
  double                    meanInlierResidual{ 0.0 };
};

/** Similarity hypotheses from 3-point samples, then a refit of the requested
 *  model on the inliers. Throws itk::ExceptionObject with fewer than 3 matches. */
RobustFitResult
FitTransform(const KeypointMatches & matches, const RobustFitOptions & options);

} // namespace brains

#endif
