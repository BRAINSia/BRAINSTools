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
#ifndef BRAINSKeypointModel_h
#define BRAINSKeypointModel_h

#include "BRAINSKeypointMatcher.h"

#include <map>
#include <string>
#include <vector>

namespace brains
{

using KeypointLandmarks = std::map<std::string, SIFT3DFeatures::PointType>;

/** Keypoints of one landmarked reference scan; all points in its physical (LPS) space. */
struct KeypointReference
{
  SIFT3DFeatures    features;
  KeypointLandmarks landmarks;
};

/** \brief Precomputed keypoints that align any scan to an alignment space.
 *
 * Each reference scan contributes its keypoints and its own landmarks; the
 * alignment space is defined by named landmarks (for example the BCD
 * template). A subject is matched to every reference, each reference's
 * landmarks are carried into the subject by that pair's affine, the
 * per-landmark median over references is taken, and one transform is fit
 * from the alignment-space landmarks to those medians. An averaged template
 * is not used because it shares few keypoints with an individual scan.
 */
struct KeypointModel
{
  double                         samplingSpacing{ 1.0 };
  double                         peakThreshold{ 0.1 };
  double                         cornerThreshold{ 0.4 };
  KeypointLandmarks              atlasLandmarks;
  std::vector<KeypointReference> references;
};

/** Descriptors are stored as 8 bits each with a per-descriptor scale. */
void
WriteKeypointModel(const KeypointModel & model, const std::string & filename);
KeypointModel
ReadKeypointModel(const std::string & filename);

struct ModelAlignment
{
  /** Maps alignment-space points to subject points: x' = matrix x + translation. */
  itk::Matrix<double, 3, 3> matrix;
  itk::Vector<double, 3>    translation;
  /** Per-landmark median over references, in subject space, for labels every used reference has. */
  KeypointLandmarks         subjectLandmarks;
  std::vector<unsigned int> inliersPerReference;
  unsigned int              referencesUsed{ 0 };
  /** Labels shared by the alignment space and the references, used for the fit. */
  std::vector<std::string> fitLabels;
};

/** referenceOptions configures each reference-to-subject fit (its model is forced to Affine);
 *  model selects the final alignment-space fit. */
ModelAlignment
AlignToModel(const KeypointModel &    model,
             const SIFT3DFeatures &   subject,
             double                   maximumRatio,
             const RobustFitOptions & referenceOptions,
             TransformModel           finalModel,
             unsigned int             minimumInliersPerReference);

} // namespace brains

#endif
