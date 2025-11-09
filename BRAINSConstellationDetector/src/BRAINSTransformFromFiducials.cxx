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
/** Some of the components in this file originated from work in the Slicer4
 * development tree.  Adaptations to work as a command line program were
 * made by Hans J. Johnson to facilitate creating a transform directly from a
 * set of landmark files.
 */

#include <itkTransformFileWriter.h>
#include <itkImage.h>
#include <itkAffineTransform.h>

#include <numeric>
#include <functional>
#include <iterator>

#include "BRAINSTransformFromFiducialsCLP.h"

#include "BRAINSThreadControl.h"

#include "Slicer3LandmarkIO.h"
#include "GenericTransformImage.h"
#include "landmarksConstellationCommon.h"

namespace
{
//// Function to convert a point from std::vector to itk::Point
//// this also performs the RAS -> LPS conversion necessary
//// from slicer -> ITK
static itk::Point<double, 3>
convertStdVectorToITKPoint(const std::vector<float> & vec)
{
  itk::Point<double, 3> p;

  // convert RAS to LPS
  p[0] = -vec[0];
  p[1] = -vec[1];
  p[2] = vec[2];
  return p;
}
} // namespace

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if (saveTransform.empty())
  {
    std::cerr << "An output transform must be specified" << std::endl;
    return EXIT_FAILURE;
  }

  if (((!fixedLandmarks.empty()) && (!fixedLandmarksFile.empty())) ||
      ((!movingLandmarks.empty()) && (!movingLandmarksFile.empty())))
  {
    std::cerr << "Program can accept landmark files with points, or directly specify points, but not both."
              << std::endl;
    return EXIT_FAILURE;
  }

  PointList fixedPoints(0);
  PointList movingPoints(0);

  if (!fixedLandmarks.empty()) //
  {
    fixedPoints.resize(fixedLandmarks.size());
    movingPoints.resize(movingLandmarks.size());
    // Convert both points lists to ITK points and convert RAS -> LPS

    std::transform(fixedLandmarks.begin(), fixedLandmarks.end(), fixedPoints.begin(), convertStdVectorToITKPoint);

    std::transform(movingLandmarks.begin(), movingLandmarks.end(), movingPoints.begin(), convertStdVectorToITKPoint);
  }
  else
  {
    // Read fcsv files and make lists
    if (!fixedLandmarksFile.empty() || !movingLandmarksFile.empty())
    {
      if (movingLandmarksFile.empty() || fixedLandmarksFile.empty())
      {
        std::cerr << "Must supply both fixed and moving landmark files" << std::endl;
        return EXIT_FAILURE;
      }
      using LocalLandmarksMapType = std::map<std::string, ImagePointType>;

      // NOTE: ReadSlicer3toITKLmk returns points in LPS system
      const LocalLandmarksMapType fixedLandmarkMap = ReadSlicer3toITKLmk(fixedLandmarksFile);
      LocalLandmarksMapType       movingLandmarkMap = ReadSlicer3toITKLmk(movingLandmarksFile);
      for (auto & fmapit : fixedLandmarkMap)
      {
        auto mmapit = movingLandmarkMap.find(fmapit.first);
        if (mmapit != movingLandmarkMap.end())
        {
          fixedPoints.push_back(fmapit.second);
          movingPoints.push_back(fmapit.second);
        }
      }
    }
  }

  if (fixedLandmarks.empty() || movingLandmarks.empty() || fixedLandmarks.size() != movingLandmarks.size())
  {
    std::cerr << "Fixed and moving landmark lists must be of the same size "
              << "and contain at least one point" << std::endl;
    return EXIT_FAILURE;
  }

  if (transformType != "Translation" && fixedLandmarks.size() < 3)
  {
    std::cerr << "At least 3 fiducual points must be specified for Rigid or Similarity transforms" << std::endl;
    return EXIT_FAILURE;
  }

  using GenericTransformType = itk::Transform<double, 3, 3>;
  GenericTransformType::Pointer genericTransform = nullptr;

  if (transformType == "Rigid")
  {
    const VersorRigidTransformType::Pointer rigidTransform =
      ComputeRigidTransformFromLandmarkLists(fixedPoints, movingPoints);
    // do nothing
    genericTransform = rigidTransform.GetPointer();
  }
  else if (transformType == "Translation")
  {
    const VersorRigidTransformType::Pointer rigidTransform =
      ComputeRigidTransformFromLandmarkLists(fixedPoints, movingPoints);
    // Clear out the computed rotation if we only requested translation
    itk::Versor<double> v;
    v.SetIdentity();
    rigidTransform->SetRotation(v);
    genericTransform = rigidTransform.GetPointer();
  }
  else if (transformType == "Similarity")
  {
    const SimilarityTransformType::Pointer similarityTransform = DoIt_Similarity(fixedPoints, movingPoints);
    genericTransform = similarityTransform.GetPointer();
  }
  /*else if (transformType == "Affine")
  {
     itk::Matrix<double, 3> a =
       computeAffineTransform(fixedPoints, movingPoints,
                              fixedCenter, movingCenter);
  }*/
  else
  {
    std::cerr << "Unsupported transform type: " << transformType << std::endl;
    genericTransform = nullptr;
    return EXIT_FAILURE;
  }

  itk::WriteTransformToDisk<double>(genericTransform.GetPointer(), saveTransform);

  return EXIT_SUCCESS;
}
