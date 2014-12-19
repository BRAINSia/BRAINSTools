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
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkSimilarity3DTransform.h"
#include "BRAINSThreadControl.h"

#include "Slicer3LandmarkIO.h"
#include "GenericTransformImage.h"

namespace
{
typedef  itk::Point<double, 3>              ImagePointType;
typedef  std::vector<ImagePointType>        PointList;
typedef itk::VersorRigid3DTransform<double> VersorRigidTransformType;
typedef itk::Similarity3DTransform<double>  SimilarityTransformType;

// Function to convert a point from std::vector to itk::Point
// this also performs the RAS -> LPS conversion necessary
// from slicer -> ITK
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

// Operator to compute the squared distance between two points
class SquaredPointDistance
{
public:
  explicit SquaredPointDistance(const itk::Point<double, 3>& ctr)
    : m_Point(ctr)
  {
  }

  double operator()(const itk::Point<double, 3>& p)
  {
    return (p - m_Point).GetSquaredNorm();
  }

private:
  itk::Point<double, 3> m_Point;
};

// Function to compute the scaling factor between two sets of points.
// This is the symmetric form given by
//    Berthold K. P. Horn (1987),
//    "Closed-form solution of absolute orientation using unit quaternions,"
//    Journal of the Optical Society of America A, 4:629-642
static double
computeSymmetricScale(const std::vector<itk::Point<double, 3> >& fixedPoints,
                      const std::vector<itk::Point<double, 3> >& movingPoints,
                      const itk::Point<double, 3>& fixedcenter,
                      const itk::Point<double, 3>& movingcenter)
{
  std::vector<double> centeredFixedPoints(fixedPoints.size(), 0.0);
  std::vector<double> centeredMovingPoints(movingPoints.size(), 0.0);

  std::transform(fixedPoints.begin(), fixedPoints.end(),
                 centeredFixedPoints.begin(),
                 SquaredPointDistance(fixedcenter) );

  std::transform(movingPoints.begin(), movingPoints.end(),
                 centeredMovingPoints.begin(),
                 SquaredPointDistance(movingcenter) );

  const double fixedmag = std::accumulate(centeredFixedPoints.begin(),
                                          centeredFixedPoints.end(),
                                          0.0);

  const double movingmag = std::accumulate(centeredMovingPoints.begin(),
                                           centeredMovingPoints.end(),
                                           0.0);

  return sqrt(movingmag / fixedmag);
}
}

static SimilarityTransformType::Pointer DoIt_Similarity(PointList fixedPoints, PointList movingPoints)
{
  // Our input into landmark based initialize will be of this form
  // The format for saving to slicer is defined later
  SimilarityTransformType::Pointer similarityTransform = SimilarityTransformType::New();

  similarityTransform->SetIdentity();
  // workaround a bug in older versions of ITK
  similarityTransform->SetScale(1.0);

  typedef itk::LandmarkBasedTransformInitializer<SimilarityTransformType,
                                                 itk::Image<short, 3>, itk::Image<short, 3> > InitializerType;
  InitializerType::Pointer initializer = InitializerType::New();

  // This expects a VersorRigid3D.  The similarity transform works because
  // it derives from that class
  initializer->SetTransform(similarityTransform);

  initializer->SetFixedLandmarks(fixedPoints);
  initializer->SetMovingLandmarks(movingPoints);
  initializer->InitializeTransform();

  // Compute the scaling factor and add that in
  itk::Point<double, 3> fixedCenter(similarityTransform->GetCenter() );
  itk::Point<double, 3> movingCenter(similarityTransform->GetCenter() + similarityTransform->GetTranslation() );

  const double s = computeSymmetricScale(fixedPoints, movingPoints, fixedCenter, movingCenter);
  similarityTransform->SetScale(s);
  return similarityTransform;
}

static VersorRigidTransformType::Pointer DoIt_Rigid(PointList fixedPoints, PointList movingPoints)
{
  // Our input into landmark based initialize will be of this form
  // The format for saving to slicer is defined later
  VersorRigidTransformType::Pointer rigidTransform = VersorRigidTransformType::New();

  rigidTransform->SetIdentity();

  typedef itk::LandmarkBasedTransformInitializer<VersorRigidTransformType,
                                                 itk::Image<short, 3>, itk::Image<short, 3> > InitializerType;
  InitializerType::Pointer initializer = InitializerType::New();

  // This expects a VersorRigid3D.  The similarity transform works because
  // it derives from that class
  initializer->SetTransform(rigidTransform);

  initializer->SetFixedLandmarks(fixedPoints);
  initializer->SetMovingLandmarks(movingPoints);
  initializer->InitializeTransform();
  return rigidTransform;
}

int main(int argc, char* argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if( saveTransform == "" )
    {
    std::cerr << "An output transform must be specified" << std::endl;
    return EXIT_FAILURE;
    }

  if(  ( ( fixedLandmarks.size() > 0 ) && ( fixedLandmarksFile != "" ) )
       || ( ( movingLandmarks.size() > 0 ) && ( movingLandmarksFile != "" ) ) )
    {
    std::cerr << "Program can accept landmark files with points, or directly specify points, but not both."
              << std::endl;
    return EXIT_FAILURE;
    }

  PointList fixedPoints(0);
  PointList movingPoints(0);

  if( fixedLandmarks.size() > 0 ) //
    {
    fixedPoints.resize(fixedLandmarks.size() );
    movingPoints.resize(movingLandmarks.size() );
    // Convert both points lists to ITK points and convert RAS -> LPS

    std::transform(fixedLandmarks.begin(), fixedLandmarks.end(),
                   fixedPoints.begin(),
                   convertStdVectorToITKPoint);

    std::transform(movingLandmarks.begin(), movingLandmarks.end(),
                   movingPoints.begin(),
                   convertStdVectorToITKPoint);
    }
  else
    {
    // Read fcsv files and make lists
    if( fixedLandmarksFile != ""  || movingLandmarksFile != "" )
      {
      if( movingLandmarksFile == "" || fixedLandmarksFile == "" )
        {
        std::cerr << "Must supply both fixed and moving landmark files" << std::endl;
        return EXIT_FAILURE;
        }
      typedef std::map<std::string, ImagePointType> LocalLandmarksMapType;

      // NOTE: ReadSlicer3toITKLmk returns points in LPS system
      LocalLandmarksMapType fixedLandmarkMap = ReadSlicer3toITKLmk( fixedLandmarksFile );
      LocalLandmarksMapType movingLandmarkMap = ReadSlicer3toITKLmk( movingLandmarksFile );
      for( LocalLandmarksMapType::const_iterator fmapit = fixedLandmarkMap.begin();
           fmapit != fixedLandmarkMap.end(); ++fmapit )
        {
        LocalLandmarksMapType::const_iterator mmapit = movingLandmarkMap.find(fmapit->first);
        if( mmapit != movingLandmarkMap.end() )
          {
          fixedPoints.push_back(fmapit->second);
          movingPoints.push_back(fmapit->second);
          }
        }
      }
    }

  if( fixedLandmarks.size() <= 0 || movingLandmarks.size() <= 0 ||
      fixedLandmarks.size() != movingLandmarks.size() )
    {
    std::cerr << "Fixed and moving landmark lists must be of the same size "
              << "and contain at least one point" << std::endl;
    return EXIT_FAILURE;
    }

  if( transformType != "Translation" && fixedLandmarks.size() < 3 )
    {
    std::cerr << "At least 3 fiducual points must be specified for Rigid or Similarity transforms"
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Transform<double, 3, 3> GenericTransformType;
  GenericTransformType::Pointer genericTransform = ITK_NULLPTR;

  if( transformType == "Rigid" )
    {
    VersorRigidTransformType::Pointer rigidTransform = DoIt_Rigid(fixedPoints, movingPoints);
    // do nothing
    genericTransform = rigidTransform.GetPointer();
    }
  else if( transformType == "Translation" )
    {
    VersorRigidTransformType::Pointer rigidTransform = DoIt_Rigid(fixedPoints, movingPoints);
    // Clear out the computed rotation if we only requested translation
    itk::Versor<double> v;
    v.SetIdentity();
    rigidTransform->SetRotation(v);
    genericTransform = rigidTransform.GetPointer();
    }
  else if( transformType == "Similarity" )
    {
    SimilarityTransformType::Pointer similarityTransform = DoIt_Similarity(fixedPoints, movingPoints);
    genericTransform = similarityTransform.GetPointer();
    }
  else if( transformType == "Affine" )
    {
    // itk::Matrix<double, 3> a =
    //   computeAffineTransform(fixedPoints, movingPoints,
    //                          fixedCenter, movingCenter);
    std::cerr << "Unsupported transform type: " << transformType << std::endl;
    genericTransform = ITK_NULLPTR;
    return EXIT_FAILURE;
    }
  else
    {
    std::cerr << "Unsupported transform type: " << transformType << std::endl;
    genericTransform = ITK_NULLPTR;
    return EXIT_FAILURE;
    }

  itk::WriteTransformToDisk<double>(genericTransform.GetPointer(), saveTransform);

  return EXIT_SUCCESS;
}
