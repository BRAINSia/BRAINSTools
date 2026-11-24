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
// Recover a known rigid motion between two samplings of one synthetic volume.
// The moving sampling is mirrored, oblique, and anisotropic, so the test also
// covers the detector's canonical grid. With "Model", the motion is recovered
// through a keypoint model of two landmarked references written to disk.
#include "BRAINSKeypointMatcher.h"
#include "BRAINSKeypointModel.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVersor.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>

namespace
{
using ImageType = brains::SIFT3DFeatureDetector::ImageType;
using PointType = itk::Point<double, 3>;
using MatrixType = itk::Matrix<double, 3, 3>;

struct Blob
{
  PointType  center;
  MatrixType precision;
  double     amplitude;
};

std::vector<Blob>
MakeBlobs(unsigned int count)
{
  std::mt19937                           rng(1);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::uniform_real_distribution<double> sd(2.5, 7.0);
  std::vector<Blob>                      blobs;
  while (blobs.size() < count)
  {
    Blob b;
    for (unsigned int d = 0; d < 3; ++d)
    {
      b.center[d] = 45.0 * u(rng);
    }
    if (b.center.GetVectorFromOrigin().GetNorm() > 45.0)
    {
      continue;
    }
    itk::Versor<double>    q;
    itk::Vector<double, 3> axis;
    axis[0] = u(rng);
    axis[1] = u(rng);
    axis[2] = u(rng) + 1e-3;
    q.Set(axis, M_PI * u(rng));
    const MatrixType R = q.GetMatrix();
    MatrixType       S;
    S.Fill(0.0);
    for (unsigned int d = 0; d < 3; ++d)
    {
      const double s = sd(rng);
      S(d, d) = 1.0 / (s * s);
    }
    b.precision = R * S * R.GetTranspose();
    b.amplitude = u(rng) > -0.3 ? 1.0 + u(rng) * 0.5 : -0.8;
    blobs.push_back(b);
  }
  return blobs;
}

double
Evaluate(const std::vector<Blob> & blobs, const PointType & x)
{
  double v = 0.0;
  for (const auto & b : blobs)
  {
    const itk::Vector<double, 3> d = x - b.center;
    if (d.GetSquaredNorm() > 30.0 * 30.0)
    {
      continue;
    }
    v += b.amplitude * std::exp(-0.5 * (d * (b.precision * d)));
  }
  return v;
}

/** Sample f(R^T (y - t)) on the given grid. */
ImageType::Pointer
Sample(const std::vector<Blob> &        blobs,
       const ImageType::SizeType &      size,
       const ImageType::SpacingType &   spacing,
       const ImageType::DirectionType & direction,
       const MatrixType &               R,
       const itk::Vector<double, 3> &   t)
{
  auto image = ImageType::New();
  image->SetRegions(size);
  image->SetSpacing(spacing);
  image->SetDirection(direction);
  itk::Vector<double, 3> half;
  for (unsigned int d = 0; d < 3; ++d)
  {
    half[d] = 0.5 * spacing[d] * (size[d] - 1);
  }
  PointType origin;
  origin.Fill(0.0);
  image->SetOrigin(origin - direction * half);
  image->Allocate();
  const MatrixType Rt = R.GetTranspose();
  for (itk::ImageRegionIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
  {
    PointType y;
    image->TransformIndexToPhysicalPoint(it.GetIndex(), y);
    PointType                    x;
    const itk::Vector<double, 3> v = Rt * (y.GetVectorFromOrigin() - t);
    for (unsigned int d = 0; d < 3; ++d)
    {
      x[d] = v[d];
    }
    it.Set(static_cast<float>(Evaluate(blobs, x)));
  }
  return image;
}
} // namespace

int
main(int argc, char * argv[])
{
  brains::RobustFitOptions options;
  const bool               useModel = argc > 1 && std::string(argv[1]) == "Model";
  if (argc > 1 && !useModel && !brains::ParseRobustEstimator(argv[1], options.estimator))
  {
    std::cerr << "unknown estimator " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }
  const std::vector<Blob> blobs = MakeBlobs(160);

  itk::Versor<double>    truth;
  itk::Vector<double, 3> axis;
  axis[0] = 1.0;
  axis[1] = 0.3;
  axis[2] = 0.2;
  truth.Set(axis, 25.0 * M_PI / 180.0);
  const MatrixType       R = truth.GetMatrix();
  itk::Vector<double, 3> t;
  t[0] = 8.0;
  t[1] = -5.0;
  t[2] = 3.0;

  MatrixType identity;
  identity.SetIdentity();
  itk::Vector<double, 3> zero;
  zero.Fill(0.0);

  ImageType::SizeType    fixedSize = { { 112, 112, 112 } };
  ImageType::SpacingType fixedSpacing;
  fixedSpacing.Fill(1.0);
  ImageType::DirectionType fixedDirection;
  fixedDirection.SetIdentity();
  const ImageType::Pointer fixed = Sample(blobs, fixedSize, fixedSpacing, fixedDirection, identity, zero);

  ImageType::SizeType    movingSize = { { 104, 104, 88 } };
  ImageType::SpacingType movingSpacing;
  movingSpacing[0] = 1.1;
  movingSpacing[1] = 1.1;
  movingSpacing[2] = 1.4;
  itk::Versor<double>    tilt;
  itk::Vector<double, 3> zAxis;
  zAxis[0] = 0.0;
  zAxis[1] = 0.0;
  zAxis[2] = 1.0;
  tilt.Set(zAxis, 10.0 * M_PI / 180.0);
  MatrixType flip;
  flip.SetIdentity();
  flip(0, 0) = -1.0;
  ImageType::DirectionType movingDirection = tilt.GetMatrix() * flip;
  const ImageType::Pointer moving = Sample(blobs, movingSize, movingSpacing, movingDirection, R, t);

  auto                         detector = brains::SIFT3DFeatureDetector::New();
  const brains::SIFT3DFeatures ff = detector->Compute(fixed);
  const brains::SIFT3DFeatures mf = detector->Compute(moving);
  brains::RobustFitResult      fit;
  size_t                       numberOfMatches = 0;
  try
  {
    if (useModel)
    {
      // Second reference: the same volume under another rigid motion, its landmarks moved with it.
      itk::Versor<double>    q2;
      itk::Vector<double, 3> axis2;
      axis2[0] = -0.2;
      axis2[1] = 1.0;
      axis2[2] = 0.4;
      q2.Set(axis2, -15.0 * M_PI / 180.0);
      const MatrixType       R2 = q2.GetMatrix();
      itk::Vector<double, 3> t2;
      t2[0] = -4.0;
      t2[1] = 6.0;
      t2[2] = 2.0;
      brains::KeypointModel model;
      model.references.resize(2);
      model.references[0].features = ff;
      model.references[1].features = detector->Compute(Sample(blobs, fixedSize, fixedSpacing, fixedDirection, R2, t2));
      for (int i = 0; i < 12; ++i)
      {
        const std::string name = "L" + std::to_string(i);
        PointType         p;
        p[0] = 30.0 * std::cos(i * 0.9);
        p[1] = 30.0 * std::sin(i * 1.3);
        p[2] = -30.0 + 5.0 * i;
        model.atlasLandmarks[name] = p;
        model.references[0].landmarks[name] = p;
        model.references[1].landmarks[name] = R2 * p + t2;
      }
      const std::string path = "BRAINSKeypointRecoveryTest_model.bkpm";
      brains::WriteKeypointModel(model, path);
      const brains::KeypointModel loaded = brains::ReadKeypointModel(path);
      std::remove(path.c_str());
      const brains::ModelAlignment al =
        brains::AlignToModel(loaded, mf, 0.9, options, brains::TransformModel::Rigid, 20);
      if (al.referencesUsed != 2 || al.fitLabels.size() != 12)
      {
        std::cerr << "model alignment used " << al.referencesUsed << " references and " << al.fitLabels.size()
                  << " landmarks" << std::endl;
        return EXIT_FAILURE;
      }
      fit.matrix = al.matrix;
      fit.translation = al.translation;
      fit.numberOfInliers = std::min(al.inliersPerReference[0], al.inliersPerReference[1]);
    }
    else
    {
      const brains::KeypointMatches matches = brains::MatchFeatures(ff, mf, 0.9);
      numberOfMatches = matches.size();
      fit = brains::FitTransform(matches, options);
    }
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << "fit failed: " << err.GetDescription() << std::endl;
    return EXIT_FAILURE;
  }

  double sum = 0.0;
  double worst = 0.0;
  int    n = 0;
  for (double x = -40; x <= 40; x += 10)
  {
    for (double y = -40; y <= 40; y += 10)
    {
      for (double z = -40; z <= 40; z += 10)
      {
        itk::Vector<double, 3> p;
        p[0] = x;
        p[1] = y;
        p[2] = z;
        if (p.GetNorm() > 40.0)
        {
          continue;
        }
        const double e = (fit.matrix * p + fit.translation - (R * p + t)).GetNorm();
        sum += e;
        worst = std::max(worst, e);
        ++n;
      }
    }
  }
  const double mean = sum / n;
  std::cout << "keypoints " << ff.size() << " / " << mf.size() << ", matches " << numberOfMatches << ", inliers "
            << fit.numberOfInliers << ", error mean " << mean << " max " << worst << " mm" << std::endl;
  if (fit.numberOfInliers < 20 || mean > 0.5 || worst > 1.5)
  {
    std::cerr << "transform not recovered" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
