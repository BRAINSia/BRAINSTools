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
#include "BRAINSKeypointModel.h"

#include "itkMacro.h"
#include "itk_eigen.h"
#include ITK_EIGEN(Dense)
#include ITK_EIGEN(Geometry)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>

namespace brains
{
namespace
{
constexpr char     kMagic[8] = { 'B', 'K', 'P', 'M', 'O', 'D', 'E', 'L' };
constexpr uint32_t kVersion = 1;
constexpr uint32_t kEndianProbe = 0x01020304;
constexpr unsigned kLen = SIFT3DFeatures::DescriptorLength;

template <typename T>
void
Put(std::ostream & os, const T & v)
{
  os.write(reinterpret_cast<const char *>(&v), sizeof(T));
}

template <typename T>
T
Get(std::istream & is)
{
  T v{};
  is.read(reinterpret_cast<char *>(&v), sizeof(T));
  if (!is)
  {
    itkGenericExceptionMacro("truncated keypoint model");
  }
  return v;
}

void
PutLandmarks(std::ostream & os, const KeypointLandmarks & lm)
{
  Put<uint32_t>(os, static_cast<uint32_t>(lm.size()));
  for (const auto & [name, p] : lm)
  {
    Put<uint32_t>(os, static_cast<uint32_t>(name.size()));
    os.write(name.data(), static_cast<std::streamsize>(name.size()));
    for (unsigned int d = 0; d < 3; ++d)
    {
      Put<double>(os, p[d]);
    }
  }
}

KeypointLandmarks
GetLandmarks(std::istream & is)
{
  KeypointLandmarks lm;
  const auto        n = Get<uint32_t>(is);
  for (uint32_t i = 0; i < n; ++i)
  {
    std::string name(Get<uint32_t>(is), '\0');
    is.read(name.data(), static_cast<std::streamsize>(name.size()));
    SIFT3DFeatures::PointType p;
    for (unsigned int d = 0; d < 3; ++d)
    {
      p[d] = Get<double>(is);
    }
    lm[name] = p;
  }
  return lm;
}

double
Median(std::vector<double> v)
{
  const size_t h = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + h, v.end());
  if (v.size() % 2)
  {
    return v[h];
  }
  return 0.5 * (v[h] + *std::max_element(v.begin(), v.begin() + h));
}
} // namespace

void
WriteKeypointModel(const KeypointModel & model, const std::string & filename)
{
  std::ofstream os(filename, std::ios::binary);
  if (!os)
  {
    itkGenericExceptionMacro("cannot write " << filename);
  }
  os.write(kMagic, sizeof(kMagic));
  Put(os, kVersion);
  Put(os, kEndianProbe);
  Put(os, model.samplingSpacing);
  Put(os, model.peakThreshold);
  Put(os, model.cornerThreshold);
  PutLandmarks(os, model.atlasLandmarks);
  Put<uint32_t>(os, static_cast<uint32_t>(model.references.size()));
  std::vector<uint8_t> q(kLen);
  for (const auto & ref : model.references)
  {
    PutLandmarks(os, ref.landmarks);
    const auto & f = ref.features;
    Put<uint32_t>(os, static_cast<uint32_t>(f.size()));
    for (size_t i = 0; i < f.size(); ++i)
    {
      for (unsigned int d = 0; d < 3; ++d)
      {
        Put<double>(os, f.locations[i][d]);
      }
      Put<double>(os, f.scales[i]);
      Put<double>(os, f.strengths[i]);
      const float * desc = f.descriptors.data() + i * kLen;
      const float   top = *std::max_element(desc, desc + kLen);
      const float   scale = top > 0.0f ? top / 255.0f : 1.0f;
      for (unsigned int k = 0; k < kLen; ++k)
      {
        q[k] = static_cast<uint8_t>(std::lround(desc[k] / scale));
      }
      Put<float>(os, scale);
      os.write(reinterpret_cast<const char *>(q.data()), kLen);
    }
  }
  if (!os)
  {
    itkGenericExceptionMacro("error writing " << filename);
  }
}

KeypointModel
ReadKeypointModel(const std::string & filename)
{
  std::ifstream is(filename, std::ios::binary);
  if (!is)
  {
    itkGenericExceptionMacro("cannot read " << filename);
  }
  char magic[sizeof(kMagic)];
  is.read(magic, sizeof(magic));
  if (!is || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0)
  {
    itkGenericExceptionMacro("" << filename << " is not a keypoint model");
  }
  if (Get<uint32_t>(is) != kVersion || Get<uint32_t>(is) != kEndianProbe)
  {
    itkGenericExceptionMacro("" << filename << ": unsupported keypoint model version or byte order");
  }
  KeypointModel model;
  model.samplingSpacing = Get<double>(is);
  model.peakThreshold = Get<double>(is);
  model.cornerThreshold = Get<double>(is);
  model.atlasLandmarks = GetLandmarks(is);
  model.references.resize(Get<uint32_t>(is));
  std::vector<uint8_t> q(kLen);
  for (auto & ref : model.references)
  {
    ref.landmarks = GetLandmarks(is);
    auto &     f = ref.features;
    const auto n = Get<uint32_t>(is);
    f.locations.resize(n);
    f.scales.resize(n);
    f.strengths.resize(n);
    f.descriptors.resize(static_cast<size_t>(n) * kLen);
    for (uint32_t i = 0; i < n; ++i)
    {
      for (unsigned int d = 0; d < 3; ++d)
      {
        f.locations[i][d] = Get<double>(is);
      }
      f.scales[i] = Get<double>(is);
      f.strengths[i] = Get<double>(is);
      const auto scale = Get<float>(is);
      is.read(reinterpret_cast<char *>(q.data()), kLen);
      float * desc = f.descriptors.data() + static_cast<size_t>(i) * kLen;
      for (unsigned int k = 0; k < kLen; ++k)
      {
        desc[k] = scale * q[k];
      }
    }
  }
  if (!is)
  {
    itkGenericExceptionMacro("truncated keypoint model " << filename);
  }
  return model;
}

ModelAlignment
AlignToModel(const KeypointModel &    model,
             const SIFT3DFeatures &   subject,
             double                   maximumRatio,
             const RobustFitOptions & referenceOptions,
             TransformModel           finalModel,
             unsigned int             minimumInliersPerReference)
{
  ModelAlignment out;
  if (model.references.empty())
  {
    itkGenericExceptionMacro("keypoint model has no references");
  }
  RobustFitOptions refOptions = referenceOptions;
  refOptions.model = TransformModel::Affine;

  std::vector<const KeypointReference *> used;
  std::vector<RobustFitResult>           fits;
  for (const auto & ref : model.references)
  {
    const KeypointMatches matches = MatchFeatures(ref.features, subject, maximumRatio);
    unsigned int          inliers = 0;
    if (matches.size() >= std::max(4u, minimumInliersPerReference))
    {
      try
      {
        RobustFitResult fit = FitTransform(matches, refOptions);
        inliers = fit.numberOfInliers;
        if (inliers >= minimumInliersPerReference)
        {
          used.push_back(&ref);
          fits.push_back(std::move(fit));
        }
      }
      catch (const itk::ExceptionObject &)
      {}
    }
    out.inliersPerReference.push_back(inliers);
  }
  out.referencesUsed = static_cast<unsigned int>(used.size());
  if (used.empty())
  {
    itkGenericExceptionMacro("no reference reached " << minimumInliersPerReference << " inliers");
  }

  for (const auto & [name, p0] : used.front()->landmarks)
  {
    std::array<std::vector<double>, 3> coords;
    bool                               everywhere = true;
    for (size_t r = 0; r < used.size() && everywhere; ++r)
    {
      const auto it = used[r]->landmarks.find(name);
      if (it == used[r]->landmarks.end())
      {
        everywhere = false;
        break;
      }
      const auto carried = fits[r].matrix * it->second + fits[r].translation;
      for (unsigned int d = 0; d < 3; ++d)
      {
        coords[d].push_back(carried[d]);
      }
    }
    if (everywhere)
    {
      SIFT3DFeatures::PointType m;
      for (unsigned int d = 0; d < 3; ++d)
      {
        m[d] = Median(coords[d]);
      }
      out.subjectLandmarks[name] = m;
    }
  }

  for (const auto & [name, p] : model.atlasLandmarks)
  {
    if (out.subjectLandmarks.count(name))
    {
      out.fitLabels.push_back(name);
    }
  }
  const auto n = static_cast<Eigen::Index>(out.fitLabels.size());
  if (n < 4)
  {
    itkGenericExceptionMacro("only " << n << " landmarks shared by the alignment space and the references");
  }
  Eigen::Matrix3Xd P(3, n);
  Eigen::Matrix3Xd Q(3, n);
  for (Eigen::Index i = 0; i < n; ++i)
  {
    const auto & a = model.atlasLandmarks.at(out.fitLabels[i]);
    const auto & b = out.subjectLandmarks.at(out.fitLabels[i]);
    for (unsigned int d = 0; d < 3; ++d)
    {
      P(d, i) = a[d];
      Q(d, i) = b[d];
    }
  }
  Eigen::Matrix3d A;
  Eigen::Vector3d t;
  if (finalModel == TransformModel::Affine)
  {
    Eigen::MatrixXd Ph(n, 4);
    Ph.leftCols<3>() = P.transpose();
    Ph.col(3).setOnes();
    const Eigen::MatrixXd X = Ph.colPivHouseholderQr().solve(Q.transpose());
    A = X.topRows<3>().transpose();
    t = X.row(3).transpose();
  }
  else
  {
    const Eigen::Matrix4d H = Eigen::umeyama(P, Q, finalModel == TransformModel::Similarity);
    A = H.topLeftCorner<3, 3>();
    t = H.topRightCorner<3, 1>();
  }
  for (unsigned int r = 0; r < 3; ++r)
  {
    out.translation[r] = t(r);
    for (unsigned int c = 0; c < 3; ++c)
    {
      out.matrix(r, c) = A(r, c);
    }
  }
  return out;
}

} // namespace brains
