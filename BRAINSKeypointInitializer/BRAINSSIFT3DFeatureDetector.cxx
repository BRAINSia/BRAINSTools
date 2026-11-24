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
/*
 * Portions adapted from SIFT3D (https://github.com/bbrister/SIFT3D),
 * Copyright (c) 2015-2016 Blaine Rister et al., MIT License; the full notice
 * is in BRAINSSIFT3DFeatureDetector.h.
 */
#include "BRAINSSIFT3DFeatureDetector.h"

#include "itkIdentityTransform.h"
#include "itkImageRegionConstIterator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMultiThreaderBase.h"
#include "itkResampleImageFilter.h"
#include "itk_eigen.h"
#include ITK_EIGEN(Dense)

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <functional>
#include <limits>
#include <mutex>
#include <tuple>

namespace brains
{
namespace
{
constexpr double kMaxEigRatio = 0.90;
constexpr double kOriGradThresh = 1e-10;
constexpr float  kBaryEps = FLT_EPSILON * 10.0f;
constexpr double kOriSigFctr = 1.5;
constexpr double kOriRadFctr = 3.0;
constexpr double kDescSigFctr = 7.071067812;
constexpr double kDescRadFctr = 2.0;
constexpr int    kHistPerDim = 4;
constexpr int    kIcosVert = 12;
constexpr int    kIcosFaces = 20;
constexpr int    kNumHist = kHistPerDim * kHistPerDim * kHistPerDim;
constexpr float  kTruncThresh = 0.2f * 128.0f / static_cast<float>(SIFT3DFeatures::DescriptorLength);
static_assert(kNumHist * kIcosVert == SIFT3DFeatures::DescriptorLength, "descriptor length");

using Vec3 = Eigen::Vector3f;

/** Dense float volume with x fastest. */
struct Volume
{
  int                nx{ 0 }, ny{ 0 }, nz{ 0 };
  std::vector<float> v;

  void
  Resize(int x, int y, int z)
  {
    nx = x;
    ny = y;
    nz = z;
    v.assign(static_cast<size_t>(x) * y * z, 0.0f);
  }
  float &
  operator()(int x, int y, int z)
  {
    return v[(static_cast<size_t>(z) * ny + y) * nx + x];
  }
  float
  operator()(int x, int y, int z) const
  {
    return v[(static_cast<size_t>(z) * ny + y) * nx + x];
  }
  Vec3
  Gradient(int x, int y, int z, float unit) const
  {
    const float f = 0.5f / unit;
    return { f * ((*this)(x + 1, y, z) - (*this)(x - 1, y, z)),
             f * ((*this)(x, y + 1, z) - (*this)(x, y - 1, z)),
             f * ((*this)(x, y, z + 1) - (*this)(x, y, z - 1)) };
  }
};

void
ParallelOverSlices(int nz, const std::function<void(int)> & fn)
{
  itk::MultiThreaderBase::New()->ParallelizeArray(
    0, static_cast<itk::SizeValueType>(nz), [&fn](itk::SizeValueType z) { fn(static_cast<int>(z)); }, nullptr);
}

/** Separable Gaussian blur in voxel units, clamping at the edges. */
void
GaussianBlur(const Volume & src, Volume & dst, double sigma)
{
  const int          half = std::max(static_cast<int>(std::ceil(sigma * 3.0)), 1);
  std::vector<float> k(2 * half + 1);
  double             acc = 0.0;
  for (int i = -half; i <= half; ++i)
  {
    const double x = i / (sigma + DBL_EPSILON);
    k[i + half] = static_cast<float>(std::exp(-0.5 * x * x));
    acc += k[i + half];
  }
  for (auto & w : k)
  {
    w = static_cast<float>(w / acc);
  }
  Volume a = src;
  Volume b;
  b.Resize(src.nx, src.ny, src.nz);
  const int dims[3] = { src.nx, src.ny, src.nz };
  for (int d = 0; d < 3; ++d)
  {
    ParallelOverSlices(src.nz, [&](int z) {
      for (int y = 0; y < src.ny; ++y)
      {
        for (int x = 0; x < src.nx; ++x)
        {
          const int c[3] = { x, y, z };
          float     s = 0.0f;
          for (int i = -half; i <= half; ++i)
          {
            int p[3] = { x, y, z };
            p[d] = std::clamp(c[d] + i, 0, dims[d] - 1);
            s += k[i + half] * a(p[0], p[1], p[2]);
          }
          b(x, y, z) = s;
        }
      }
    });
    std::swap(a, b);
  }
  dst = std::move(a);
}

void
Downsample2x(const Volume & src, Volume & dst)
{
  dst.Resize(src.nx / 2, src.ny / 2, src.nz / 2);
  for (int z = 0; z < dst.nz; ++z)
  {
    for (int y = 0; y < dst.ny; ++y)
    {
      for (int x = 0; x < dst.nx; ++x)
      {
        dst(x, y, z) = src(2 * x, 2 * y, 2 * z);
      }
    }
  }
}

/** Faces of a regular icosahedron, vertices on the unit sphere, wound outward. */
struct Icosahedron
{
  std::array<std::array<Vec3, 3>, kIcosFaces> v;
  std::array<std::array<int, 3>, kIcosFaces>  idx;

  Icosahedron()
  {
    const float gr = 1.6180339887f;
    const float vert[kIcosVert][3] = { { 0, 1, gr }, { 0, -1, gr }, { 0, 1, -gr }, { 0, -1, -gr },
                                       { 1, gr, 0 }, { -1, gr, 0 }, { 1, -gr, 0 }, { -1, -gr, 0 },
                                       { gr, 0, 1 }, { -gr, 0, 1 }, { gr, 0, -1 }, { -gr, 0, -1 } };
    const int   faces[kIcosFaces][3] = { { 0, 1, 8 },  { 0, 8, 4 },  { 0, 4, 5 },  { 0, 5, 9 },  { 0, 9, 1 },
                                         { 1, 6, 8 },  { 8, 6, 10 }, { 8, 10, 4 }, { 4, 10, 2 }, { 4, 2, 5 },
                                         { 5, 2, 11 }, { 5, 11, 9 }, { 9, 11, 7 }, { 9, 7, 1 },  { 1, 7, 6 },
                                         { 3, 6, 7 },  { 3, 7, 11 }, { 3, 11, 2 }, { 3, 2, 10 }, { 3, 10, 6 } };
    for (int i = 0; i < kIcosFaces; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        idx[i][j] = faces[i][j];
        v[i][j] = Vec3(vert[faces[i][j]][0], vert[faces[i][j]][1], vert[faces[i][j]][2]).normalized();
      }
      if (((v[i][2] - v[i][1]).cross(v[i][1] - v[i][0])).dot(v[i][0]) < 0)
      {
        std::swap(v[i][0], v[i][1]);
        std::swap(idx[i][0], idx[i][1]);
      }
    }
  }

  /** Face hit by the ray along g, with barycentric weights (Moller-Trumbore). */
  bool
  Bin(const Vec3 & g, int & face, Vec3 & bary) const
  {
    if (g.squaredNorm() < kBaryEps)
    {
      return false;
    }
    for (int i = 0; i < kIcosFaces; ++i)
    {
      const Vec3  e1 = v[i][1] - v[i][0];
      const Vec3  e2 = v[i][2] - v[i][0];
      const Vec3  p = g.cross(e2);
      const float det = e1.dot(p);
      if (std::fabs(det) < kBaryEps)
      {
        continue;
      }
      const float inv = 1.0f / det;
      const Vec3  t = -v[i][0];
      const Vec3  q = t.cross(e1);
      const float by = inv * t.dot(p);
      const float bz = inv * g.dot(q);
      const float bx = 1.0f - by - bz;
      const float k = inv * e2.dot(q);
      if (bx < -kBaryEps || by < -kBaryEps || bz < -kBaryEps || k < 0)
      {
        continue;
      }
      face = i;
      bary = Vec3(bx, by, bz);
      return true;
    }
    return false;
  }
};

struct Candidate
{
  int             o, s, x, y, z;
  double          sd;
  double          strength;
  Eigen::Matrix3f R;
};

template <typename Fn>
void
ForSphere(const Volume & im, const Vec3 & c, float rad, float unit, Fn && fn)
{
  const int x0 = std::max(static_cast<int>(std::floor(c[0] - rad / unit)), 1);
  const int x1 = std::min(static_cast<int>(std::ceil(c[0] + rad / unit)), im.nx - 2);
  const int y0 = std::max(static_cast<int>(std::floor(c[1] - rad / unit)), 1);
  const int y1 = std::min(static_cast<int>(std::ceil(c[1] + rad / unit)), im.ny - 2);
  const int z0 = std::max(static_cast<int>(std::floor(c[2] - rad / unit)), 1);
  const int z1 = std::min(static_cast<int>(std::ceil(c[2] + rad / unit)), im.nz - 2);
  for (int z = z0; z <= z1; ++z)
  {
    for (int y = y0; y <= y1; ++y)
    {
      for (int x = x0; x <= x1; ++x)
      {
        const Vec3  disp = (Vec3(x, y, z) - c) * unit;
        const float sq = disp.squaredNorm();
        if (sq <= rad * rad)
        {
          fn(x, y, z, disp, sq);
        }
      }
    }
  }
}

/** Orientation from the windowed structure tensor; false rejects the keypoint. */
bool
AssignOrientation(const Volume & im, const Vec3 & c, double sigma, float unit, double cornerThresh, Eigen::Matrix3f & R)
{
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  Eigen::Vector3d win = Eigen::Vector3d::Zero();
  ForSphere(im, c, static_cast<float>(sigma * kOriRadFctr), unit, [&](int x, int y, int z, const Vec3 &, float sq) {
    const double          w = std::exp(-0.5 * sq / (sigma * sigma));
    const Eigen::Vector3d g = im.Gradient(x, y, z, unit).cast<double>();
    A += w * g * g.transpose();
    win += w * g;
  });
  if (win.squaredNorm() < kOriGradThresh)
  {
    return false;
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(A);
  const Eigen::Vector3d &                              L = eig.eigenvalues(); // ascending
  for (int i = 0; i < 2; ++i)
  {
    if (!(std::fabs(L[i + 1]) > 0.0) || std::fabs(L[i] / L[i + 1]) > kMaxEigRatio)
    {
      return false;
    }
  }
  double                         corner = DBL_MAX;
  std::array<Eigen::Vector3d, 2> axes;
  for (int i = 0; i < 2; ++i)
  {
    Eigen::Vector3d vr = eig.eigenvectors().col(2 - i);
    const double    d = win.dot(vr);
    corner = std::min(corner, std::fabs(d / (vr.norm() * win.norm())));
    axes[i] = d > 0.0 ? vr : -vr;
  }
  if (corner < cornerThresh)
  {
    return false;
  }
  R.col(0) = axes[0].cast<float>();
  R.col(1) = axes[1].cast<float>();
  R.col(2) = axes[0].cross(axes[1]).cast<float>();
  return true;
}

void
ExtractDescriptor(const Volume & im, const Candidate & key, float unit, const Icosahedron & ico, float * out)
{
  const float           sigma = static_cast<float>(key.sd * kDescSigFctr);
  const float           rad = static_cast<float>(kDescRadFctr) * sigma;
  const float           halfWidth = rad / std::sqrt(2.0f);
  const float           binFactor = kHistPerDim / (2.0f * halfWidth);
  const Eigen::Matrix3f Rt = key.R.transpose();
  std::fill(out, out + SIFT3DFeatures::DescriptorLength, 0.0f);

  ForSphere(im, Vec3(key.x, key.y, key.z), rad, unit, [&](int x, int y, int z, const Vec3 & disp, float sq) {
    const Vec3 bins = (Rt * disp + Vec3::Constant(halfWidth)) * binFactor;
    if ((bins.array() < 0.0f).any() || (bins.array() >= static_cast<float>(kHistPerDim)).any())
    {
      return;
    }
    const Vec3 grad = Rt * (im.Gradient(x, y, z, unit) * std::exp(-0.5f * sq / (sigma * sigma)));
    int        face;
    Vec3       bary;
    if (!ico.Bin(grad, face, bary))
    {
      return;
    }
    const float mag = grad.norm();
    const Vec3  frac = bins - bins.array().floor().matrix();
    for (int dz = 0; dz < 2; ++dz)
    {
      for (int dy = 0; dy < 2; ++dy)
      {
        for (int dx = 0; dx < 2; ++dx)
        {
          const int hx = static_cast<int>(bins[0]) + dx;
          const int hy = static_cast<int>(bins[1]) + dy;
          const int hz = static_cast<int>(bins[2]) + dz;
          if (hx >= kHistPerDim || hy >= kHistPerDim || hz >= kHistPerDim)
          {
            continue;
          }
          const float w =
            (dx ? frac[0] : 1.0f - frac[0]) * (dy ? frac[1] : 1.0f - frac[1]) * (dz ? frac[2] : 1.0f - frac[2]);
          float * hist = out + (hx + kHistPerDim * (hy + kHistPerDim * hz)) * kIcosVert;
          for (int j = 0; j < 3; ++j)
          {
            hist[ico.idx[face][j]] += mag * w * bary[j];
          }
        }
      }
    }
  });

  auto normalize = [out]() {
    double n = 0.0;
    for (unsigned int i = 0; i < SIFT3DFeatures::DescriptorLength; ++i)
    {
      n += static_cast<double>(out[i]) * out[i];
    }
    const float inv = static_cast<float>(1.0 / (std::sqrt(n) + DBL_EPSILON));
    for (unsigned int i = 0; i < SIFT3DFeatures::DescriptorLength; ++i)
    {
      out[i] *= inv;
    }
  };
  normalize();
  for (unsigned int i = 0; i < SIFT3DFeatures::DescriptorLength; ++i)
  {
    out[i] = std::min(out[i], kTruncThresh);
  }
  normalize();
}
} // namespace

SIFT3DFeatures
SIFT3DFeatureDetector::Compute(const ImageType * image) const
{
  // Canonical grid: identity direction, isotropic spacing, covering the input.
  const auto &         region = image->GetLargestPossibleRegion();
  ImageType::PointType lo, hi;
  lo.Fill(std::numeric_limits<double>::max());
  hi.Fill(std::numeric_limits<double>::lowest());
  for (unsigned int corner = 0; corner < 8; ++corner)
  {
    itk::ContinuousIndex<double, 3> ci;
    for (unsigned int d = 0; d < 3; ++d)
    {
      ci[d] = region.GetIndex(d) - 0.5 + ((corner >> d) & 1U ? region.GetSize(d) : 0);
    }
    ImageType::PointType p;
    image->TransformContinuousIndexToPhysicalPoint(ci, p);
    for (unsigned int d = 0; d < 3; ++d)
    {
      lo[d] = std::min(lo[d], p[d]);
      hi[d] = std::max(hi[d], p[d]);
    }
  }
  const double           h = m_SamplingSpacing;
  ImageType::SizeType    size;
  ImageType::PointType   origin;
  ImageType::SpacingType spacing;
  spacing.Fill(h);
  for (unsigned int d = 0; d < 3; ++d)
  {
    size[d] = static_cast<itk::SizeValueType>(std::floor((hi[d] - lo[d]) / h));
    origin[d] = lo[d] + 0.5 * h;
  }
  using ResampleType = itk::ResampleImageFilter<ImageType, ImageType>;
  auto resample = ResampleType::New();
  resample->SetInput(image);
  resample->SetTransform(itk::IdentityTransform<double, 3>::New());
  resample->SetInterpolator(itk::LinearInterpolateImageFunction<ImageType, double>::New());
  resample->SetOutputOrigin(origin);
  resample->SetOutputSpacing(spacing);
  resample->SetSize(size);
  resample->SetDefaultPixelValue(0.0f);
  resample->Update();
  const ImageType * grid = resample->GetOutput();

  Volume base;
  base.Resize(size[0], size[1], size[2]);
  std::copy(grid->GetBufferPointer(), grid->GetBufferPointer() + base.v.size(), base.v.begin());
  float maxAbs = 0.0f;
  for (const float x : base.v)
  {
    maxAbs = std::max(maxAbs, std::fabs(x));
  }
  if (maxAbs > 0.0f)
  {
    for (auto & x : base.v)
    {
      x /= maxAbs;
    }
  }

  const int minDim = std::min({ base.nx, base.ny, base.nz });
  const int numOctaves = static_cast<int>(std::log2(static_cast<double>(minDim))) - 3 + 1;
  if (numOctaves < 1)
  {
    itkExceptionMacro("image too small: need at least 8 voxels in each dimension at " << h << " mm");
  }
  const int n = static_cast<int>(m_NumberOfKeypointLevels);
  const int numGauss = n + 3; // levels s = -1 .. n + 1
  auto      sigmaAt = [&](int o, int s) { return m_Sigma0 * std::pow(2.0, o + static_cast<double>(s) / n); };
  if (sigmaAt(0, -1) <= m_SigmaN)
  {
    itkExceptionMacro("SigmaN must be smaller than the first pyramid scale " << sigmaAt(0, -1));
  }

  std::vector<std::vector<Volume>> gpyr(numOctaves, std::vector<Volume>(numGauss));
  std::vector<std::vector<Volume>> dog(numOctaves, std::vector<Volume>(numGauss - 1));
  for (int o = 0; o < numOctaves; ++o)
  {
    if (o == 0)
    {
      GaussianBlur(base, gpyr[0][0], std::sqrt(sigmaAt(0, -1) * sigmaAt(0, -1) - m_SigmaN * m_SigmaN));
    }
    else
    {
      Downsample2x(gpyr[o - 1][n], gpyr[o][0]); // level s = n - 1 has the scale of s = -1 one octave up
    }
    for (int k = 1; k < numGauss; ++k)
    {
      const double a = sigmaAt(0, k - 2);
      const double b = sigmaAt(0, k - 1);
      GaussianBlur(gpyr[o][k - 1], gpyr[o][k], std::sqrt(b * b - a * a));
    }
    for (int k = 0; k + 1 < numGauss; ++k)
    {
      Volume & d = dog[o][k];
      d.Resize(gpyr[o][k].nx, gpyr[o][k].ny, gpyr[o][k].nz);
      for (size_t i = 0; i < d.v.size(); ++i)
      {
        d.v[i] = gpyr[o][k].v[i] - gpyr[o][k + 1].v[i];
      }
    }
  }

  // DoG extrema: 6-neighbors in the level and the same voxel one level up and down.
  std::vector<Candidate> cand;
  std::mutex             mtx;
  for (int o = 0; o < numOctaves; ++o)
  {
    for (int k = 1; k <= n; ++k)
    {
      const Volume & prev = dog[o][k - 1];
      const Volume & cur = dog[o][k];
      const Volume & next = dog[o][k + 1];
      float          dmax = 0.0f;
      for (const float x : cur.v)
      {
        dmax = std::max(dmax, std::fabs(x));
      }
      const float thresh = static_cast<float>(m_PeakThreshold) * dmax;
      if (cur.nx < 3 || cur.ny < 3 || cur.nz < 3)
      {
        continue;
      }
      ParallelOverSlices(cur.nz - 2, [&](int zi) {
        const int              z = zi + 1;
        std::vector<Candidate> local;
        for (int y = 1; y < cur.ny - 1; ++y)
        {
          for (int x = 1; x < cur.nx - 1; ++x)
          {
            const float p = cur(x, y, z);
            if (!(p > thresh || p < -thresh))
            {
              continue;
            }
            const float nb[8] = { prev(x, y, z),    next(x, y, z),    cur(x + 1, y, z), cur(x - 1, y, z),
                                  cur(x, y + 1, z), cur(x, y - 1, z), cur(x, y, z + 1), cur(x, y, z - 1) };
            const bool  isMax = std::all_of(nb, nb + 8, [p](float q) { return p > q; });
            const bool  isMin = std::all_of(nb, nb + 8, [p](float q) { return p < q; });
            if (isMax || isMin)
            {
              local.push_back({ o, k - 1, x, y, z, sigmaAt(o, k - 1), std::fabs(p), Eigen::Matrix3f::Identity() });
            }
          }
        }
        const std::lock_guard<std::mutex> lock(mtx);
        cand.insert(cand.end(), local.begin(), local.end());
      });
    }
  }
  std::sort(cand.begin(), cand.end(), [](const Candidate & a, const Candidate & b) {
    return std::tie(a.o, a.s, a.z, a.y, a.x) < std::tie(b.o, b.s, b.z, b.y, b.x);
  });

  // Orientation; rejected candidates are dropped.
  std::vector<char> keep(cand.size(), 0);
  itk::MultiThreaderBase::New()->ParallelizeArray(
    0,
    cand.size(),
    [&](itk::SizeValueType i) {
      Candidate & c = cand[i];
      const float unit = static_cast<float>(1 << c.o);
      keep[i] =
        AssignOrientation(gpyr[c.o][c.s + 1], Vec3(c.x, c.y, c.z), kOriSigFctr * c.sd, unit, m_CornerThreshold, c.R);
    },
    nullptr);
  std::vector<Candidate> keys;
  for (size_t i = 0; i < cand.size(); ++i)
  {
    if (keep[i])
    {
      keys.push_back(cand[i]);
    }
  }
  if (m_MaximumNumberOfKeypoints > 0 && keys.size() > m_MaximumNumberOfKeypoints)
  {
    std::stable_sort(
      keys.begin(), keys.end(), [](const Candidate & a, const Candidate & b) { return a.strength > b.strength; });
    keys.resize(m_MaximumNumberOfKeypoints);
  }

  static const Icosahedron ico;
  SIFT3DFeatures           out;
  out.locations.resize(keys.size());
  out.scales.resize(keys.size());
  out.strengths.resize(keys.size());
  out.descriptors.resize(keys.size() * SIFT3DFeatures::DescriptorLength);
  itk::MultiThreaderBase::New()->ParallelizeArray(
    0,
    keys.size(),
    [&](itk::SizeValueType i) {
      const Candidate & c = keys[i];
      const float       unit = static_cast<float>(1 << c.o);
      ExtractDescriptor(
        gpyr[c.o][c.s + 1], c, unit, ico, out.descriptors.data() + i * SIFT3DFeatures::DescriptorLength);
      itk::ContinuousIndex<double, 3> ci;
      ci[0] = c.x * unit;
      ci[1] = c.y * unit;
      ci[2] = c.z * unit;
      grid->TransformContinuousIndexToPhysicalPoint(ci, out.locations[i]);
      out.scales[i] = c.sd * h;
      out.strengths[i] = c.strength;
    },
    nullptr);
  return out;
}

void
SIFT3DFeatureDetector::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "SamplingSpacing: " << m_SamplingSpacing << '\n'
     << indent << "PeakThreshold: " << m_PeakThreshold << '\n'
     << indent << "CornerThreshold: " << m_CornerThreshold << '\n'
     << indent << "NumberOfKeypointLevels: " << m_NumberOfKeypointLevels << '\n'
     << indent << "Sigma0: " << m_Sigma0 << '\n'
     << indent << "SigmaN: " << m_SigmaN << '\n'
     << indent << "MaximumNumberOfKeypoints: " << m_MaximumNumberOfKeypoints << '\n';
}

} // namespace brains
