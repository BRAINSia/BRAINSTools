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

#include "itkMacro.h"
#include "itkMultiThreaderBase.h"
#include "itk_eigen.h"
#include ITK_EIGEN(Dense)
#include ITK_EIGEN(Geometry)

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <random>

namespace brains
{
namespace
{
using Points = Eigen::Matrix<double, 3, Eigen::Dynamic>;

struct Model
{
  Eigen::Matrix3d A{ Eigen::Matrix3d::Identity() };
  Eigen::Vector3d t{ Eigen::Vector3d::Zero() };
};

Model
FromHomogeneous(const Eigen::Matrix4d & T)
{
  Model m;
  m.A = T.topLeftCorner<3, 3>();
  m.t = T.topRightCorner<3, 1>();
  return m;
}

Model
Fit(TransformModel model, const Points & P, const Points & Q, const std::vector<int> & idx)
{
  const auto n = static_cast<Eigen::Index>(idx.size());
  Points     p(3, n);
  Points     q(3, n);
  for (Eigen::Index i = 0; i < n; ++i)
  {
    p.col(i) = P.col(idx[i]);
    q.col(i) = Q.col(idx[i]);
  }
  if (model == TransformModel::Affine && n >= 12)
  {
    Eigen::MatrixXd M(n, 4);
    M.leftCols<3>() = p.transpose();
    M.col(3).setOnes();
    const Eigen::MatrixXd X = M.colPivHouseholderQr().solve(Eigen::MatrixXd(q.transpose()));
    Model                 m;
    m.A = X.topRows<3>().transpose();
    m.t = X.row(3).transpose();
    return m;
  }
  return FromHomogeneous(Eigen::umeyama(p, q, model != TransformModel::Rigid));
}

Eigen::VectorXd
Residuals(const Model & m, const Points & P, const Points & Q)
{
  return ((m.A * P).colwise() + m.t - Q).colwise().norm().transpose();
}

/** Lower is better. */
double
Score(RobustEstimator estimator, const Eigen::VectorXd & r, double tol)
{
  if (estimator == RobustEstimator::RANSAC)
  {
    return -static_cast<double>((r.array() < tol).count());
  }
  return r.array().square().min(tol * tol).sum();
}

std::vector<int>
Inliers(const Eigen::VectorXd & r, double tol)
{
  std::vector<int> idx;
  for (Eigen::Index i = 0; i < r.size(); ++i)
  {
    if (r[i] < tol)
    {
      idx.push_back(static_cast<int>(i));
    }
  }
  return idx;
}

/** Refit on inliers until the inlier set stops changing. */
Model
Refit(TransformModel model, const Points & P, const Points & Q, Model m, double tol)
{
  std::vector<int> inl = Inliers(Residuals(m, P, Q), tol);
  for (int it = 0; it < 3 && inl.size() >= 3; ++it)
  {
    m = Fit(model, P, Q, inl);
    std::vector<int> next = Inliers(Residuals(m, P, Q), tol);
    if (next.size() < 4 || next == inl)
    {
      break;
    }
    inl = std::move(next);
  }
  return m;
}
} // namespace

KeypointMatches
MatchFeatures(const SIFT3DFeatures & fixed, const SIFT3DFeatures & moving, double maximumRatio)
{
  using RowMat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  constexpr auto  L = static_cast<Eigen::Index>(SIFT3DFeatures::DescriptorLength);
  const auto      nf = static_cast<Eigen::Index>(fixed.size());
  const auto      nm = static_cast<Eigen::Index>(moving.size());
  KeypointMatches out;
  if (nf < 2 || nm < 1)
  {
    return out;
  }
  const Eigen::Map<const RowMat> F(fixed.descriptors.data(), nf, L);
  const Eigen::Map<const RowMat> M(moving.descriptors.data(), nm, L);
  const Eigen::VectorXf          fn = F.rowwise().squaredNorm();

  std::vector<int>       best(nm, -1);
  std::vector<double>    ratio(nm, 1.0);
  constexpr Eigen::Index chunk = 256;
  const auto             numChunks = static_cast<itk::SizeValueType>((nm + chunk - 1) / chunk);
  itk::MultiThreaderBase::New()->ParallelizeArray(
    0,
    numChunks,
    [&](itk::SizeValueType c) {
      const Eigen::Index r0 = static_cast<Eigen::Index>(c) * chunk;
      const Eigen::Index len = std::min(chunk, nm - r0);
      const RowMat       D = (-2.0f * M.middleRows(r0, len) * F.transpose()).rowwise() + fn.transpose();
      for (Eigen::Index r = 0; r < len; ++r)
      {
        const float mn = M.row(r0 + r).squaredNorm();
        float       d1 = std::numeric_limits<float>::max();
        float       d2 = d1;
        int         j1 = -1;
        for (Eigen::Index j = 0; j < nf; ++j)
        {
          const float d = D(r, j) + mn;
          if (d < d1)
          {
            d2 = d1;
            d1 = d;
            j1 = static_cast<int>(j);
          }
          else if (d < d2)
          {
            d2 = d;
          }
        }
        const double rr = std::sqrt(std::max(d1, 0.0f) / std::max(d2, std::numeric_limits<float>::min()));
        if (rr < maximumRatio)
        {
          best[r0 + r] = j1;
          ratio[r0 + r] = rr;
        }
      }
    },
    nullptr);

  for (Eigen::Index i = 0; i < nm; ++i)
  {
    if (best[i] >= 0)
    {
      out.fixed.push_back(fixed.locations[best[i]]);
      out.moving.push_back(moving.locations[i]);
      out.ratio.push_back(ratio[i]);
    }
  }
  return out;
}

bool
ParseTransformModel(const std::string & name, TransformModel & model)
{
  if (name == "Rigid")
  {
    model = TransformModel::Rigid;
  }
  else if (name == "Similarity")
  {
    model = TransformModel::Similarity;
  }
  else if (name == "Affine")
  {
    model = TransformModel::Affine;
  }
  else
  {
    return false;
  }
  return true;
}

bool
ParseRobustEstimator(const std::string & name, RobustEstimator & estimator)
{
  if (name == "RANSAC")
  {
    estimator = RobustEstimator::RANSAC;
  }
  else if (name == "MSAC")
  {
    estimator = RobustEstimator::MSAC;
  }
  else if (name == "PROSAC")
  {
    estimator = RobustEstimator::PROSAC;
  }
  else
  {
    return false;
  }
  return true;
}

RobustFitResult
FitTransform(const KeypointMatches & matches, const RobustFitOptions & options)
{
  const auto n = static_cast<Eigen::Index>(matches.size());
  if (n < 3)
  {
    itkGenericExceptionMacro("need at least 3 keypoint matches, found " << n);
  }
  Points P(3, n);
  Points Q(3, n);
  for (Eigen::Index i = 0; i < n; ++i)
  {
    for (int d = 0; d < 3; ++d)
    {
      P(d, i) = matches.fixed[i][d];
      Q(d, i) = matches.moving[i][d];
    }
  }

  // PROSAC draws from the best-ratio prefix, which grows to all matches by half the iterations.
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  if (options.estimator == RobustEstimator::PROSAC)
  {
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return matches.ratio[a] < matches.ratio[b]; });
  }
  const double tol = options.inlierTolerance;
  const auto   scoring = options.estimator == RobustEstimator::RANSAC ? RobustEstimator::RANSAC : RobustEstimator::MSAC;

  std::mt19937 rng(options.seed);
  Model        best;
  double       bestScore = std::numeric_limits<double>::max();
  const double growth = std::max(1.0, options.iterations / 2.0);
  for (unsigned int it = 0; it < options.iterations; ++it)
  {
    Eigen::Index pool = n;
    if (options.estimator == RobustEstimator::PROSAC)
    {
      pool = std::min<Eigen::Index>(n, 3 + static_cast<Eigen::Index>((n - 3) * (it / growth)));
    }
    std::uniform_int_distribution<Eigen::Index> pick(0, pool - 1);
    std::array<int, 3>                          s{};
    s[0] = order[pick(rng)];
    do
    {
      s[1] = order[pick(rng)];
    } while (s[1] == s[0]);
    do
    {
      s[2] = order[pick(rng)];
    } while (s[2] == s[0] || s[2] == s[1]);
    const Eigen::Vector3d e1 = P.col(s[1]) - P.col(s[0]);
    const Eigen::Vector3d e2 = P.col(s[2]) - P.col(s[0]);
    if (e1.cross(e2).norm() < 1.0)
    {
      continue; // nearly collinear sample
    }
    Model        m = Fit(TransformModel::Similarity, P, Q, { s[0], s[1], s[2] });
    const double scale = m.A.col(0).norm();
    if (scale < options.minimumScale || scale > options.maximumScale)
    {
      continue;
    }
    double score = Score(scoring, Residuals(m, P, Q), tol);
    if (score < bestScore && options.localOptimization)
    {
      const Model  lo = Refit(options.model, P, Q, m, tol);
      const double loScore = Score(scoring, Residuals(lo, P, Q), tol);
      if (loScore < score)
      {
        m = lo;
        score = loScore;
      }
    }
    if (score < bestScore)
    {
      bestScore = score;
      best = m;
    }
  }
  if (bestScore == std::numeric_limits<double>::max())
  {
    itkGenericExceptionMacro("no valid hypothesis among " << n << " matches");
  }
  best = Refit(options.model, P, Q, best, tol);

  RobustFitResult       out;
  const Eigen::VectorXd r = Residuals(best, P, Q);
  out.inliers.resize(n);
  double sum = 0.0;
  for (Eigen::Index i = 0; i < n; ++i)
  {
    out.inliers[i] = r[i] < tol;
    if (out.inliers[i])
    {
      ++out.numberOfInliers;
      sum += r[i];
    }
  }
  out.meanInlierResidual = out.numberOfInliers ? sum / out.numberOfInliers : 0.0;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      out.matrix(i, j) = best.A(i, j);
    }
    out.translation[i] = best.t[i];
  }
  return out;
}

} // namespace brains
