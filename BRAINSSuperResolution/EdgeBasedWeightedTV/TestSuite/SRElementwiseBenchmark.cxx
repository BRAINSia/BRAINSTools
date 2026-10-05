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
/// \file SRElementwiseBenchmark.cxx
/// \brief Times a plain loop and the Eigen formulation of out = c*(a*x + y) and a full OpWeightedL2 run.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <vector>

#include "itk_eigen.h"
#include ITK_EIGEN(Core)

#include "itkGaussianImageSource.h"
#include "itkImageBufferRange.h"
#include "itkTimeProbe.h"
#include "OpWeightedL2.h"
#include "SRElementwiseOps.h"

namespace
{
double
MedianSeconds(std::vector<double> & t)
{
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}

void
BenchmarkAxpyScale(const size_t N, const int reps)
{
  std::mt19937                          gen(42);
  std::uniform_real_distribution<float> dist(-1.0F, 1.0F);
  std::vector<float>                    x(N), y0(N), yLoop(N), yEigen(N);
  for (size_t i = 0; i < N; ++i)
  {
    x[i] = dist(gen);
    y0[i] = dist(gen);
  }
  const float         a = 0.7F;
  const float         c = 1.3F;
  std::vector<double> tLoop, tEigen;
  for (int r = 0; r < reps; ++r)
  {
    yLoop = y0;
    const auto t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < N; ++i)
    {
      yLoop[i] = c * (a * x[i] + yLoop[i]);
    }
    tLoop.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());

    yEigen = y0;
    const auto                        t1 = std::chrono::steady_clock::now();
    Eigen::Map<Eigen::VectorXf>       yv(yEigen.data(), static_cast<Eigen::Index>(N));
    Eigen::Map<const Eigen::VectorXf> xv(x.data(), static_cast<Eigen::Index>(N));
    yv = c * (a * xv + yv);
    tEigen.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count());
  }
  float maxDiff = 0.0F;
  for (size_t i = 0; i < N; ++i)
  {
    maxDiff = std::max(maxDiff, std::abs(yLoop[i] - yEigen[i]));
  }
  const double mb = MedianSeconds(tLoop);
  const double me = MedianSeconds(tEigen);
  std::printf("N=%10zu floats (%6.1f MB)  %s %.4f ms  eigen %.4f ms  ratio %.2f  max|diff| %.3g\n",
              N,
              N * 4.0 / 1e6,
              "loop ",
              mb * 1e3,
              me * 1e3,
              me / mb,
              maxDiff);
}

FloatImageType::Pointer
MakeBlob(unsigned int side, double sigma)
{
  auto                     source = itk::GaussianImageSource<FloatImageType>::New();
  FloatImageType::SizeType size;
  size.Fill(side);
  source->SetSize(size);
  FloatImageType::SpacingType spacing;
  spacing.Fill(1.0);
  source->SetSpacing(spacing);
  itk::FixedArray<double, 3> sig;
  sig.Fill(sigma);
  source->SetSigma(sig);
  itk::FixedArray<double, 3> mean;
  mean.Fill(side / 2.0);
  source->SetMean(mean);
  source->SetScale(1.0);
  source->SetNormalized(false);
  source->Update();
  return source->GetOutput();
}
} // namespace

int
main(int argc, char * argv[])
{
  const int reps = (argc > 1) ? std::atoi(argv[1]) : 9;
  for (const size_t side : { 64, 128, 192 })
  {
    BenchmarkAxpyScale(3 * side * side * side, reps);
  }

  const unsigned int lowSide = (argc > 2) ? std::atoi(argv[2]) : 24;
  auto               lowres = MakeBlob(lowSide, lowSide / 4.0);
  auto               edge = MakeBlob(2 * lowSide, lowSide / 2.0);
  itk::TimeProbe     probe;
  probe.Start();
  FloatImageType::Pointer result = OpWeightedL2(lowres, edge);
  probe.Stop();
  double sum = 0.0;
  double sumSq = 0.0;
  for (const float v : itk::ImageBufferRange<FloatImageType>(*result))
  {
    sum += v;
    sumSq += static_cast<double>(v) * v;
  }
  std::printf("OpWeightedL2 %ux%ux%u -> %ux%ux%u  %.3f s  sum %.9g  sumsq %.9g\n",
              lowSide,
              lowSide,
              lowSide,
              2 * lowSide,
              2 * lowSide,
              2 * lowSide,
              probe.GetTotal(),
              sum,
              sumSq);
  return 0;
}
