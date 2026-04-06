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
/// \file AverageImageFilterGTest.cxx
/// \brief GTest port of AverageImageFilterTest.cxx (issue #403 — GTest migration).
///
/// Verifies that itkAverageImageFilter produces the correct pixel-wise
/// arithmetic mean of its inputs.  The approach mirrors the original:
/// generate N random images, compute a manual average, subtract the filter
/// output, and check that the error statistics are within floating-point
/// round-off.

#include <gtest/gtest.h>

#include <itkImage.h>
#include <itkRandomImageSource.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkSubtractImageFilter.h>
#include <itkStatisticsImageFilter.h>

#include "itkAverageImageFilter.h"

TEST(AverageImageFilter, AverageMatchesManualComputation)
{
  const unsigned int numTestImages = 4;
  const unsigned int testImageDim = 4;

  using FloatImage2DType = itk::Image<float, 2>;
  using AverageImageFilterType = itk::AverageImageFilter<FloatImage2DType, FloatImage2DType>;
  using SubtractFilterType = itk::SubtractImageFilter<FloatImage2DType>;
  using StatFilterType = itk::StatisticsImageFilter<FloatImage2DType>;
  using FloatImageVector = std::vector<FloatImage2DType::Pointer>;
  using FloatImageConstIterator = itk::ImageRegionConstIterator<FloatImage2DType>;
  using FloatImageIterator = itk::ImageRegionIterator<FloatImage2DType>;
  using FloatImageConstIteratorVector = std::vector<FloatImageConstIterator>;

  FloatImage2DType::SizeType randomSize;
  randomSize[1] = randomSize[0] = testImageDim;

  FloatImageVector              inputImages;
  FloatImageConstIteratorVector inputIterators;

  for (unsigned int i = 0; i < numTestImages; ++i)
  {
    const auto random = itk::RandomImageSource<FloatImage2DType>::New();
    random->SetMin(0.0);
    random->SetMax(1000.0);
    random->SetSize(randomSize);
    random->Update();
    inputImages.emplace_back(random->GetOutput());
    inputIterators.emplace_back(inputImages[i], inputImages[i]->GetLargestPossibleRegion());
  }

  // Manual average — same algorithm as the original test
  const FloatImage2DType::Pointer testAvg = FloatImage2DType::New();
  testAvg->CopyInformation(inputImages[0]);
  testAvg->SetRegions(randomSize);
  testAvg->Allocate();

  FloatImageIterator avgIt(testAvg, testAvg->GetLargestPossibleRegion());
  for (avgIt.GoToBegin(); !avgIt.IsAtEnd(); ++avgIt)
  {
    double sum = 0.0;
    for (unsigned int i = 0; i < numTestImages; ++i)
    {
      sum += inputIterators[i].Value();
      ++inputIterators[i];
    }
    avgIt.Set(static_cast<float>(sum / numTestImages));
  }

  // Run the filter
  const auto avgFilter = AverageImageFilterType::New();
  for (unsigned int i = 0; i < numTestImages; ++i)
  {
    avgFilter->SetInput(i, inputImages[i]);
  }

  // Subtract manual average from filter output; stats should be ~zero
  const auto subFilter = SubtractFilterType::New();
  subFilter->SetInput1(testAvg);
  subFilter->SetInput2(avgFilter->GetOutput());

  const auto statFilter = StatFilterType::New();
  statFilter->SetInput(subFilter->GetOutput());
  statFilter->Update();

  const double eps = 1e-8;
  EXPECT_NEAR(statFilter->GetMinimum(), 0.0, eps);
  EXPECT_NEAR(statFilter->GetMaximum(), 0.0, eps);
  EXPECT_NEAR(statFilter->GetMean(), 0.0, eps);
}
