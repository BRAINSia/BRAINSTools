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
/// \file SRElementwiseOpsGTest.cxx
/// \brief Pins the semantics of the SRElementwiseOps buffer helpers.

#include <gtest/gtest.h>
#include <cmath>

#include "itkGaussianImageSource.h"
#include "itkImageBufferRange.h"
#include "OpWeightedL2.h"
#include "SRElementwiseOps.h"

namespace
{
constexpr double kReferenceSum = 6987.67037; // 16^3 blob -> 32^3, same to 1e-6 with the former BLAS implementation
template <typename TImage>
typename TImage::Pointer
MakeImage(const typename TImage::PixelType & value, unsigned int side = 4)
{
  auto                      image = TImage::New();
  typename TImage::SizeType size;
  size.Fill(side);
  typename TImage::IndexType start;
  start.Fill(0);
  image->SetRegions(typename TImage::RegionType(start, size));
  image->Allocate();
  image->FillBuffer(value);
  return image;
}
} // namespace

TEST(SRElementwiseOps, AddAllElementsScalesAxpyInPlace)
{
  auto                    x = MakeImage<FloatImageType>(2.0F);
  auto                    y = MakeImage<FloatImageType>(5.0F);
  FloatImageType::Pointer out = y;
  AddAllElements(out, 3.0F, x, y, 0.5F); // out = 0.5 * (3*2 + 5)
  EXPECT_EQ(out.GetPointer(), y.GetPointer());
  for (const float v : itk::ImageBufferRange<FloatImageType>(*out))
  {
    EXPECT_FLOAT_EQ(v, 5.5F);
  }
  for (const float v : itk::ImageBufferRange<FloatImageType>(*x))
  {
    EXPECT_FLOAT_EQ(v, 2.0F);
  }
}

TEST(SRElementwiseOps, AddAllElementsUnitScaleSkipsScaling)
{
  auto                    x = MakeImage<FloatImageType>(-1.5F);
  auto                    y = MakeImage<FloatImageType>(4.0F);
  FloatImageType::Pointer out = y;
  AddAllElements(out, 2.0F, x, y); // out = 2*(-1.5) + 4
  for (const float v : itk::ImageBufferRange<FloatImageType>(*out))
  {
    EXPECT_FLOAT_EQ(v, 1.0F);
  }
}

TEST(SRElementwiseOps, AddAllElementsSwapsDistinctOutput)
{
  auto                         x = MakeImage<FloatImageType>(1.0F);
  auto                         y = MakeImage<FloatImageType>(1.0F);
  auto                         out = MakeImage<FloatImageType>(0.0F);
  const FloatImageType * const originalY = y.GetPointer();
  const FloatImageType * const originalOut = out.GetPointer();
  AddAllElements(out, 1.0F, x, y, 2.0F); // 2*(1+1) lands in the buffer formerly held by y
  EXPECT_EQ(out.GetPointer(), originalY);
  EXPECT_EQ(y.GetPointer(), originalOut);
  for (const float v : itk::ImageBufferRange<FloatImageType>(*out))
  {
    EXPECT_FLOAT_EQ(v, 4.0F);
  }
}

TEST(SRElementwiseOps, AddAllElementsCoversAllCovariantComponents)
{
  CVType xv;
  xv[0] = 1.0F;
  xv[1] = 2.0F;
  xv[2] = 3.0F;
  CVType yv;
  yv.Fill(10.0F);
  auto                 x = MakeImage<CVImageType>(xv);
  auto                 y = MakeImage<CVImageType>(yv);
  CVImageType::Pointer out = y;
  AddAllElements(out, 2.0F, x, y, 0.5F); // 0.5*(2*x + 10)
  for (const CVType & v : itk::ImageBufferRange<CVImageType>(*out))
  {
    EXPECT_FLOAT_EQ(v[0], 6.0F);
    EXPECT_FLOAT_EQ(v[1], 7.0F);
    EXPECT_FLOAT_EQ(v[2], 8.0F);
  }
}

TEST(SRElementwiseOps, MultiplyVectorsIsElementwise)
{
  CVType xv;
  xv[0] = 1.0F;
  xv[1] = 2.0F;
  xv[2] = 3.0F;
  CVType yv;
  yv[0] = 10.0F;
  yv[1] = 20.0F;
  yv[2] = 30.0F;
  auto x = MakeImage<CVImageType>(xv);
  auto y = MakeImage<CVImageType>(yv);
  auto out = MakeImage<CVImageType>(CVType());
  MultiplyVectors(out, x, y);
  for (const CVType & v : itk::ImageBufferRange<CVImageType>(*out))
  {
    EXPECT_FLOAT_EQ(v[0], 10.0F);
    EXPECT_FLOAT_EQ(v[1], 40.0F);
    EXPECT_FLOAT_EQ(v[2], 90.0F);
  }
}

TEST(SRElementwiseOps, DuplicateCopiesBuffer)
{
  auto src = MakeImage<FloatImageType>(7.25F);
  auto dst = MakeImage<FloatImageType>(0.0F);
  Duplicate(src, dst);
  for (const float v : itk::ImageBufferRange<FloatImageType>(*dst))
  {
    EXPECT_FLOAT_EQ(v, 7.25F);
  }
}

TEST(OpWeightedL2, SyntheticBlobReproducesReference)
{
  auto makeBlob = [](unsigned int side, double sigma) {
    auto                     source = itk::GaussianImageSource<FloatImageType>::New();
    FloatImageType::SizeType size;
    size.Fill(side);
    source->SetSize(size);
    itk::FixedArray<double, 3> sig;
    sig.Fill(sigma);
    source->SetSigma(sig);
    itk::FixedArray<double, 3> mean;
    mean.Fill(side / 2.0);
    source->SetMean(mean);
    source->SetScale(1.0);
    source->SetNormalized(false);
    source->Update();
    return FloatImageType::Pointer(source->GetOutput());
  };
  constexpr unsigned int lowSide = 16;
  const auto             result = OpWeightedL2(makeBlob(lowSide, lowSide / 4.0), makeBlob(2 * lowSide, lowSide / 2.0));
  ASSERT_TRUE(result.IsNotNull());
  EXPECT_EQ(result->GetLargestPossibleRegion().GetSize()[0], 2 * lowSide);
  double sum = 0.0;
  for (const float v : itk::ImageBufferRange<FloatImageType>(*result))
  {
    ASSERT_TRUE(std::isfinite(v));
    sum += v;
  }
  EXPECT_NEAR(sum, kReferenceSum, 1e-3 * std::abs(kReferenceSum));
}
