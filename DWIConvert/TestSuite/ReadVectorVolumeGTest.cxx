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
/// \file ReadVectorVolumeGTest.cxx
/// \brief ReadVectorVolume rescales a non-short vector image into the short range.

#include <gtest/gtest.h>
#include <string>

#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "DWIConvertUtils.h"

namespace
{
using FloatVectorImageType = itk::VectorImage<float, 3>;
using ShortVectorImageType = itk::VectorImage<short, 3>;

FloatVectorImageType::Pointer
MakeFloatImage(const float lowValue, const float highValue)
{
  auto img = FloatVectorImageType::New();
  img->SetVectorLength(2);
  FloatVectorImageType::SizeType size;
  size.Fill(2);
  img->SetRegions(size);
  FloatVectorImageType::SpacingType spacing;
  spacing[0] = 1.5;
  spacing[1] = 2.0;
  spacing[2] = 2.5;
  img->SetSpacing(spacing);
  img->Allocate();

  itk::ImageRegionIterator<FloatVectorImageType> it(img, img->GetLargestPossibleRegion());
  unsigned int                                   n = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++n)
  {
    FloatVectorImageType::PixelType px(2);
    px[0] = (n == 0) ? lowValue : highValue;
    px[1] = (n == 1) ? lowValue : highValue;
    it.Set(px);
  }
  return img;
}

std::string
WriteFloatImage(const FloatVectorImageType * img, const std::string & name)
{
  const std::string fname = std::string(::testing::TempDir()) + name;
  auto              writer = itk::ImageFileWriter<FloatVectorImageType>::New();
  writer->SetInput(img);
  writer->SetFileName(fname);
  writer->Update();
  return fname;
}
} // namespace

TEST(ReadVectorVolume, FloatInputIsRescaledToShortRange)
{
  const std::string fname = WriteFloatImage(MakeFloatImage(-10.0F, 30.0F), "ReadVectorVolumeGTest_float.nrrd");
  ShortVectorImageType::Pointer img;
  ASSERT_EQ(ReadVectorVolume<ShortVectorImageType>(img, fname, true), EXIT_SUCCESS);
  ASSERT_TRUE(img.IsNotNull());

  EXPECT_EQ(img->GetVectorLength(), 2U);
  EXPECT_EQ(img->GetLargestPossibleRegion().GetSize()[2], 2U);
  EXPECT_DOUBLE_EQ(img->GetSpacing()[0], 1.5);
  EXPECT_DOUBLE_EQ(img->GetSpacing()[1], 2.0);
  EXPECT_DOUBLE_EQ(img->GetSpacing()[2], 2.5);

  // Data range [-10, 30] maps to [0, 32766]: the minimum to 0, the maximum to 32766.
  const ShortVectorImageType::PixelType first = img->GetPixel({ 0, 0, 0 });
  const ShortVectorImageType::PixelType second = img->GetPixel({ 1, 0, 0 });
  EXPECT_EQ(first[0], 0);
  EXPECT_EQ(first[1], 32766);
  EXPECT_EQ(second[0], 32766);
  EXPECT_EQ(second[1], 0);
}

TEST(ReadVectorVolume, ConstantFloatInputIsFinite)
{
  const std::string fname = WriteFloatImage(MakeFloatImage(5.0F, 5.0F), "ReadVectorVolumeGTest_constant.nrrd");
  ShortVectorImageType::Pointer img;
  ASSERT_EQ(ReadVectorVolume<ShortVectorImageType>(img, fname, true), EXIT_SUCCESS);
  ASSERT_TRUE(img.IsNotNull());
  const ShortVectorImageType::PixelType px = img->GetPixel({ 0, 0, 0 });
  EXPECT_EQ(px[0], 0);
  EXPECT_EQ(px[1], 0);
}
