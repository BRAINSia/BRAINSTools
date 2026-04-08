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
/// \file itkIOGTest.cxx
/// \brief GTest coverage for itkIO.h utilities (issue #403).
///
/// Tests ReadImage, WriteImage, CopyImage, AllocateImageFromExample,
/// and TypeCast with non-trivial pixel values and image metadata.

#include <gtest/gtest.h>
#include <filesystem>
#include <string>

#include "itkIO.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"

// ---------------------------------------------------------------------------
// Helper: build a 5x6x7 float image with non-trivial origin/spacing/pixels
// ---------------------------------------------------------------------------
template <typename TPixel>
static typename itk::Image<TPixel, 3>::Pointer
MakeTestImage(std::function<TPixel(unsigned int, unsigned int, unsigned int)> pixelFn)
{
  using ImageType = itk::Image<TPixel, 3>;
  auto image = ImageType::New();

  typename ImageType::IndexType start;
  start.Fill(0);
  typename ImageType::SizeType size;
  size[0] = 5;
  size[1] = 6;
  size[2] = 7;
  typename ImageType::RegionType region;
  region.SetIndex(start);
  region.SetSize(size);

  typename ImageType::PointType origin;
  origin[0] = -12.5;
  origin[1] = 33.7;
  origin[2] = 8.1;

  typename ImageType::SpacingType spacing;
  spacing[0] = 0.75;
  spacing[1] = 1.25;
  spacing[2] = 2.0;

  image->SetRegions(region);
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
  image->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> it(image, region);
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    const auto idx = it.GetIndex();
    it.Set(pixelFn(idx[0], idx[1], idx[2]));
  }
  return image;
}

// Helper: unique temp path for each test
static std::string
UniqueTempPath(const std::string & suffix)
{
  auto p = std::filesystem::temp_directory_path() / ("itkIOGTest_" + suffix);
  return p.string();
}

// ---------------------------------------------------------------------------
// Float write/read roundtrip
// ---------------------------------------------------------------------------
TEST(itkIO, FloatWriteReadRoundtrip)
{
  using ImageType = itk::Image<float, 3>;
  const std::string path = UniqueTempPath("float.nrrd");

  auto image = MakeTestImage<float>([](unsigned int x, unsigned int y, unsigned int z) -> float {
    return 42.5f * static_cast<float>(x) - 17.3f * static_cast<float>(y) + 0.001f * static_cast<float>(z) + 100.0f;
  });

  itkUtil::WriteImage<ImageType>(image, path);
  auto loaded = itkUtil::ReadImage<ImageType>(path);

  // Verify origin and spacing
  for (unsigned int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(loaded->GetOrigin()[d], image->GetOrigin()[d], 1e-4) << "origin mismatch dim " << d;
    EXPECT_NEAR(loaded->GetSpacing()[d], image->GetSpacing()[d], 1e-4) << "spacing mismatch dim " << d;
  }

  // Verify ALL pixels
  itk::ImageRegionIteratorWithIndex<ImageType> itOrig(image, image->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> itLoad(loaded, loaded->GetLargestPossibleRegion());
  for (itOrig.GoToBegin(), itLoad.GoToBegin(); !itOrig.IsAtEnd(); ++itOrig, ++itLoad)
  {
    EXPECT_FLOAT_EQ(itLoad.Get(), itOrig.Get()) << "pixel mismatch at " << itOrig.GetIndex();
  }

  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// Double write/read roundtrip
// ---------------------------------------------------------------------------
TEST(itkIO, DoubleWriteReadRoundtrip)
{
  using ImageType = itk::Image<double, 3>;
  const std::string path = UniqueTempPath("double.nrrd");

  auto image = MakeTestImage<double>([](unsigned int x, unsigned int y, unsigned int z) -> double {
    return 42.5 * static_cast<double>(x) - 17.3 * static_cast<double>(y) + 0.001 * static_cast<double>(z) + 100.0;
  });

  itkUtil::WriteImage<ImageType>(image, path);
  auto loaded = itkUtil::ReadImage<ImageType>(path);

  for (unsigned int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(loaded->GetOrigin()[d], image->GetOrigin()[d], 1e-4);
    EXPECT_NEAR(loaded->GetSpacing()[d], image->GetSpacing()[d], 1e-4);
  }

  itk::ImageRegionIteratorWithIndex<ImageType> itOrig(image, image->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> itLoad(loaded, loaded->GetLargestPossibleRegion());
  for (itOrig.GoToBegin(), itLoad.GoToBegin(); !itOrig.IsAtEnd(); ++itOrig, ++itLoad)
  {
    EXPECT_DOUBLE_EQ(itLoad.Get(), itOrig.Get()) << "pixel mismatch at " << itOrig.GetIndex();
  }

  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// Short write/read roundtrip
// ---------------------------------------------------------------------------
TEST(itkIO, ShortWriteReadRoundtrip)
{
  using ImageType = itk::Image<short, 3>;
  const std::string path = UniqueTempPath("short.nrrd");

  auto image = MakeTestImage<short>([](unsigned int x, unsigned int y, unsigned int z) -> short {
    return static_cast<short>(-17 * static_cast<int>(x) + 200);
  });

  itkUtil::WriteImage<ImageType>(image, path);
  auto loaded = itkUtil::ReadImage<ImageType>(path);

  itk::ImageRegionIteratorWithIndex<ImageType> itOrig(image, image->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> itLoad(loaded, loaded->GetLargestPossibleRegion());
  for (itOrig.GoToBegin(), itLoad.GoToBegin(); !itOrig.IsAtEnd(); ++itOrig, ++itLoad)
  {
    EXPECT_EQ(itLoad.Get(), itOrig.Get()) << "pixel mismatch at " << itOrig.GetIndex();
  }

  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// UChar write/read roundtrip
// ---------------------------------------------------------------------------
TEST(itkIO, UCharWriteReadRoundtrip)
{
  using ImageType = itk::Image<unsigned char, 3>;
  const std::string path = UniqueTempPath("uchar.nrrd");

  auto image = MakeTestImage<unsigned char>([](unsigned int x, unsigned int y, unsigned int z) -> unsigned char {
    return static_cast<unsigned char>((13 * x + 7 * y + 3 * z) % 256);
  });

  itkUtil::WriteImage<ImageType>(image, path);
  auto loaded = itkUtil::ReadImage<ImageType>(path);

  itk::ImageRegionIteratorWithIndex<ImageType> itOrig(image, image->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> itLoad(loaded, loaded->GetLargestPossibleRegion());
  for (itOrig.GoToBegin(), itLoad.GoToBegin(); !itOrig.IsAtEnd(); ++itOrig, ++itLoad)
  {
    EXPECT_EQ(itLoad.Get(), itOrig.Get()) << "pixel mismatch at " << itOrig.GetIndex();
  }

  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// Read non-existent file throws
// ---------------------------------------------------------------------------
TEST(itkIO, ReadNonExistentFileThrows)
{
  using ImageType = itk::Image<float, 3>;
  EXPECT_THROW(itkUtil::ReadImage<ImageType>("/this/path/does/not/exist.nrrd"), itk::ExceptionObject);
}

// ---------------------------------------------------------------------------
// CopyImage: pixels match and images are distinct objects
// ---------------------------------------------------------------------------
TEST(itkIO, CopyImagePixelsMatchAndDistinctPointers)
{
  using ImageType = itk::Image<float, 3>;

  auto original = MakeTestImage<float>([](unsigned int x, unsigned int y, unsigned int z) -> float {
    return 42.5f * static_cast<float>(x) - 17.3f * static_cast<float>(y) + 0.001f * static_cast<float>(z) + 100.0f;
  });

  auto copy = itkUtil::CopyImage<ImageType>(original);

  // Distinct objects
  EXPECT_NE(original.GetPointer(), copy.GetPointer());
  EXPECT_NE(original->GetBufferPointer(), copy->GetBufferPointer());

  // Same pixel values
  itk::ImageRegionIteratorWithIndex<ImageType> itOrig(original, original->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> itCopy(copy, copy->GetLargestPossibleRegion());
  for (itOrig.GoToBegin(), itCopy.GoToBegin(); !itOrig.IsAtEnd(); ++itOrig, ++itCopy)
  {
    EXPECT_FLOAT_EQ(itCopy.Get(), itOrig.Get());
  }
}

// ---------------------------------------------------------------------------
// AllocateImageFromExample: region/spacing/origin copied from template
// ---------------------------------------------------------------------------
TEST(itkIO, AllocateImageFromExampleCopiesMetadata)
{
  using ImageType = itk::Image<float, 3>;

  auto templateImage = MakeTestImage<float>([](unsigned int, unsigned int, unsigned int) -> float { return 0.0f; });

  auto allocated = itkUtil::AllocateImageFromExample<ImageType>(templateImage);

  // Verify region
  EXPECT_EQ(allocated->GetLargestPossibleRegion().GetSize(), templateImage->GetLargestPossibleRegion().GetSize());
  EXPECT_EQ(allocated->GetLargestPossibleRegion().GetIndex(), templateImage->GetLargestPossibleRegion().GetIndex());

  // Verify spacing and origin
  for (unsigned int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(allocated->GetSpacing()[d], templateImage->GetSpacing()[d], 1e-4);
    EXPECT_NEAR(allocated->GetOrigin()[d], templateImage->GetOrigin()[d], 1e-4);
  }

  // Verify direction
  EXPECT_EQ(allocated->GetDirection(), templateImage->GetDirection());
}

// ---------------------------------------------------------------------------
// TypeCast float->double: pixel values preserved
// ---------------------------------------------------------------------------
TEST(itkIO, TypeCastFloatToDouble)
{
  using FloatImageType = itk::Image<float, 3>;
  using DoubleImageType = itk::Image<double, 3>;

  auto floatImage = MakeTestImage<float>([](unsigned int x, unsigned int y, unsigned int z) -> float {
    return 42.5f * static_cast<float>(x) - 17.3f * static_cast<float>(y) + 0.001f * static_cast<float>(z) + 100.0f;
  });

  auto doubleImage = itkUtil::TypeCast<FloatImageType, DoubleImageType>(floatImage);

  itk::ImageRegionIteratorWithIndex<FloatImageType>  itF(floatImage, floatImage->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<DoubleImageType> itD(doubleImage, doubleImage->GetLargestPossibleRegion());
  for (itF.GoToBegin(), itD.GoToBegin(); !itF.IsAtEnd(); ++itF, ++itD)
  {
    EXPECT_NEAR(itD.Get(), static_cast<double>(itF.Get()), 1e-6) << "pixel mismatch at " << itF.GetIndex();
  }
}
