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
/// \file Slicer3LandmarkIOGTest.cxx
/// \brief GTest port of Slicer3LandmarkIOTest.cxx (issue #403 — GTest migration).
///
/// Verifies that ReadSlicer3toITKLmk() and ReadLandmarkWeights() throw
/// itk::ExceptionObject when given non-existent file paths.  The original
/// test used manual try/catch; GTest's EXPECT_THROW makes the intent explicit
/// and produces named CTest entries.
///
/// Extended with write/read roundtrip tests for V3/V4 FCSV, negative
/// coordinate preservation, empty/comment-only file rejection, and
/// landmark weight IO.

#include <gtest/gtest.h>

#include "Slicer3LandmarkIO.h"

#include <filesystem>
#include <fstream>
#include <string>

static constexpr double kCoordTol = 1e-4;

// Helper: unique temp file path
static std::string
UniqueTempPath(const std::string & suffix)
{
  auto p = std::filesystem::temp_directory_path() / ("Slicer3LmkGTest_" + suffix);
  return p.string();
}

// ---------------------------------------------------------------------------
// Original exception tests
// ---------------------------------------------------------------------------

TEST(Slicer3LandmarkIO, ThrowsOnNonExistentLandmarkFile)
{
  // Expect an itk::ExceptionObject (or any derived class) when the file
  // path does not exist.  This mirrors the original catch block that
  // treated the exception as proof of correct error handling.
  EXPECT_THROW(ReadSlicer3toITKLmk("/this/file/does/not/exist.fcsv"), itk::ExceptionObject);
}

TEST(Slicer3LandmarkIO, ThrowsOnNonExistentLandmarkWeightFile)
{
  EXPECT_THROW(ReadLandmarkWeights("/this/file/does/not/exist/either.fcsv"), itk::ExceptionObject);
}

// ---------------------------------------------------------------------------
// V4 FCSV write/read roundtrip with non-trivial LPS coordinates
// ---------------------------------------------------------------------------
TEST(Slicer3LandmarkIO, V4FcsvWriteReadRoundtrip)
{
  const std::string path = UniqueTempPath("v4_roundtrip.fcsv");

  LandmarksMapType  landmarks;
  LandmarkPointType ac;
  ac[0] = -127.456;
  ac[1] = 42.5;
  ac[2] = -0.001;
  landmarks["AC"] = ac;

  LandmarkPointType pc;
  pc[0] = 88.123;
  pc[1] = -33.7;
  pc[2] = 255.0;
  landmarks["PC"] = pc;

  WriteITKtoSlicer3Lmk(path, landmarks, SLICER_V4_FCSV);
  const LandmarksMapType loaded = ReadSlicer3toITKLmk(path);

  ASSERT_EQ(loaded.size(), landmarks.size());
  for (const auto & [name, pt] : landmarks)
  {
    auto it = loaded.find(name);
    ASSERT_NE(it, loaded.end()) << "landmark " << name << " not found in loaded file";
    for (unsigned int d = 0; d < 3; ++d)
    {
      EXPECT_NEAR(it->second[d], pt[d], kCoordTol) << "landmark " << name << " component " << d;
    }
  }

  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// V3 FCSV write/read roundtrip
// ---------------------------------------------------------------------------
TEST(Slicer3LandmarkIO, V3FcsvWriteReadRoundtrip)
{
  const std::string path = UniqueTempPath("v3_roundtrip.fcsv");

  LandmarksMapType  landmarks;
  LandmarkPointType ac;
  ac[0] = -127.456;
  ac[1] = 42.5;
  ac[2] = -0.001;
  landmarks["AC"] = ac;

  LandmarkPointType pc;
  pc[0] = 88.123;
  pc[1] = -33.7;
  pc[2] = 255.0;
  landmarks["PC"] = pc;

  WriteITKtoSlicer3Lmk(path, landmarks, SLICER_V3_FCSV);
  const LandmarksMapType loaded = ReadSlicer3toITKLmk(path);

  ASSERT_EQ(loaded.size(), landmarks.size());
  for (const auto & [name, pt] : landmarks)
  {
    auto it = loaded.find(name);
    ASSERT_NE(it, loaded.end()) << "landmark " << name << " not found";
    for (unsigned int d = 0; d < 3; ++d)
    {
      EXPECT_NEAR(it->second[d], pt[d], kCoordTol) << "landmark " << name << " component " << d;
    }
  }

  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// Negative coordinates preserved through double LPS<->RAS conversion
// ---------------------------------------------------------------------------
TEST(Slicer3LandmarkIO, NegativeCoordinatesPreservedThroughDoubleConversion)
{
  const std::string path = UniqueTempPath("neg_coords.fcsv");

  LandmarksMapType  landmarks;
  LandmarkPointType pt;
  pt[0] = -42.5;
  pt[1] = -17.3;
  pt[2] = -0.001;
  landmarks["NEG"] = pt;

  // Write as V4 (LPS -> RAS negation of first two components)
  WriteITKtoSlicer3Lmk(path, landmarks, SLICER_V4_FCSV);
  // Read back (RAS -> LPS negation of first two components)
  const LandmarksMapType loaded = ReadSlicer3toITKLmk(path);

  auto it = loaded.find("NEG");
  ASSERT_NE(it, loaded.end());
  for (unsigned int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(it->second[d], pt[d], kCoordTol) << "component " << d;
  }

  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// Empty file reading throws
// ---------------------------------------------------------------------------
TEST(Slicer3LandmarkIO, EmptyFileThrows)
{
  const std::string path = UniqueTempPath("empty.fcsv");
  {
    std::ofstream ofs(path);
    // Write nothing — empty file
  }

  EXPECT_THROW(ReadSlicer3toITKLmk(path), itk::ExceptionObject);
  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// Comment-only file throws (no recognizable header)
// ---------------------------------------------------------------------------
TEST(Slicer3LandmarkIO, CommentOnlyFileThrows)
{
  const std::string path = UniqueTempPath("comments.fcsv");
  {
    std::ofstream ofs(path);
    ofs << "# This is just a comment\n";
    ofs << "# Another comment\n";
  }

  EXPECT_THROW(ReadSlicer3toITKLmk(path), itk::ExceptionObject);
  std::filesystem::remove(path);
}

// ---------------------------------------------------------------------------
// ReadLandmarkWeights roundtrip
// ---------------------------------------------------------------------------
TEST(Slicer3LandmarkIO, LandmarkWeightsRoundtrip)
{
  const std::string path = UniqueTempPath("weights.csv");
  {
    std::ofstream ofs(path);
    ofs << "# Header comment\n";
    ofs << "Landmark,Weight\n"; // header row (skipped by reader)
    ofs << "AC,0.75\n";
    ofs << "PC,1.25\n";
    ofs << "RP,0.001\n";
  }

  const LandmarksWeightMapType weights = ReadLandmarkWeights(path);

  ASSERT_EQ(weights.size(), 3u);
  EXPECT_NEAR(weights.at("AC"), 0.75, kCoordTol);
  EXPECT_NEAR(weights.at("PC"), 1.25, kCoordTol);
  EXPECT_NEAR(weights.at("RP"), 0.001, kCoordTol);

  std::filesystem::remove(path);
}
