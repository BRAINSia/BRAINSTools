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
/// \file DWIMetaDataDictionaryValidatorGTest.cxx
/// \brief GTest coverage for DWIMetaDataDictionaryValidator (issue #403).
///
/// All tests are in-memory metadata operations; no file IO is needed.

#include <gtest/gtest.h>
#include <cmath>
#include "DWIMetaDataDictionaryValidator.h"

static constexpr double kTol = 1e-4;

// ---------------------------------------------------------------------------
// Default constructor creates an empty dictionary
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, DefaultConstructorCreatesEmptyDict)
{
  DWIMetaDataDictionaryValidator validator;
  EXPECT_EQ(validator.GetGradientCount(), 0);
  // GetMetaDataDictionary() should return a reference to the internal dict;
  // the dict should have no keys.
  const auto & dict = validator.GetMetaDataDictionary();
  EXPECT_TRUE(dict.GetKeys().empty());
}

// ---------------------------------------------------------------------------
// Modality roundtrip
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, ModalityRoundtrip)
{
  DWIMetaDataDictionaryValidator validator;
  validator.SetModality("DWMRI");
  EXPECT_EQ(validator.GetModality(), "DWMRI");
}

// ---------------------------------------------------------------------------
// BValue roundtrip with non-trivial value
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, BValueRoundtripNonTrivial)
{
  DWIMetaDataDictionaryValidator validator;
  const double                   bval = 3333.333;
  validator.SetBValue(bval);
  EXPECT_NEAR(validator.GetBValue(), bval, kTol);
}

// ---------------------------------------------------------------------------
// BValue zero
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, BValueZero)
{
  DWIMetaDataDictionaryValidator validator;
  validator.SetBValue(0.0);
  EXPECT_NEAR(validator.GetBValue(), 0.0, kTol);
}

// ---------------------------------------------------------------------------
// Measurement frame roundtrip (30-degree rotation about Z)
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, MeasurementFrameRoundtripRotation)
{
  DWIMetaDataDictionaryValidator                     validator;
  DWIMetaDataDictionaryValidator::RotationMatrixType mf;

  const double angle = 30.0 * M_PI / 180.0;
  mf(0, 0) = std::cos(angle);
  mf(0, 1) = -std::sin(angle);
  mf(0, 2) = 0.0;
  mf(1, 0) = std::sin(angle);
  mf(1, 1) = std::cos(angle);
  mf(1, 2) = 0.0;
  mf(2, 0) = 0.0;
  mf(2, 1) = 0.0;
  mf(2, 2) = 1.0;

  validator.SetMeasurementFrame(mf);
  const auto result = validator.GetMeasurementFrame();
  for (unsigned int r = 0; r < 3; ++r)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      EXPECT_NEAR(result(r, c), mf(r, c), kTol) << "mismatch at (" << r << "," << c << ")";
    }
  }
}

// ---------------------------------------------------------------------------
// Measurement frame identity
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, MeasurementFrameIdentity)
{
  DWIMetaDataDictionaryValidator                     validator;
  DWIMetaDataDictionaryValidator::RotationMatrixType mf;
  mf.SetIdentity();

  validator.SetMeasurementFrame(mf);
  const auto result = validator.GetMeasurementFrame();
  for (unsigned int r = 0; r < 3; ++r)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      EXPECT_NEAR(result(r, c), mf(r, c), kTol);
    }
  }
}

// ---------------------------------------------------------------------------
// Gradient table roundtrip (8 gradients)
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, GradientTableRoundtrip)
{
  DWIMetaDataDictionaryValidator                    validator;
  DWIMetaDataDictionaryValidator::GradientTableType table;

  // Zero vector
  DWIMetaDataDictionaryValidator::GradientDirectionType g0;
  g0[0] = 0.0;
  g0[1] = 0.0;
  g0[2] = 0.0;
  table.push_back(g0);

  // Unit vectors along axes
  DWIMetaDataDictionaryValidator::GradientDirectionType g1;
  g1[0] = 1.0;
  g1[1] = 0.0;
  g1[2] = 0.0;
  table.push_back(g1);

  DWIMetaDataDictionaryValidator::GradientDirectionType g2;
  g2[0] = 0.0;
  g2[1] = 1.0;
  g2[2] = 0.0;
  table.push_back(g2);

  DWIMetaDataDictionaryValidator::GradientDirectionType g3;
  g3[0] = 0.0;
  g3[1] = 0.0;
  g3[2] = 1.0;
  table.push_back(g3);

  // Non-unit vector
  DWIMetaDataDictionaryValidator::GradientDirectionType g4;
  g4[0] = 0.577;
  g4[1] = 0.577;
  g4[2] = 0.577;
  table.push_back(g4);

  // Negative components
  DWIMetaDataDictionaryValidator::GradientDirectionType g5;
  g5[0] = -0.707;
  g5[1] = 0.707;
  g5[2] = 0.0;
  table.push_back(g5);

  DWIMetaDataDictionaryValidator::GradientDirectionType g6;
  g6[0] = 0.3;
  g6[1] = -0.8;
  g6[2] = 0.5;
  table.push_back(g6);

  DWIMetaDataDictionaryValidator::GradientDirectionType g7;
  g7[0] = -0.123;
  g7[1] = -0.456;
  g7[2] = 0.789;
  table.push_back(g7);

  validator.SetGradientTable(table);
  EXPECT_EQ(validator.GetGradientCount(), 8);

  const auto result = validator.GetGradientTable();
  ASSERT_EQ(result.size(), table.size());
  for (size_t i = 0; i < table.size(); ++i)
  {
    for (unsigned int d = 0; d < 3; ++d)
    {
      EXPECT_NEAR(result[i][d], table[i][d], kTol) << "gradient " << i << " component " << d;
    }
  }
}

// ---------------------------------------------------------------------------
// Gradient table overwrite (set 8 then set 4, verify count=4)
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, GradientTableOverwrite)
{
  DWIMetaDataDictionaryValidator                        validator;
  DWIMetaDataDictionaryValidator::GradientTableType     table8;
  DWIMetaDataDictionaryValidator::GradientDirectionType g;
  g[0] = 1.0;
  g[1] = 0.0;
  g[2] = 0.0;
  for (int i = 0; i < 8; ++i)
  {
    table8.push_back(g);
  }
  validator.SetGradientTable(table8);
  EXPECT_EQ(validator.GetGradientCount(), 8);

  // Now set only 4
  DWIMetaDataDictionaryValidator::GradientTableType table4;
  g[0] = 0.0;
  g[1] = 1.0;
  g[2] = 0.0;
  for (int i = 0; i < 4; ++i)
  {
    table4.push_back(g);
  }
  validator.SetGradientTable(table4);
  EXPECT_EQ(validator.GetGradientCount(), 4);

  const auto result = validator.GetGradientTable();
  ASSERT_EQ(result.size(), 4u);
}

// ---------------------------------------------------------------------------
// Delete gradient table
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, DeleteGradientTable)
{
  DWIMetaDataDictionaryValidator                        validator;
  DWIMetaDataDictionaryValidator::GradientTableType     table;
  DWIMetaDataDictionaryValidator::GradientDirectionType g;
  g[0] = 1.0;
  g[1] = 0.0;
  g[2] = 0.0;
  table.push_back(g);
  table.push_back(g);
  table.push_back(g);
  validator.SetGradientTable(table);
  EXPECT_EQ(validator.GetGradientCount(), 3);

  validator.DeleteGradientTable();
  EXPECT_EQ(validator.GetGradientCount(), 0);
}

// ---------------------------------------------------------------------------
// Single gradient set/get with {42.5, -17.3, 0.001}
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, SingleGradientSetGet)
{
  DWIMetaDataDictionaryValidator                     validator;
  DWIMetaDataDictionaryValidator::Double3x1ArrayType grad;
  grad[0] = 42.5;
  grad[1] = -17.3;
  grad[2] = 0.001;

  validator.SetGradient(0, grad);
  const auto result = validator.GetGradient(0);
  EXPECT_NEAR(result[0], 42.5, kTol);
  EXPECT_NEAR(result[1], -17.3, kTol);
  EXPECT_NEAR(result[2], 0.001, kTol);
}

// ---------------------------------------------------------------------------
// Centerings roundtrip
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, CenteringsRoundtrip)
{
  DWIMetaDataDictionaryValidator validator;
  std::vector<std::string>       centerings = { "cell", "cell", "cell", "???" };
  validator.SetCenterings(centerings);

  const auto result = validator.GetCenterings();
  ASSERT_EQ(result.size(), 4u);
  EXPECT_EQ(result[0], "cell");
  EXPECT_EQ(result[1], "cell");
  EXPECT_EQ(result[2], "cell");
  // The 4th element "???" is skipped by SetCenterings (only non-"???" are stored),
  // and GetCenterings returns "???" as the default for missing keys.
  EXPECT_EQ(result[3], "???");
}

// ---------------------------------------------------------------------------
// SetMetaDataDictionary copies all fields
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, SetMetaDataDictionaryCopiesAllFields)
{
  DWIMetaDataDictionaryValidator src;
  src.SetModality("DWMRI");
  src.SetBValue(1500.5);

  DWIMetaDataDictionaryValidator::RotationMatrixType mf;
  mf.SetIdentity();
  src.SetMeasurementFrame(mf);

  DWIMetaDataDictionaryValidator::Double3x1ArrayType grad;
  grad[0] = 0.5;
  grad[1] = 0.5;
  grad[2] = 0.707;
  src.SetGradient(0, grad);

  // Copy via SetMetaDataDictionary
  DWIMetaDataDictionaryValidator dst;
  dst.SetMetaDataDictionary(src.GetMetaDataDictionary());

  EXPECT_EQ(dst.GetModality(), "DWMRI");
  EXPECT_NEAR(dst.GetBValue(), 1500.5, kTol);
  EXPECT_EQ(dst.GetGradientCount(), 1);

  const auto dstGrad = dst.GetGradient(0);
  EXPECT_NEAR(dstGrad[0], 0.5, kTol);
  EXPECT_NEAR(dstGrad[1], 0.5, kTol);
  EXPECT_NEAR(dstGrad[2], 0.707, kTol);

  const auto dstMF = dst.GetMeasurementFrame();
  for (unsigned int r = 0; r < 3; ++r)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      EXPECT_NEAR(dstMF(r, c), mf(r, c), kTol);
    }
  }
}

// ---------------------------------------------------------------------------
// Gradient with negative values preserves sign
// ---------------------------------------------------------------------------
TEST(DWIMetaDataDictionaryValidator, GradientNegativeValuesPreserveSign)
{
  DWIMetaDataDictionaryValidator                     validator;
  DWIMetaDataDictionaryValidator::Double3x1ArrayType grad;
  grad[0] = -42.5;
  grad[1] = -17.3;
  grad[2] = -0.001;

  validator.SetGradient(0, grad);
  const auto result = validator.GetGradient(0);
  EXPECT_NEAR(result[0], -42.5, kTol);
  EXPECT_NEAR(result[1], -17.3, kTol);
  EXPECT_NEAR(result[2], -0.001, kTol);
}
