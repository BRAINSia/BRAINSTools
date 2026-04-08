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
/// \file ConvertToRigidAffineGTest.cxx
/// \brief GTest unit tests for ConvertToRigidAffine.h transform conversions (issue #403).
///
/// Tests: identity versor->affine, 90-degree Z rotation, translation-only,
/// affine<->vnl 4x4 roundtrip, extract versor from ScaleVersor3D,
/// extract versor from non-orthogonal affine (verify orthogonality).

#include <gtest/gtest.h>
#include <cmath>

#include "ConvertToRigidAffine.h"

static constexpr double kTol = 1.0e-10;

// -----------------------------------------------------------------------
// Identity versor -> affine
// -----------------------------------------------------------------------
TEST(ConvertToRigidAffine, IdentityVersorToAffine)
{
  auto versor = AssignRigid::VersorRigid3DTransformType::New();
  versor->SetIdentity();

  auto affine = AssignRigid::AffineTransformType::New();
  affine->SetIdentity();

  AssignRigid::AssignConvertedTransform(affine, AssignRigid::VersorRigid3DTransformType::ConstPointer(versor));

  const auto & matrix = affine->GetMatrix();
  for (unsigned int r = 0; r < 3; ++r)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      const double expected = (r == c) ? 1.0 : 0.0;
      EXPECT_NEAR(matrix(r, c), expected, kTol) << "row=" << r << " col=" << c;
    }
  }
  const auto & offset = affine->GetOffset();
  for (unsigned int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(offset[i], 0.0, kTol);
  }
}

// -----------------------------------------------------------------------
// 90-degree Z rotation versor -> affine
// Expected matrix: [0,-1,0; 1,0,0; 0,0,1]
// -----------------------------------------------------------------------
TEST(ConvertToRigidAffine, Rotation90DegZVersorToAffine)
{
  auto versor = AssignRigid::VersorRigid3DTransformType::New();
  versor->SetIdentity();

  // 90 degrees about Z: versor = (0, 0, sin(pi/4)) with w = cos(pi/4)
  AssignRigid::VersorType v;
  using VectorType = itk::Vector<double, 3>;
  VectorType axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  v.Set(axis, M_PI / 2.0);
  versor->SetRotation(v);

  auto affine = AssignRigid::AffineTransformType::New();
  affine->SetIdentity();

  AssignRigid::AssignConvertedTransform(affine, AssignRigid::VersorRigid3DTransformType::ConstPointer(versor));

  const auto & m = affine->GetMatrix();
  // Row 0: [0, -1, 0]
  EXPECT_NEAR(m(0, 0), 0.0, kTol);
  EXPECT_NEAR(m(0, 1), -1.0, kTol);
  EXPECT_NEAR(m(0, 2), 0.0, kTol);
  // Row 1: [1, 0, 0]
  EXPECT_NEAR(m(1, 0), 1.0, kTol);
  EXPECT_NEAR(m(1, 1), 0.0, kTol);
  EXPECT_NEAR(m(1, 2), 0.0, kTol);
  // Row 2: [0, 0, 1]
  EXPECT_NEAR(m(2, 0), 0.0, kTol);
  EXPECT_NEAR(m(2, 1), 0.0, kTol);
  EXPECT_NEAR(m(2, 2), 1.0, kTol);
}

// -----------------------------------------------------------------------
// Translation-only versor -> affine preserves translation
// -----------------------------------------------------------------------
TEST(ConvertToRigidAffine, TranslationOnlyVersorToAffine)
{
  auto versor = AssignRigid::VersorRigid3DTransformType::New();
  versor->SetIdentity();

  // Non-trivial translation values
  AssignRigid::VersorRigid3DTransformType::OutputVectorType translation;
  translation[0] = 42.5;
  translation[1] = -17.3;
  translation[2] = 0.001;
  versor->SetTranslation(translation);

  auto affine = AssignRigid::AffineTransformType::New();
  affine->SetIdentity();

  AssignRigid::AssignConvertedTransform(affine, AssignRigid::VersorRigid3DTransformType::ConstPointer(versor));

  // Matrix should be identity
  const auto & m = affine->GetMatrix();
  for (unsigned int r = 0; r < 3; ++r)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      EXPECT_NEAR(m(r, c), (r == c) ? 1.0 : 0.0, kTol);
    }
  }

  // Translation should be preserved
  const auto & t = affine->GetTranslation();
  EXPECT_NEAR(t[0], 42.5, kTol);
  EXPECT_NEAR(t[1], -17.3, kTol);
  EXPECT_NEAR(t[2], 0.001, kTol);
}

// -----------------------------------------------------------------------
// Affine -> VnlTransformMatrixType44 -> Affine roundtrip
// -----------------------------------------------------------------------
TEST(ConvertToRigidAffine, AffineToVnl44Roundtrip)
{
  auto affineOrig = AssignRigid::AffineTransformType::New();
  affineOrig->SetIdentity();

  // Set a non-trivial rotation via versor, then convert
  auto versor = AssignRigid::VersorRigid3DTransformType::New();
  versor->SetIdentity();
  AssignRigid::VersorType v;
  using VectorType = itk::Vector<double, 3>;
  VectorType axis;
  axis[0] = 1.0;
  axis[1] = 0.0;
  axis[2] = 0.0;
  v.Set(axis, M_PI / 3.0); // 60 degrees about X
  versor->SetRotation(v);

  AssignRigid::VersorRigid3DTransformType::OutputVectorType translation;
  translation[0] = 42.5;
  translation[1] = -17.3;
  translation[2] = 255.0;
  versor->SetTranslation(translation);

  AssignRigid::AssignConvertedTransform(affineOrig, AssignRigid::VersorRigid3DTransformType::ConstPointer(versor));

  // Affine -> Vnl 4x4
  AssignRigid::VnlTransformMatrixType44 vnlMat;
  vnlMat.set_identity();
  AssignRigid::AssignConvertedTransform(vnlMat, AssignRigid::AffineTransformType::ConstPointer(affineOrig));

  // Vnl 4x4 -> Affine
  auto affineBack = AssignRigid::AffineTransformType::New();
  affineBack->SetIdentity();
  AssignRigid::AssignConvertedTransform(affineBack, vnlMat);

  // Compare matrices
  const auto & mOrig = affineOrig->GetMatrix();
  const auto & mBack = affineBack->GetMatrix();
  for (unsigned int r = 0; r < 3; ++r)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      EXPECT_NEAR(mBack(r, c), mOrig(r, c), kTol) << "row=" << r << " col=" << c;
    }
  }

  // Compare offsets
  const auto & oOrig = affineOrig->GetOffset();
  const auto & oBack = affineBack->GetOffset();
  for (unsigned int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(oBack[i], oOrig[i], kTol) << "offset[" << i << "]";
  }
}

// -----------------------------------------------------------------------
// Extract versor from ScaleVersor3DTransform
// -----------------------------------------------------------------------
TEST(ConvertToRigidAffine, ExtractVersorFromScaleVersor)
{
  auto scaleVersor = AssignRigid::ScaleVersor3DTransformType::New();
  scaleVersor->SetIdentity();

  // Set a rotation: 45 degrees about Y
  AssignRigid::VersorType v;
  using VectorType = itk::Vector<double, 3>;
  VectorType axis;
  axis[0] = 0.0;
  axis[1] = 1.0;
  axis[2] = 0.0;
  v.Set(axis, M_PI / 4.0);
  scaleVersor->SetRotation(v);

  // Set non-uniform scale
  AssignRigid::ScaleVersor3DTransformType::ScaleVectorType scale;
  scale[0] = 1.5;
  scale[1] = 2.0;
  scale[2] = 0.75;
  scaleVersor->SetScale(scale);

  // Set non-trivial translation
  AssignRigid::ScaleVersor3DTransformType::OutputVectorType translation;
  translation[0] = -17.3;
  translation[1] = 42.5;
  translation[2] = 0.001;
  scaleVersor->SetTranslation(translation);

  auto rigid = AssignRigid::VersorRigid3DTransformType::New();
  rigid->SetIdentity();

  AssignRigid::ExtractVersorRigid3DTransform(rigid, AssignRigid::ScaleVersor3DTransformType::ConstPointer(scaleVersor));

  // The extracted versor should match the original versor
  const auto & vExtracted = rigid->GetVersor();
  EXPECT_NEAR(vExtracted.GetX(), v.GetX(), kTol);
  EXPECT_NEAR(vExtracted.GetY(), v.GetY(), kTol);
  EXPECT_NEAR(vExtracted.GetZ(), v.GetZ(), kTol);
  EXPECT_NEAR(vExtracted.GetW(), v.GetW(), kTol);

  // Translation should be preserved
  const auto & t = rigid->GetTranslation();
  EXPECT_NEAR(t[0], -17.3, kTol);
  EXPECT_NEAR(t[1], 42.5, kTol);
  EXPECT_NEAR(t[2], 0.001, kTol);
}

// -----------------------------------------------------------------------
// Extract versor from non-orthogonal affine (verify result is orthogonal)
// -----------------------------------------------------------------------
TEST(ConvertToRigidAffine, ExtractVersorFromNonOrthogonalAffine)
{
  auto affine = AssignRigid::AffineTransformType::New();
  affine->SetIdentity();

  // Create a non-orthogonal matrix (scale + shear)
  AssignRigid::MatrixType nonOrthog;
  nonOrthog(0, 0) = 1.5;
  nonOrthog(0, 1) = 0.3;
  nonOrthog(0, 2) = 0.0;
  nonOrthog(1, 0) = 0.0;
  nonOrthog(1, 1) = 2.0;
  nonOrthog(1, 2) = -0.1;
  nonOrthog(2, 0) = 0.0;
  nonOrthog(2, 1) = 0.0;
  nonOrthog(2, 2) = 0.75;
  affine->SetMatrix(nonOrthog);

  AssignRigid::AffineTransformType::OutputVectorType translation;
  translation[0] = 42.5;
  translation[1] = -17.3;
  translation[2] = 255.0;
  affine->SetTranslation(translation);

  auto rigid = AssignRigid::VersorRigid3DTransformType::New();
  rigid->SetIdentity();

  AssignRigid::ExtractVersorRigid3DTransform(rigid, AssignRigid::AffineTransformType::ConstPointer(affine));

  // Verify the resulting rotation matrix is orthogonal: R^T * R = I
  const auto & m = rigid->GetMatrix();
  for (unsigned int r = 0; r < 3; ++r)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      double dot = 0.0;
      for (unsigned int k = 0; k < 3; ++k)
      {
        dot += m(k, r) * m(k, c); // column r dot column c
      }
      const double expected = (r == c) ? 1.0 : 0.0;
      EXPECT_NEAR(dot, expected, 1.0e-6) << "R^T*R(" << r << "," << c << ")";
    }
  }

  // Verify determinant is +1 (proper rotation, not reflection)
  const double det = m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
                     m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
                     m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
  EXPECT_NEAR(det, 1.0, 1.0e-6);
}
