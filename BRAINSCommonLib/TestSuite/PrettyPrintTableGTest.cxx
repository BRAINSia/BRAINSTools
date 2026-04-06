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
/// \file PrettyPrintTableGTest.cxx
/// \brief GTest port of PrettyPrintTableTest.cxx (issue #403 — GTest migration).
///
/// Verifies that PrettyPrintTable produces non-empty output and that the
/// BRAINSTools::Version API returns non-empty strings.  The original test
/// had no assertions beyond EXIT_SUCCESS; the GTest version makes the
/// smoke-test intent explicit.

#include <gtest/gtest.h>
#include <sstream>

#include "PrettyPrintTable.h"
#include "BRAINSToolsVersion.h"

// ---------------------------------------------------------------------------
// PrettyPrintTable smoke tests
// ---------------------------------------------------------------------------

TEST(PrettyPrintTable, PrintsNonEmptyOutput)
{
  PrettyPrintTable p;
  p.setTablePad(5);
  p.add(0, 0, "String");
  p.add(0, 1, "SecondColumn");
  p.add(0, 2, "4C");
  p.add(0, 3, "5C");
  p.add(1, 0, "Integers");
  p.add(1, 1, 1, "%d");
  p.add(1, 2, 2, "%d");
  p.add(1, 3, 3, "%d");
  p.add(1, 0, "ZeroPadInt");
  p.add(1, 1, 1, "%02d");
  p.add(1, 2, 2, "%02d");
  p.add(1, 3, 3, "%02d");
  p.add(2, 0, "FloatingPoint");
  p.add(2, 1, 1.0F, "%+5.2f");
  p.add(2, 2, 2.0F, "%+5.2f");
  p.add(2, 3, 3.0F, "%+5.2f");

  std::ostringstream oss;
  EXPECT_NO_THROW(p.Print(oss));
  EXPECT_FALSE(oss.str().empty());
}

TEST(PrettyPrintTable, RightJustifyProducesNonEmptyOutput)
{
  PrettyPrintTable p;
  p.setTablePad(5);
  p.add(0, 0, "String");
  p.add(0, 1, "SecondColumn");
  p.add(1, 0, "Integers");
  p.add(1, 1, 1, "%d");

  p.rightJustify();

  std::ostringstream oss;
  EXPECT_NO_THROW(p.Print(oss));
  EXPECT_FALSE(oss.str().empty());
}

// ---------------------------------------------------------------------------
// BRAINSTools::Version API smoke tests
// ---------------------------------------------------------------------------

TEST(BRAINSToolsVersion, VersionStringIsNonEmpty) { EXPECT_FALSE(BRAINSTools::Version::VersionString().empty()); }

TEST(BRAINSToolsVersion, ExtendedVersionStringIsNonEmpty)
{
  EXPECT_FALSE(BRAINSTools::Version::ExtendedVersionString().empty());
}

TEST(BRAINSToolsVersion, ToStringIsNonEmpty) { EXPECT_FALSE(BRAINSTools::Version::ToString().empty()); }

TEST(BRAINSToolsVersion, BuildDateIsNonEmpty) { EXPECT_FALSE(BRAINSTools::Version::BuildDate().empty()); }
