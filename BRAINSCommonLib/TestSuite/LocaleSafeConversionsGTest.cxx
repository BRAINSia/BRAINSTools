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
/// \file LocaleSafeConversionsGTest.cxx
/// \brief GTest-based unit tests for BRAINSCommonLib/LocaleSafeConversions.h
///
/// This file is the GTest port of LocaleSafeConversionsTest.cxx.
/// Every test category from the original custom-harness file is preserved;
/// the locale-independence suite additionally uses GTest parameterized tests
/// (TEST_P) so that each locale becomes a named, individually-skippable
/// CTest entry rather than a single monolithic test.
///
/// Test suites:
///   LocaleSafeConversions/RoundTrip*      — basic correctness
///   LocaleSafeConversions/DicomPadding*   — Siemens/Philips header patterns
///   LocaleSafeConversions/SentryBug*      — std::ws + eofbit regression
///   LocaleSafeConversions/ErrorHandling*  — invalid inputs throw
///   LocaleSafeConversions/StodMatch*      — behaviours matching std::stod
///   LocaleSafeConversions/StoiMatch*      — behaviours matching std::stoi
///   LocaleSafeConversions/StodDiff*       — intentional differences from std::stod/atof
///   LocaleSafeConversions/StoiDiff*       — intentional differences from std::stoi/atoi
///   LocaleSafeConversions/SafeStoui*      — safe_stoui-specific behaviour
///   LocaleIndependenceTest/<locale>       — one named test per locale candidate

#include "LocaleSafeConversions.h"

#include <climits>
#include <locale>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

// ---------------------------------------------------------------------------
// Helper: construct a std::string containing embedded null bytes.
// sizeof(arr)-1 strips the implicit C-string terminator.
// ---------------------------------------------------------------------------
#define STR_WITH_NULLS(arr) std::string((arr), sizeof(arr) - 1)

// ===========================================================================
// 1. Round-trip correctness
// ===========================================================================

TEST(LocaleSafeConversions, RoundTripDouble)
{
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14159"), 3.14159, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("-2.718"), -2.718, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("0.0"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("1e10"), 1e10, 1.0);
  EXPECT_NEAR(BRAINSTools::safe_stod("-1.5e-3"), -1.5e-3, 1e-10);
  // BCD real-world value from DWI metadata
  EXPECT_NEAR(BRAINSTools::safe_stod("1.2303285706216762"), 1.2303285706216762, 1e-12);
}

TEST(LocaleSafeConversions, RoundTripFloat)
{
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5"), 1.5f, 1e-5f);
  EXPECT_NEAR(BRAINSTools::safe_stof("-0.25"), -0.25f, 1e-5f);
}

TEST(LocaleSafeConversions, RoundTripInt)
{
  EXPECT_EQ(BRAINSTools::safe_stoi("42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("-7"), -7);
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);
}

TEST(LocaleSafeConversions, RoundTripUInt)
{
  EXPECT_EQ(BRAINSTools::safe_stoui("100"), 100u);
  EXPECT_EQ(BRAINSTools::safe_stoui("0"), 0u);
}

// ===========================================================================
// 2. Locale independence — parameterized, one test per locale candidate
//
// Each locale in the INSTANTIATE list gets its own named CTest entry:
//   LocaleIndependenceTest/CommaDecimalLocales/de_DE.UTF-8
//   LocaleIndependenceTest/CommaDecimalLocales/fr_FR.UTF-8
//   ...
// Locales not installed on the host system are individually SKIPPED, leaving
// a clear record of which were exercised and which were absent.
// ===========================================================================

class LocaleIndependenceTest : public ::testing::TestWithParam<std::string>
{
protected:
  // Save and restore the global locale around each test case so that one
  // locale failure cannot contaminate subsequent test cases.
  std::locale original_;

  void
  SetUp() override
  {
    original_ = std::locale::global(std::locale::classic());
  }

  void
  TearDown() override
  {
    std::locale::global(original_);
  }
};

TEST_P(LocaleIndependenceTest, ParseDotDecimalUnderCommaLocale)
{
  const std::string name = GetParam();
  try
  {
    std::locale::global(std::locale(name.c_str()));
  }
  catch (const std::runtime_error &)
  {
    GTEST_SKIP() << "locale " << name << " not installed on this system";
  }

  // Under any comma-decimal locale, safe_ must still parse dot-decimal.
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5"), 1.5f, 1e-5f);
  EXPECT_EQ(BRAINSTools::safe_stoi("99"), 99);
  EXPECT_EQ(BRAINSTools::safe_stoui("7"), 7u);

  // Comma-decimal strings must still be rejected.
  EXPECT_THROW(BRAINSTools::safe_stod("3,14"), std::invalid_argument);
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(
  CommaDecimalLocales,
  LocaleIndependenceTest,
  ::testing::Values(
    // Germanic
    std::string("de_DE.UTF-8"), std::string("de_DE"),
    std::string("de_AT.UTF-8"), std::string("de_AT"),
    std::string("de_CH.UTF-8"), std::string("de_CH"),
    std::string("nl_NL.UTF-8"), std::string("nl_NL"),
    std::string("nl_BE.UTF-8"), std::string("nl_BE"),
    std::string("af_ZA.UTF-8"), std::string("af_ZA"),
    // Romance
    std::string("fr_FR.UTF-8"), std::string("fr_FR"),
    std::string("fr_BE.UTF-8"), std::string("fr_BE"),
    std::string("fr_CA.UTF-8"), std::string("fr_CA"),
    std::string("fr_CH.UTF-8"), std::string("fr_CH"),
    std::string("es_ES.UTF-8"), std::string("es_ES"),
    std::string("es_AR.UTF-8"), std::string("es_AR"),
    std::string("it_IT.UTF-8"), std::string("it_IT"),
    std::string("it_CH.UTF-8"), std::string("it_CH"),
    std::string("pt_PT.UTF-8"), std::string("pt_PT"),
    std::string("pt_BR.UTF-8"), std::string("pt_BR"),
    std::string("ca_ES.UTF-8"), std::string("ca_ES"),
    std::string("ro_RO.UTF-8"), std::string("ro_RO"),
    // Nordic
    std::string("da_DK.UTF-8"), std::string("da_DK"),
    std::string("sv_SE.UTF-8"), std::string("sv_SE"),
    std::string("fi_FI.UTF-8"), std::string("fi_FI"),
    std::string("nb_NO.UTF-8"), std::string("nb_NO"),
    // Slavic
    std::string("pl_PL.UTF-8"),
    std::string("cs_CZ.UTF-8"), std::string("cs_CZ"),
    std::string("sk_SK.UTF-8"), std::string("sk_SK"),
    std::string("hr_HR.UTF-8"), std::string("hr_HR"),
    std::string("bg_BG.UTF-8"), std::string("bg_BG"),
    std::string("ru_RU.UTF-8"),
    std::string("uk_UA.UTF-8"), std::string("uk_UA"),
    // Baltic
    std::string("lv_LV.UTF-8"),
    std::string("lt_LT.UTF-8"),
    std::string("et_EE.UTF-8"),
    // Other
    std::string("hu_HU.UTF-8"), std::string("hu_HU"),
    std::string("el_GR.UTF-8"),
    std::string("tr_TR.UTF-8")
  )
);

// Dot-decimal locales: English, East Asian, SE Asian, Middle Eastern, South Asian.
// Under these locales both std::stod and safe_stod naturally parse dot-decimal,
// but it is still important to verify no unexpected interaction with the
// "C"-locale-imbued istringstream inside safe_.
// clang-format off
INSTANTIATE_TEST_SUITE_P(
  DotDecimalLocales,
  LocaleIndependenceTest,
  ::testing::Values(
    // English
    std::string("en_US.UTF-8"), std::string("en_US"),
    std::string("en_GB.UTF-8"), std::string("en_GB"),
    std::string("en_AU.UTF-8"),
    std::string("en_CA.UTF-8"),
    std::string("en_IE.UTF-8"),
    std::string("en_NZ.UTF-8"),
    std::string("en_ZA.UTF-8"),
    // East Asian
    std::string("ja_JP.UTF-8"),
    std::string("zh_CN.UTF-8"),
    std::string("zh_TW.UTF-8"),
    std::string("ko_KR.UTF-8"),
    // Southeast Asian
    std::string("th_TH.UTF-8"),
    std::string("ms_MY.UTF-8"),
    std::string("id_ID.UTF-8"),
    std::string("vi_VN.UTF-8"),
    // Middle Eastern
    std::string("ar_SA.UTF-8"),
    std::string("ar_EG.UTF-8"),
    std::string("he_IL.UTF-8"),
    // South Asian
    std::string("hi_IN.UTF-8")
  )
);
// clang-format on

// ===========================================================================
// 3. DICOM CSA header padding patterns
//    Siemens: ASCII padded with '\n' then '\0' bytes to 4-byte boundary
//    Philips: integer fields stored as "0000\0\0\0\0"
// ===========================================================================

TEST(LocaleSafeConversions, DicomPaddingTrailingNewline)
{
  // "0.00000000\n" — Siemens CSA float value
  EXPECT_NEAR(BRAINSTools::safe_stod("0.00000000\n"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14 "), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5\n"), 1.5f, 1e-5f);
  EXPECT_EQ(BRAINSTools::safe_stoi("42\n"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoui("7\t"), 7u);
}

TEST(LocaleSafeConversions, DicomPaddingNewlineFollowedByNulls)
{
  // "0.00000000\n\0\0" — Siemens CSA, 4-byte aligned; null bytes must be
  // stripped by detail::null_truncate() before parsing.
  EXPECT_NEAR(BRAINSTools::safe_stod(STR_WITH_NULLS("0.00000000\n\0\0")), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod(STR_WITH_NULLS("0.00000000\n\0\0\0")), 0.0, 1e-15);
}

TEST(LocaleSafeConversions, DicomPaddingNullPaddedIntegers)
{
  // Philips header: "0000\0\0\0\0"
  EXPECT_EQ(BRAINSTools::safe_stoi(STR_WITH_NULLS("0000\0\0\0\0")), 0);
  EXPECT_EQ(BRAINSTools::safe_stoui(STR_WITH_NULLS("0000\0")), 0u);
  // SpaceThicknessDiff field: "2\0\0\0"
  EXPECT_EQ(BRAINSTools::safe_stoi(STR_WITH_NULLS("2\0\0\0")), 2);
}

TEST(LocaleSafeConversions, DicomPaddingExactSingleChars)
{
  // Strings exactly consumed to EOF — sentry-bug regression (see suite 4).
  EXPECT_NEAR(BRAINSTools::safe_stod("0"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.2303285706216762"), 1.2303285706216762, 1e-12);
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("0000"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("2"), 2);
  EXPECT_EQ(BRAINSTools::safe_stoui("0"), 0u);
}

// ===========================================================================
// 4. Sentry-bug regression
//    Calling std::ws when eofbit is set causes the sentry to call
//    setstate(failbit), making valid inputs like "3.14" throw.
//    detail::only_whitespace_remains() guards with !iss.eof().
// ===========================================================================

TEST(LocaleSafeConversions, SentryBugRegression)
{
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14159"), 3.14159, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("-2.718"), -2.718, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("0"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5"), 1.5f, 1e-5f);
  EXPECT_EQ(BRAINSTools::safe_stoi("42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("0000"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoui("100"), 100u);
}

// ===========================================================================
// 5. Error handling — invalid inputs must throw std::invalid_argument
// ===========================================================================

TEST(LocaleSafeConversions, ErrorEmptyString)
{
  EXPECT_THROW(BRAINSTools::safe_stod(""), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stof(""), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi(""), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoui(""), std::invalid_argument);
}

TEST(LocaleSafeConversions, ErrorNullOnlyString)
{
  EXPECT_THROW(BRAINSTools::safe_stod(STR_WITH_NULLS("\0\0\0")), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi(STR_WITH_NULLS("\0")), std::invalid_argument);
}

TEST(LocaleSafeConversions, ErrorWhitespaceOnly)
{
  EXPECT_THROW(BRAINSTools::safe_stod("   "), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("   "), std::invalid_argument);
}

TEST(LocaleSafeConversions, ErrorTrailingNonWhitespace)
{
  EXPECT_THROW(BRAINSTools::safe_stod("3.14abc"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stof("1.0f"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("42xyz"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoui("5u"), std::invalid_argument);
  // Second decimal point counts as trailing non-whitespace
  EXPECT_THROW(BRAINSTools::safe_stod("3.14.15"), std::invalid_argument);
  // Embedded space between digits
  EXPECT_THROW(BRAINSTools::safe_stoi("4 2"), std::invalid_argument);
}

TEST(LocaleSafeConversions, ErrorNonNumeric)
{
  EXPECT_THROW(BRAINSTools::safe_stod("nan_not_accepted"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("abc"), std::invalid_argument);
}

// ===========================================================================
// 6. Behaviours matching std::stod (C++11)
//
//    safe_stod preserves all std::stod behaviours except partial parse.
//    See suite 8 for intentional differences.
// ===========================================================================

TEST(LocaleSafeConversions, StodMatchLeadingWhitespace)
{
  EXPECT_NEAR(BRAINSTools::safe_stod("  3.14"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("\t-1.0"), -1.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("\n2.5"), 2.5, 1e-5);
}

TEST(LocaleSafeConversions, StodMatchPlusSign)
{
  EXPECT_NEAR(BRAINSTools::safe_stod("+3.14"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stof("+1.5"), 1.5f, 1e-5f);
}

TEST(LocaleSafeConversions, StodMatchScientificNotation)
{
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5e3"), 1500.0, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5E3"), 1500.0, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5e+3"), 1500.0, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5e-3"), 0.0015, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1e0"), 1.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stof("2.5E2"), 250.0f, 1e-3f);
}

TEST(LocaleSafeConversions, StodMatchLeadingZeros)
{
  // Leading zeros in mantissa — decimal digits, not octal
  EXPECT_NEAR(BRAINSTools::safe_stod("007.5"), 7.5, 1e-10);
  EXPECT_NEAR(BRAINSTools::safe_stod("0042.5"), 42.5, 1e-10);
  EXPECT_NEAR(BRAINSTools::safe_stod("00.5"), 0.5, 1e-10);
}

TEST(LocaleSafeConversions, StodMatchTrailingWhitespace)
{
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14 "), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14\t"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14\n"), 3.14, 1e-5);
}

// ===========================================================================
// 7. Behaviours matching std::stoi (C++11)
//
//    safe_stoi preserves all std::stoi behaviours except partial parse.
//    See suite 9 for intentional differences.
//
//    Notable: leading zeros are DECIMAL digits (base 10), not octal.
//    "042" -> 42 (decimal), NOT 34 (octal).
// ===========================================================================

TEST(LocaleSafeConversions, StoiMatchLeadingWhitespace)
{
  EXPECT_EQ(BRAINSTools::safe_stoi("  42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("\t-7"), -7);
  EXPECT_EQ(BRAINSTools::safe_stoi("\n100"), 100);
}

TEST(LocaleSafeConversions, StoiMatchPlusSign)
{
  EXPECT_EQ(BRAINSTools::safe_stoi("+42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("+0"), 0);
}

TEST(LocaleSafeConversions, StoiMatchLeadingZerosAreDecimal)
{
  // Base 10: leading zeros are NOT octal.
  // std::stoi("042") == 42 (decimal), NOT 34 (octal 042).
  EXPECT_EQ(BRAINSTools::safe_stoi("007"), 7);
  EXPECT_EQ(BRAINSTools::safe_stoi("042"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("0042"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("00"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("000"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("0000"), 0);
}

TEST(LocaleSafeConversions, StoiMatchBoundaryValues)
{
  EXPECT_EQ(BRAINSTools::safe_stoi("2147483647"), INT_MAX);
  EXPECT_EQ(BRAINSTools::safe_stoi("-2147483648"), INT_MIN);
}

TEST(LocaleSafeConversions, StoiMatchTrailingWhitespace)
{
  EXPECT_EQ(BRAINSTools::safe_stoi("42 "), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("42\n"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("42\t"), 42);
}

// ===========================================================================
// 8. Intentional differences from std::stod / atof
//
//    safe_stod is STRICTER in these ways:
//      a) Partial parse rejected:  std::stod("3.14abc") -> 3.14;
//                                  safe_stod("3.14abc") -> throws
//      b) Empty / whitespace:      atof("") -> 0.0;
//                                  safe_stod("") -> throws
//      c) inf / nan:               std::stod("inf") -> +inf;
//                                  safe_stod("inf") -> throws
//                                  (istringstream on Apple Clang/libc++ does
//                                  not parse "inf"/"nan" via locale::classic())
//      d) Overflow exception type: std::stod("1e9999") -> std::out_of_range;
//                                  safe_stod("1e9999") -> std::invalid_argument
//                                  (istringstream sets failbit on overflow)
// ===========================================================================

TEST(LocaleSafeConversions, StodDiffPartialParseRejected)
{
  // std::stod("3.14abc") returns 3.14; safe_stod must throw.
  EXPECT_THROW(BRAINSTools::safe_stod("3.14abc"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("3.14.15"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stof("1.5x"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StodDiffEmptyAndNonNumeric)
{
  // atof("") -> 0.0; safe_stod("") -> throws
  EXPECT_THROW(BRAINSTools::safe_stod(""), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("   "), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("abc"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StodDiffInfNanNotParsedOnAppleClang)
{
  // std::stod("inf") returns +infinity on all platforms.
  // safe_stod("inf") throws on Apple Clang / libc++ because
  // istringstream::operator>> does not parse "inf"/"nan" tokens when
  // imbued with locale::classic() on that platform.
  // This is a documented limitation of the istringstream approach.
  EXPECT_THROW(BRAINSTools::safe_stod("inf"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("INF"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("nan"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("NAN"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("infinity"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StodDiffOverflowThrowsInvalidArgument)
{
  // std::stod("1e9999") throws std::out_of_range.
  // safe_stod("1e9999") throws std::invalid_argument (istringstream sets
  // failbit on overflow and we always re-throw as invalid_argument).
  EXPECT_THROW(BRAINSTools::safe_stod("1e9999"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stod("-1e9999"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StodDiffUnderflowZeroOrThrows)
{
  // safe_stod("1e-9999") extreme-underflow behaviour is platform-dependent:
  //   Apple Clang/libc++: istringstream sets failbit -> throws std::invalid_argument
  //   GCC/libstdc++:      flushes to 0.0 without failbit -> returns 0.0
  // Both are acceptable.  NaN, Inf, or any non-zero result is not.
  bool   threw = false;
  double result = 0.0;
  try
  {
    result = BRAINSTools::safe_stod("1e-9999");
  }
  catch (const std::invalid_argument &)
  {
    threw = true;
  }
  if (!threw)
  {
    EXPECT_DOUBLE_EQ(result, 0.0) << "underflow must flush to zero, not " << result;
  }
}

// NOTE: Hexadecimal floating-point ("0x1.8p+1") is intentionally NOT tested
// here because the behaviour of istringstream::operator>> for hex floats is
// implementation-defined:
//   Apple Clang / libc++: parses hex floats -> safe_stod("0x1.8p+1") == 3.0
//   GCC / libstdc++:      does NOT parse   -> safe_stod("0x1.8p+1") throws
// Neither EXPECT_NEAR nor EXPECT_THROW would be portable across both
// compilers.  Callers that need to parse hexadecimal floating-point literals
// must use std::stod or std::from_chars directly.

// ===========================================================================
// 9. Intentional differences from std::stoi / atoi
//
//    safe_stoi is STRICTER in these ways:
//      a) Partial parse:    std::stoi("42abc") -> 42;  safe_stoi -> throws
//      b) Scientific:       std::stoi("1e3")   -> 1;   safe_stoi -> throws
//      c) Hex at base=10:   std::stoi("0x1a")  -> 0;   safe_stoi -> throws
//      d) Empty/non-numeric: atoi("") -> 0;             safe_stoi -> throws
//      e) Overflow:          atoi has UB;    safe_stoi throws invalid_argument
//                           (std::stoi throws out_of_range, but safe_stoi
//                            throws invalid_argument — different exception type)
// ===========================================================================

TEST(LocaleSafeConversions, StoiDiffPartialParseRejected)
{
  // std::stoi("42abc") returns 42; safe_stoi must throw.
  EXPECT_THROW(BRAINSTools::safe_stoi("42abc"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("42xyz"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("42.0"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StoiDiffScientificNotationRejected)
{
  // std::stoi("1e3") returns 1 (stops at 'e'); safe_stoi must throw.
  EXPECT_THROW(BRAINSTools::safe_stoi("1e3"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("1E3"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StoiDiffHexPrefixRejectedAtBase10)
{
  // std::stoi("0x1a") with default base=10 returns 0 (stops at 'x').
  // safe_stoi sees 'x' as trailing non-whitespace -> throws.
  EXPECT_THROW(BRAINSTools::safe_stoi("0x1a"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("0xFF"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StoiDiffEmptyAndNonNumericThrow)
{
  // atoi("") -> 0;   safe_stoi("") -> throws
  // atoi("abc") -> 0; safe_stoi("abc") -> throws
  EXPECT_THROW(BRAINSTools::safe_stoi(""), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("abc"), std::invalid_argument);
}

TEST(LocaleSafeConversions, StoiDiffOverflowThrowsInvalidArgument)
{
  // atoi has undefined behaviour on overflow.
  // safe_stoi throws (istringstream sets failbit; we throw invalid_argument).
  EXPECT_THROW(BRAINSTools::safe_stoi("9999999999"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("-9999999999"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoi("2147483648"), std::invalid_argument);  // INT_MAX + 1
  EXPECT_THROW(BRAINSTools::safe_stoi("-2147483649"), std::invalid_argument); // INT_MIN - 1
}

// ===========================================================================
// 10. safe_stoui-specific behaviour
//
//    f) Negative input "-1" wraps to UINT_MAX via unsigned operator>>
//       (same as strtoul("-1", ..., 10) on all POSIX platforms).
//       Callers that may receive signed strings should validate before calling.
// ===========================================================================

TEST(LocaleSafeConversions, SafeStouiLeadingWhitespaceAndSign)
{
  EXPECT_EQ(BRAINSTools::safe_stoui("  100"), 100u);
  EXPECT_EQ(BRAINSTools::safe_stoui("\t7"), 7u);
  EXPECT_EQ(BRAINSTools::safe_stoui("+42"), 42u);
}

TEST(LocaleSafeConversions, SafeStouiLeadingZerosAreDecimal)
{
  EXPECT_EQ(BRAINSTools::safe_stoui("007"), 7u);
  EXPECT_EQ(BRAINSTools::safe_stoui("042"), 42u);
}

TEST(LocaleSafeConversions, SafeStouiBoundaryValues) { EXPECT_EQ(BRAINSTools::safe_stoui("4294967295"), UINT_MAX); }

TEST(LocaleSafeConversions, SafeStouiOverflowThrows)
{
  EXPECT_THROW(BRAINSTools::safe_stoui("4294967296"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoui("9999999999"), std::invalid_argument);
}

TEST(LocaleSafeConversions, SafeStouiNegativeWrapsToUintMax)
{
  // Unsigned operator>> follows strtoul semantics: -1 wraps to UINT_MAX.
  // Callers must guard against negative inputs if wraparound is undesired.
  EXPECT_EQ(BRAINSTools::safe_stoui("-1"), static_cast<unsigned int>(-1));
}

TEST(LocaleSafeConversions, SafeStouiTrailingWhitespace)
{
  EXPECT_EQ(BRAINSTools::safe_stoui("7 "), 7u);
  EXPECT_EQ(BRAINSTools::safe_stoui("7\n"), 7u);
}

TEST(LocaleSafeConversions, SafeStouiInvalidInputThrows)
{
  EXPECT_THROW(BRAINSTools::safe_stoui(""), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoui("abc"), std::invalid_argument);
  EXPECT_THROW(BRAINSTools::safe_stoui("7abc"), std::invalid_argument);
}
