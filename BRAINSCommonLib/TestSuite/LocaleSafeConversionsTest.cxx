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
/// \file LocaleSafeConversionsTest.cxx
/// \brief Unit tests for BRAINSCommonLib/LocaleSafeConversions.h
///
/// Test categories:
///   1. Round-trip correctness for all four conversion functions.
///   2. Locale independence: the functions always use the "C" locale even
///      when the global locale has been changed to one that uses comma as the
///      decimal separator.
///   3. Robustness against DICOM CSA header padding patterns:
///      a. Trailing newline:  "0.00000000\n"    (Siemens CSA format)
///      b. Trailing newline + null: "0.00000000\n\0\0"
///      c. Null-padded integers:   "0000\0\0\0\0"  (Philips headers)
///      d. Exact single chars:     "0", "2", "7"   (no trailing bytes)
///   4. Sentry-bug regression: strings consumed exactly to EOF must not fail.
///      Calling std::ws when eofbit is set causes the sentry to set failbit.
///   5. Error handling: invalid inputs must throw std::invalid_argument.
///   6. Behaviours matching std::stod (C++11) — safe_stod preserves these.
///   7. Behaviours matching std::stoi (C++11) — safe_stoi preserves these.
///   8. Intentional differences from std::stod / atof — safe_ is stricter.
///   9. Intentional differences from std::stoi / atoi — safe_ is stricter.
///  10. safe_stoui-specific behaviour (unsigned wraparound, limits).
///
/// ## Mapping from old locale-unsafe calls to safe_ equivalents
///
/// The following old calls were replaced in the BRAINSTools toolkit:
///
///   std::stod(str)    -> BRAINSTools::safe_stod(str)   (34 call sites)
///   std::stoi(str)    -> BRAINSTools::safe_stoi(str)   (26 call sites)
///   atof(str.c_str()) -> BRAINSTools::safe_stod(str)   (5 call sites)
///   atoi(str.c_str()) -> BRAINSTools::safe_stoi(str)   (1 call site)
///   std::atoi(...)    -> BRAINSTools::safe_stoi(str)   (1 call site)
///
/// The safe_ wrappers differ from the old calls in two important ways:
///   (a) They always use the C locale -- locale-independent by design.
///   (b) They reject partial parses (trailing non-whitespace throws).
///       std::stod("3.14abc") returns 3.14; safe_stod("3.14abc") throws.
///       atof("") returns 0.0;             safe_stod("") throws.
///   (c) They do NOT accept "inf"/"nan" on Apple Clang / libc++ because
///       std::istringstream::operator>> does not parse those tokens when
///       imbued with locale::classic() on that platform.

#include "LocaleSafeConversions.h"

#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <locale>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// Minimal test harness
// ---------------------------------------------------------------------------
static int g_failures = 0;

#define EXPECT_EQ(a, b)                                                                                       \
  do                                                                                                          \
  {                                                                                                           \
    if (!((a) == (b)))                                                                                        \
    {                                                                                                         \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  expected " << (b) << "  got " << (a) << "\n"; \
      ++g_failures;                                                                                           \
    }                                                                                                         \
  } while (false)

#define EXPECT_NEAR(a, b, tol)                                                                                 \
  do                                                                                                           \
  {                                                                                                            \
    if (std::fabs((double)(a) - (double)(b)) > (tol))                                                          \
    {                                                                                                          \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  expected ~" << (b) << "  got " << (a) << "\n"; \
      ++g_failures;                                                                                            \
    }                                                                                                          \
  } while (false)

#define EXPECT_THROWS(expr, ExcType)                                                                 \
  do                                                                                                 \
  {                                                                                                  \
    bool _threw = false;                                                                             \
    try                                                                                              \
    {                                                                                                \
      (expr);                                                                                        \
    }                                                                                                \
    catch (const ExcType &)                                                                          \
    {                                                                                                \
      _threw = true;                                                                                 \
    }                                                                                                \
    if (!_threw)                                                                                     \
    {                                                                                                \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  expected " #ExcType " not thrown\n"; \
      ++g_failures;                                                                                  \
    }                                                                                                \
  } while (false)

// Construct a std::string with embedded null bytes from a char array literal.
// sizeof(arr)-1 excludes the implicit C-string terminator that the compiler
// appends, so the resulting string contains exactly the bytes written.
#define STR_WITH_NULLS(arr) std::string((arr), sizeof(arr) - 1)

// ---------------------------------------------------------------------------
// 1. Round-trip correctness
// ---------------------------------------------------------------------------
static void
TestRoundTrip()
{
  // safe_stod
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14159"), 3.14159, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("-2.718"), -2.718, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("0.0"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("1e10"), 1e10, 1.0);
  EXPECT_NEAR(BRAINSTools::safe_stod("-1.5e-3"), -1.5e-3, 1e-10);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.2303285706216762"), 1.2303285706216762, 1e-12); // BCD real-world value

  // safe_stof
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5"), 1.5f, 1e-5f);
  EXPECT_NEAR(BRAINSTools::safe_stof("-0.25"), -0.25f, 1e-5f);

  // safe_stoi
  EXPECT_EQ(BRAINSTools::safe_stoi("42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("-7"), -7);
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);

  // safe_stoui
  EXPECT_EQ(BRAINSTools::safe_stoui("100"), 100u);
  EXPECT_EQ(BRAINSTools::safe_stoui("0"), 0u);
}

// ---------------------------------------------------------------------------
// 2. Locale independence
//
// Locales that use ',' as the decimal separator (and '.' or narrow-no-break
// space as the thousands separator) include every major European, Latin
// American, and Slavic locale.  A representative survey:
//
//   Germanic:  de_DE, de_AT, de_CH, de_LI, nl_NL, nl_BE, af_ZA
//   Romance:   fr_FR, fr_BE, fr_CA, fr_CH, es_ES, es_AR, it_IT, it_CH,
//              pt_PT, pt_BR, ca_ES, ro_RO
//   Nordic:    da_DK, sv_SE, sv_FI, fi_FI, nb_NO
//   Slavic:    pl_PL, cs_CZ, sk_SK, hr_HR, sl_SI, bg_BG, ru_RU, uk_UA
//   Baltic:    lv_LV, lt_LT, et_EE
//   Other:     hu_HU, tr_TR, el_GR, id_ID, ms_MY, vi_VN
//
// The test iterates ALL candidates in the list below, exercising each
// available locale independently so that no single locale acts as a proxy
// for all the others.  Unavailable locales are silently skipped.
// ---------------------------------------------------------------------------
static void
TestLocaleIndependence()
{
  // Comprehensive list of locales known to use ',' as decimal separator.
  // Both UTF-8 and legacy encodings are included so the test works on
  // systems that only have one or the other installed.
  static const char * const kCommaDecimalLocales[] = {
    // Germanic
    "de_DE.UTF-8",
    "de_DE",
    "de_AT.UTF-8",
    "de_AT",
    "de_CH.UTF-8",
    "de_CH",
    "nl_NL.UTF-8",
    "nl_NL",
    "nl_BE.UTF-8",
    "nl_BE",
    "af_ZA.UTF-8",
    "af_ZA",
    // Romance
    "fr_FR.UTF-8",
    "fr_FR",
    "fr_BE.UTF-8",
    "fr_BE",
    "fr_CA.UTF-8",
    "fr_CA",
    "fr_CH.UTF-8",
    "fr_CH",
    "es_ES.UTF-8",
    "es_ES",
    "es_AR.UTF-8",
    "es_AR",
    "it_IT.UTF-8",
    "it_IT",
    "it_CH.UTF-8",
    "it_CH",
    "pt_PT.UTF-8",
    "pt_PT",
    "pt_BR.UTF-8",
    "pt_BR",
    "ca_ES.UTF-8",
    "ca_ES",
    "ro_RO.UTF-8",
    "ro_RO",
    // Nordic
    "da_DK.UTF-8",
    "da_DK",
    "sv_SE.UTF-8",
    "sv_SE",
    "fi_FI.UTF-8",
    "fi_FI",
    "nb_NO.UTF-8",
    "nb_NO",
    // Slavic
    "pl_PL.UTF-8",
    "cs_CZ.UTF-8",
    "cs_CZ",
    "sk_SK.UTF-8",
    "sk_SK",
    "hr_HR.UTF-8",
    "hr_HR",
    "bg_BG.UTF-8",
    "bg_BG",
    "ru_RU.UTF-8",
    "uk_UA.UTF-8",
    "uk_UA",
    // Baltic
    "lv_LV.UTF-8",
    "lt_LT.UTF-8",
    "et_EE.UTF-8",
    // Other
    "hu_HU.UTF-8",
    "hu_HU",
    "el_GR.UTF-8",
    "tr_TR.UTF-8",
  };

  const std::locale original = std::locale::global(std::locale::classic());

  int testedCount = 0;

  for (const char * name : kCommaDecimalLocales)
  {
    // Skip duplicate encoding variants once the base locale has been tested.
    // We identify a "base" by the part before the first '.'; if we've already
    // tested that base under a different encoding, skip to avoid redundancy
    // while still exercising both encodings when both are present.
    try
    {
      std::locale::global(std::locale(name));
    }
    catch (const std::runtime_error &)
    {
      continue; // locale not installed on this system -- skip
    }

    ++testedCount;
    std::cout << "  [locale independence] testing " << name << "\n";

    // Under a comma-decimal locale, safe_ must still parse dot-decimal.
    EXPECT_NEAR(BRAINSTools::safe_stod("3.14"), 3.14, 1e-5);
    EXPECT_NEAR(BRAINSTools::safe_stof("1.5"), 1.5f, 1e-5f);
    EXPECT_EQ(BRAINSTools::safe_stoi("99"), 99);
    EXPECT_EQ(BRAINSTools::safe_stoui("7"), 7u);

    // Comma-decimal strings must fail (safe_ only accepts dot-decimal)
    EXPECT_THROWS(BRAINSTools::safe_stod("3,14"), std::invalid_argument);

    // Restore classic locale before trying the next candidate
    std::locale::global(std::locale::classic());
  }

  if (testedCount == 0)
    std::cout << "  [locale independence] WARNING: no comma-decimal locale "
                 "available on this system; locale guard untested\n";
  else
    std::cout << "  [locale independence] tested " << testedCount << " comma-decimal locale(s)\n";

  std::locale::global(original);
}

// ---------------------------------------------------------------------------
// 3. DICOM CSA header padding patterns
//    Siemens: values stored as ASCII padded with '\n' then '\0' to 4-byte boundary
//    Philips: integer fields stored as "0000\0\0\0\0"
// ---------------------------------------------------------------------------
static void
TestDicomPadding()
{
  // --- 3a: trailing newline (Siemens CSA float values) ---
  // "0.00000000\n" -- 11 bytes including newline
  EXPECT_NEAR(BRAINSTools::safe_stod("0.00000000\n"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14 "), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5\n"), 1.5f, 1e-5f);
  EXPECT_EQ(BRAINSTools::safe_stoi("42\n"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoui("7\t"), 7u);

  // --- 3b: newline followed by null-padding (Siemens CSA, 4-byte aligned) ---
  // "0.00000000\n\0\0" -- 13 bytes; null bytes are alignment padding
  EXPECT_NEAR(BRAINSTools::safe_stod(STR_WITH_NULLS("0.00000000\n\0\0")), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod(STR_WITH_NULLS("0.00000000\n\0\0\0")), 0.0, 1e-15);

  // --- 3c: null-padded integers (Philips header format) ---
  // "0000\0\0\0\0" -- 8 bytes; integer value is "0000"
  EXPECT_EQ(BRAINSTools::safe_stoi(STR_WITH_NULLS("0000\0\0\0\0")), 0);
  EXPECT_EQ(BRAINSTools::safe_stoui(STR_WITH_NULLS("0000\0")), 0u);
  // SpaceThicknessDiff case: "2\0\0\0"
  EXPECT_EQ(BRAINSTools::safe_stoi(STR_WITH_NULLS("2\0\0\0")), 2);

  // --- 3d: exact single chars with NO trailing bytes ---
  // These triggered the std::ws sentry bug: parsing "0" or "2" consumed the
  // entire stream, setting eofbit.  Calling std::ws then caused the sentry
  // to set failbit (is.good() returns false when eofbit is set).
  EXPECT_NEAR(BRAINSTools::safe_stod("0"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.2303285706216762"), 1.2303285706216762, 1e-12);
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("0000"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("2"), 2);
  EXPECT_EQ(BRAINSTools::safe_stoui("0"), 0u);
}

// ---------------------------------------------------------------------------
// 4. Sentry-bug regression
//    Verify that strings exactly consumed to EOF pass, not fail.
// ---------------------------------------------------------------------------
static void
TestSentryBugRegression()
{
  // All of these were failing with the naive `iss >> std::ws` approach because
  // eofbit was set after parsing, and std::ws's sentry then set failbit.
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14159"), 3.14159, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("-2.718"), -2.718, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("0"), 0.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5"), 1.5f, 1e-5f);
  EXPECT_EQ(BRAINSTools::safe_stoi("42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("0000"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoui("100"), 100u);
}

// ---------------------------------------------------------------------------
// 5. Error handling
// ---------------------------------------------------------------------------
static void
TestErrorHandling()
{
  // Empty string (after null truncation, nothing to parse)
  EXPECT_THROWS(BRAINSTools::safe_stod(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stof(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoui(""), std::invalid_argument);

  // Null-only string (truncated to empty string, nothing to parse)
  EXPECT_THROWS(BRAINSTools::safe_stod(STR_WITH_NULLS("\0\0\0")), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi(STR_WITH_NULLS("\0")), std::invalid_argument);

  // Whitespace only -- parsing fails because no number found
  EXPECT_THROWS(BRAINSTools::safe_stod("   "), std::invalid_argument);

  // Trailing non-whitespace characters must still throw
  EXPECT_THROWS(BRAINSTools::safe_stod("3.14abc"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stof("1.0f"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("42xyz"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoui("5u"), std::invalid_argument);

  // Entirely non-numeric
  EXPECT_THROWS(BRAINSTools::safe_stod("nan_not_accepted"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("abc"), std::invalid_argument);

  // Multiple decimal points -- second '.' is trailing non-whitespace
  EXPECT_THROWS(BRAINSTools::safe_stod("3.14.15"), std::invalid_argument);

  // Embedded space -- digit then space then digit is not a single number
  EXPECT_THROWS(BRAINSTools::safe_stoi("4 2"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 6. Behaviours matching std::stod (C++11)
//
// std::stod documentation (cppreference):
//   - Discards leading whitespace (std::isspace).
//   - Accepts optional sign (+/-).
//   - Accepts scientific notation (e/E with optional sign in exponent).
//   - Accepts hexadecimal floating-point (0x prefix, binary exponent p/P).
//   - Accepts leading zeros in mantissa (treated as decimal digits).
//   - Throws std::invalid_argument if no conversion performed.
//   - Throws std::out_of_range on overflow.
//   - Accepts "inf"/"infinity"/"nan" on most platforms (NOT via istringstream
//     on Apple Clang/libc++ -- see TestStdStodDifferences below).
//
// safe_stod preserves all behaviours listed above EXCEPT:
//   - Partial parse: std::stod("3.14abc") -> 3.14; safe_stod throws.
//   - inf/nan: not parsed by istringstream on Apple Clang/libc++.
// ---------------------------------------------------------------------------
static void
TestStdStodBehaviorMatch()
{
  // Leading whitespace (std::stod skips it; safe_stod must too)
  EXPECT_NEAR(BRAINSTools::safe_stod("  3.14"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("\t-1.0"), -1.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stod("\n2.5"), 2.5, 1e-5);

  // Optional plus sign (std::stod accepts; safe_stod must too)
  EXPECT_NEAR(BRAINSTools::safe_stod("+3.14"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stof("+1.5"), 1.5f, 1e-5f);

  // Scientific notation -- all four forms (e/E, optional sign in exponent)
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5e3"), 1500.0, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5E3"), 1500.0, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5e+3"), 1500.0, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1.5e-3"), 0.0015, 1e-9);
  EXPECT_NEAR(BRAINSTools::safe_stod("1e0"), 1.0, 1e-15);
  EXPECT_NEAR(BRAINSTools::safe_stof("2.5E2"), 250.0f, 1e-3f);

  // Hexadecimal floating-point (C99/C++11: 0x1.8p+1 = 1.5 * 2^1 = 3.0)
  // std::stod supports these; istringstream with locale::classic() also does.
  EXPECT_NEAR(BRAINSTools::safe_stod("0x1.8p+1"), 3.0, 1e-12);
  EXPECT_NEAR(BRAINSTools::safe_stod("0x0.8p+0"), 0.5, 1e-12);

  // Leading zeros in mantissa (decimal digits, not octal)
  EXPECT_NEAR(BRAINSTools::safe_stod("007.5"), 7.5, 1e-10);
  EXPECT_NEAR(BRAINSTools::safe_stod("0042.5"), 42.5, 1e-10);
  EXPECT_NEAR(BRAINSTools::safe_stod("00.5"), 0.5, 1e-10);

  // Trailing whitespace accepted (consistent with strtod semantics)
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14 "), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14\t"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14\n"), 3.14, 1e-5);
}

// ---------------------------------------------------------------------------
// 7. Behaviours matching std::stoi (C++11)
//
// std::stoi documentation (cppreference):
//   - Discards leading whitespace (std::isspace).
//   - Accepts optional sign (+/-) for signed types.
//   - Base 10 by default; "0x" prefix does NOT trigger hex (stops at 'x').
//   - Leading zeros are decimal digits, NOT octal (unlike C integer literals).
//     "007" -> 7  (decimal, NOT 7 as octal -- they happen to equal here)
//     "042" -> 42 (decimal, NOT 34 as octal)
//   - Throws std::invalid_argument if no conversion performed.
//   - Throws std::out_of_range if result exceeds int range.
//   - Accepts trailing non-numeric content (partial parse): "42abc" -> 42.
//     safe_stoi REJECTS partial parses -- this is an intentional improvement.
//
// safe_stoi preserves all behaviours listed above EXCEPT partial parse.
// ---------------------------------------------------------------------------
static void
TestStdStoiBehaviorMatch()
{
  // Leading whitespace (std::stoi skips it; safe_stoi must too)
  EXPECT_EQ(BRAINSTools::safe_stoi("  42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("\t-7"), -7);
  EXPECT_EQ(BRAINSTools::safe_stoi("\n100"), 100);

  // Optional plus sign
  EXPECT_EQ(BRAINSTools::safe_stoi("+42"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("+0"), 0);

  // Leading zeros -- base 10, NOT octal
  // std::stoi("007") == 7 (decimal)
  // std::stoi("042") == 42 (decimal, NOT 34 which is octal 042)
  EXPECT_EQ(BRAINSTools::safe_stoi("007"), 7);
  EXPECT_EQ(BRAINSTools::safe_stoi("042"), 42); // NOT 34 (octal 042)
  EXPECT_EQ(BRAINSTools::safe_stoi("0042"), 42);

  // Zero variants -- all should return 0
  EXPECT_EQ(BRAINSTools::safe_stoi("0"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("00"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("000"), 0);
  EXPECT_EQ(BRAINSTools::safe_stoi("0000"), 0);

  // Boundary values -- match std::stoi exactly
  EXPECT_EQ(BRAINSTools::safe_stoi("2147483647"), INT_MAX);
  EXPECT_EQ(BRAINSTools::safe_stoi("-2147483648"), INT_MIN);

  // Trailing whitespace accepted
  EXPECT_EQ(BRAINSTools::safe_stoi("42 "), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("42\n"), 42);
  EXPECT_EQ(BRAINSTools::safe_stoi("42\t"), 42);
}

// ---------------------------------------------------------------------------
// 8. Intentional differences from std::stod / atof
//
// safe_stod is STRICTER than std::stod and atof in these ways:
//
//   a) Partial parse: std::stod("3.14abc") -> 3.14 (stops, no error)
//                     atof("3.14abc")      -> 3.14 (stops, no error)
//                     safe_stod("3.14abc") -> throws std::invalid_argument
//
//   b) Empty / whitespace / non-numeric:
//                     atof("")             -> 0.0 (no error)
//                     atof("abc")          -> 0.0 (no error)
//                     safe_stod("")        -> throws std::invalid_argument
//                     safe_stod("abc")     -> throws std::invalid_argument
//
//   c) inf / nan: std::stod("inf")  -> +infinity (platform-standard)
//                 safe_stod("inf")  -> throws (istringstream on Apple Clang
//                 does not parse "inf"/"nan" when imbued with locale::classic())
//                 This is a known platform limitation of the istringstream
//                 approach; none of the migrated call sites encounter
//                 "inf"/"nan" values in practice.
//
//   d) Overflow: std::stod("1e9999") -> throws std::out_of_range
//                safe_stod("1e9999") -> throws std::invalid_argument
//                (istringstream sets failbit on overflow; we always throw
//                invalid_argument regardless of the failure cause)
// ---------------------------------------------------------------------------
static void
TestStdStodDifferences()
{
  // --- 8a: partial parse rejected (safer than std::stod / atof) ---
  // std::stod("3.14abc") returns 3.14; safe_stod must throw.
  EXPECT_THROWS(BRAINSTools::safe_stod("3.14abc"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("3.14.15"), std::invalid_argument); // second '.' is trailing
  EXPECT_THROWS(BRAINSTools::safe_stof("1.5x"), std::invalid_argument);

  // --- 8b: empty / whitespace / non-numeric (atof returns 0; safe_ throws) ---
  // atof("") -> 0.0 -- safe_stod("") -> throws
  EXPECT_THROWS(BRAINSTools::safe_stod(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("   "), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("abc"), std::invalid_argument);

  // --- 8c: inf / nan not parsed by istringstream on Apple Clang / libc++ ---
  // std::stod("inf") returns +infinity.
  // safe_stod("inf") throws on platforms where istringstream does not
  // parse "inf"/"nan" (Apple Clang, and potentially other implementations).
  // This is a documented platform limitation of the istringstream approach.
  EXPECT_THROWS(BRAINSTools::safe_stod("inf"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("INF"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("nan"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("NAN"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("infinity"), std::invalid_argument);

  // --- 8d: overflow / underflow (exception type differs from std::stod) ---
  // std::stod("1e9999") throws std::out_of_range.
  // safe_stod("1e9999") throws std::invalid_argument (istringstream sets
  // failbit on overflow; we always re-throw as invalid_argument).
  EXPECT_THROWS(BRAINSTools::safe_stod("1e9999"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("-1e9999"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stod("1e-9999"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 9. Intentional differences from std::stoi / atoi
//
// safe_stoi is STRICTER than std::stoi and atoi in these ways:
//
//   a) Partial parse: std::stoi("42abc") -> 42 (stops at 'a', no error)
//                     atoi("42abc")      -> 42 (stops at 'a', no error)
//                     safe_stoi("42abc") -> throws std::invalid_argument
//
//   b) Embedded '.' or 'e' (scientific notation treated as partial):
//                     std::stoi("1e3")   -> 1 (stops at 'e', no error)
//                     safe_stoi("1e3")   -> throws ('e' is trailing non-ws)
//
//   c) Hex prefix with base=10:
//                     std::stoi("0x1a")  -> 0 (stops at 'x', returns 0)
//                     safe_stoi("0x1a")  -> throws ('x' is trailing non-ws)
//
//   d) Empty / non-numeric:
//                     atoi("")           -> 0 (no error, no throw)
//                     atoi("abc")        -> 0 (no error, no throw)
//                     safe_stoi("")      -> throws std::invalid_argument
//                     safe_stoi("abc")   -> throws std::invalid_argument
//
//   e) Overflow: atoi has UNDEFINED BEHAVIOUR on overflow.
//                safe_stoi throws (inheriting istringstream overflow detection).
//                std::stoi throws std::out_of_range; safe_stoi throws
//                std::invalid_argument (different exception type).
// ---------------------------------------------------------------------------
static void
TestStdStoiDifferences()
{
  // --- 9a: partial parse rejected ---
  // std::stoi("42abc") returns 42; safe_stoi must throw.
  EXPECT_THROWS(BRAINSTools::safe_stoi("42abc"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("42xyz"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("42.0"), std::invalid_argument); // '.' is trailing

  // --- 9b: scientific notation rejected (not an integer) ---
  // std::stoi("1e3") returns 1 (stops at 'e'); safe_stoi must throw.
  EXPECT_THROWS(BRAINSTools::safe_stoi("1e3"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("1E3"), std::invalid_argument);

  // --- 9c: hex prefix rejected at base=10 ---
  // std::stoi("0x1a") with default base=10 returns 0 (stops at 'x').
  // safe_stoi sees 'x' as trailing non-whitespace -> throws.
  EXPECT_THROWS(BRAINSTools::safe_stoi("0x1a"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("0xFF"), std::invalid_argument);

  // --- 9d: empty / non-numeric (atoi returns 0; safe_ throws) ---
  EXPECT_THROWS(BRAINSTools::safe_stoi(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("abc"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("   "), std::invalid_argument);

  // --- 9e: overflow (atoi: undefined; safe_stoi: throws invalid_argument) ---
  // 9999999999 exceeds INT_MAX (2147483647); safe_stoi must throw.
  EXPECT_THROWS(BRAINSTools::safe_stoi("9999999999"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("-9999999999"), std::invalid_argument);
  // Edge: INT_MAX+1 overflows
  EXPECT_THROWS(BRAINSTools::safe_stoi("2147483648"), std::invalid_argument);
  // Edge: INT_MIN-1 overflows
  EXPECT_THROWS(BRAINSTools::safe_stoi("-2147483649"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// 10. safe_stoui-specific behaviour
//
// safe_stoui wraps std::istringstream::operator>> for unsigned int.
// Key behaviours derived from the underlying C++ standard semantics:
//
//   a) Leading whitespace accepted (consistent with strtoul).
//   b) Plus sign accepted.
//   c) Leading zeros are decimal (same as strtoul with base=10).
//   d) UINT_MAX ("4294967295") accepted.
//   e) Overflow (4294967296+) throws std::invalid_argument.
//   f) Negative input ("-1"): std::istringstream operator>> for unsigned int
//      follows strtoul semantics -- wraps to UINT_MAX.  This matches the
//      behaviour of strtoul("-1", ..., 10) on all POSIX platforms.
//      Callers that may receive signed-looking strings should validate
//      input before calling safe_stoui.
// ---------------------------------------------------------------------------
static void
TestSafeStouiBehavior()
{
  // Leading whitespace
  EXPECT_EQ(BRAINSTools::safe_stoui("  100"), 100u);
  EXPECT_EQ(BRAINSTools::safe_stoui("\t7"), 7u);

  // Plus sign
  EXPECT_EQ(BRAINSTools::safe_stoui("+42"), 42u);

  // Leading zeros (decimal)
  EXPECT_EQ(BRAINSTools::safe_stoui("007"), 7u);
  EXPECT_EQ(BRAINSTools::safe_stoui("042"), 42u);

  // UINT_MAX boundary
  EXPECT_EQ(BRAINSTools::safe_stoui("4294967295"), UINT_MAX);

  // Overflow -- one above UINT_MAX must throw
  EXPECT_THROWS(BRAINSTools::safe_stoui("4294967296"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoui("9999999999"), std::invalid_argument);

  // Negative input -- wraps to UINT_MAX via unsigned operator>> (strtoul semantics).
  // This is an inherent property of the unsigned extraction operator; callers
  // should guard against negative inputs if wraparound is not desired.
  EXPECT_EQ(BRAINSTools::safe_stoui("-1"), static_cast<unsigned int>(-1));

  // Trailing whitespace accepted
  EXPECT_EQ(BRAINSTools::safe_stoui("7 "), 7u);
  EXPECT_EQ(BRAINSTools::safe_stoui("7\n"), 7u);

  // Partial / non-numeric -- must throw
  EXPECT_THROWS(BRAINSTools::safe_stoui(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoui("abc"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoui("7abc"), std::invalid_argument);
}

// ---------------------------------------------------------------------------

int
main(int, char **)
{
  std::cout << "LocaleSafeConversionsTest: running...\n";

  TestRoundTrip();
  std::cout << "  round-trip: done\n";

  TestLocaleIndependence();
  std::cout << "  locale independence: done\n";

  TestDicomPadding();
  std::cout << "  DICOM padding patterns: done\n";

  TestSentryBugRegression();
  std::cout << "  sentry-bug regression: done\n";

  TestErrorHandling();
  std::cout << "  error handling: done\n";

  TestStdStodBehaviorMatch();
  std::cout << "  std::stod behavior match: done\n";

  TestStdStoiBehaviorMatch();
  std::cout << "  std::stoi behavior match: done\n";

  TestStdStodDifferences();
  std::cout << "  std::stod intentional differences: done\n";

  TestStdStoiDifferences();
  std::cout << "  std::stoi intentional differences: done\n";

  TestSafeStouiBehavior();
  std::cout << "  safe_stoui-specific behaviour: done\n";

  if (g_failures > 0)
  {
    std::cerr << "FAILED: " << g_failures << " assertion(s) failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASSED\n";
  return EXIT_SUCCESS;
}
