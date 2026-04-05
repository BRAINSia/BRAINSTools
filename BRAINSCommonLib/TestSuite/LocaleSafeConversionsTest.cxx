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
/// Tests cover:
///   1. Round-trip correctness for all four conversion functions.
///   2. Locale independence: the functions always use the "C" locale even
///      when the global locale has been changed to one that uses comma as the
///      decimal separator.
///   3. Error handling: std::invalid_argument is thrown for empty strings,
///      strings with trailing non-numeric characters, and strings that are
///      entirely non-numeric.

#include "LocaleSafeConversions.h"

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

// ---------------------------------------------------------------------------
// Test cases
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

static void
TestLocaleIndependence()
{
  // Change the global locale to one where ',' is the decimal separator if
  // available.  If the locale is unavailable the test still validates the
  // C-locale parsing path.
  const std::locale original = std::locale::global(std::locale::classic());

  bool changed = false;
  for (const char * name : { "de_DE.UTF-8", "de_DE", "fr_FR.UTF-8", "fr_FR" })
  {
    try
    {
      std::locale::global(std::locale(name));
      changed = true;
      break;
    }
    catch (const std::runtime_error &)
    {
      // locale not available on this system — try the next one
    }
  }

  if (changed)
  {
    std::cout << "  [locale independence] non-C locale active\n";
  }
  else
  {
    std::cout << "  [locale independence] no comma-decimal locale available; "
                 "testing C locale path only\n";
  }

  // Even under a comma-decimal locale the functions must parse dot-decimal.
  EXPECT_NEAR(BRAINSTools::safe_stod("3.14"), 3.14, 1e-5);
  EXPECT_NEAR(BRAINSTools::safe_stof("1.5"), 1.5f, 1e-5f);
  EXPECT_EQ(BRAINSTools::safe_stoi("99"), 99);
  EXPECT_EQ(BRAINSTools::safe_stoui("7"), 7u);

  // Restore original locale.
  std::locale::global(original);
}

static void
TestErrorHandling()
{
  // Empty string
  EXPECT_THROWS(BRAINSTools::safe_stod(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stof(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi(""), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoui(""), std::invalid_argument);

  // Whitespace only
  EXPECT_THROWS(BRAINSTools::safe_stod("   "), std::invalid_argument);

  // Trailing non-numeric characters
  EXPECT_THROWS(BRAINSTools::safe_stod("3.14abc"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stof("1.0f"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("42xyz"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoui("5u"), std::invalid_argument);

  // Entirely non-numeric
  EXPECT_THROWS(BRAINSTools::safe_stod("nan_not_accepted"), std::invalid_argument);
  EXPECT_THROWS(BRAINSTools::safe_stoi("abc"), std::invalid_argument);

  // Comma-decimal (must fail — functions only accept dot-decimal)
  EXPECT_THROWS(BRAINSTools::safe_stod("3,14"), std::invalid_argument);
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

  TestErrorHandling();
  std::cout << "  error handling: done\n";

  if (g_failures > 0)
  {
    std::cerr << "FAILED: " << g_failures << " assertion(s) failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASSED\n";
  return EXIT_SUCCESS;
}
