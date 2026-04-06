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
#ifndef LocaleSafeConversions_h
#define LocaleSafeConversions_h

/// \file LocaleSafeConversions.h
/// \brief Locale-independent string-to-number conversion utilities.
///
/// The standard functions std::stod(), std::stof(), std::stoi(), and atof()
/// are locale-dependent: they interpret ',' vs '.' as the decimal separator
/// based on the process locale.  International users whose locale uses comma
/// as the decimal separator will see silent data corruption or parse failures.
///
/// This header provides drop-in replacements that always parse numbers using
/// the "C" (POSIX) locale, regardless of the user's environment.
///
/// Implementation uses std::istringstream imbued with std::locale::classic(),
/// which is portable across all C++17 compilers (GCC, Clang, Apple Clang,
/// MSVC) without requiring platform-specific extensions such as strtod_l or
/// std::from_chars for floating-point (which Apple Clang did not support
/// until very recent versions).
///
/// ## Robustness design notes
///
/// ### Trailing whitespace (DICOM CSA header data)
/// Siemens DICOM CSA header fields store values as fixed-width ASCII strings
/// padded with trailing newlines and/or null bytes, e.g. "0.00000000\n\0\0".
/// We must accept trailing whitespace (matching strtod/stod semantics) while
/// still rejecting genuinely invalid trailing non-whitespace characters.
///
/// ### Null-byte truncation
/// C-string semantics stop at the first '\0'.  Binary DICOM data may embed
/// null bytes as alignment padding after the meaningful content.  We truncate
/// the input at the first null byte before parsing, matching the behaviour of
/// strtod(s.c_str(), ...) which naturally stops at '\0'.
///
/// ### std::ws sentry pitfall
/// Calling `iss >> std::ws` when `eofbit` is already set causes the sentry
/// inside std::ws to call `setstate(failbit)` (because `is.good()` returns
/// false whenever any error bit is set, even eofbit alone).  This would make
/// valid inputs like "3.14" fail after exact consumption.  We guard the
/// `std::ws` call with `if (!iss.eof())` to avoid this.
///
/// Addresses GitHub issue #403.

#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

namespace BRAINSTools
{

namespace detail
{
/// Truncate \p s at the first embedded null byte.
/// DICOM CSA headers null-pad values to 4-byte alignment; C-string strtod
/// stops at '\0' naturally.  We replicate that behaviour.
inline std::string
null_truncate(const std::string & s)
{
  const auto pos = s.find('\0');
  return (pos != std::string::npos) ? s.substr(0, pos) : s;
}

/// After successfully extracting a numeric value, verify that only
/// whitespace (or nothing) remains.  Calling `iss >> std::ws` directly
/// when eofbit is already set triggers the sentry to set failbit, so we
/// guard with !iss.eof() first.
inline bool
only_whitespace_remains(std::istringstream & iss)
{
  if (iss.eof())
  {
    return true; // entire string was consumed — nothing left to check
  }
  iss >> std::ws; // sentry is safe here: eofbit was not set entering this call
  return iss.eof();
}
} // namespace detail

/// \brief Locale-independent string to double conversion.
/// \param s The string to parse (null bytes and trailing whitespace tolerated).
/// \return The parsed double value.
/// \throws std::invalid_argument if the string contains non-numeric content
///         beyond optional trailing whitespace.
inline double
safe_stod(const std::string & s)
{
  const std::string  trimmed = detail::null_truncate(s);
  std::istringstream iss(trimmed);
  iss.imbue(std::locale::classic());
  double value{};
  iss >> value;
  if (iss.fail() || !detail::only_whitespace_remains(iss))
  {
    throw std::invalid_argument("safe_stod: cannot parse '" + trimmed + "'");
  }
  return value;
}

/// \brief Locale-independent string to float conversion.
/// \param s The string to parse (null bytes and trailing whitespace tolerated).
/// \return The parsed float value.
/// \throws std::invalid_argument if the string contains non-numeric content
///         beyond optional trailing whitespace.
inline float
safe_stof(const std::string & s)
{
  const std::string  trimmed = detail::null_truncate(s);
  std::istringstream iss(trimmed);
  iss.imbue(std::locale::classic());
  float value{};
  iss >> value;
  if (iss.fail() || !detail::only_whitespace_remains(iss))
  {
    throw std::invalid_argument("safe_stof: cannot parse '" + trimmed + "'");
  }
  return value;
}

/// \brief Locale-independent string to int conversion.
/// \param s The string to parse (null bytes and trailing whitespace tolerated).
/// \return The parsed int value.
/// \throws std::invalid_argument if the string contains non-numeric content
///         beyond optional trailing whitespace.
inline int
safe_stoi(const std::string & s)
{
  const std::string  trimmed = detail::null_truncate(s);
  std::istringstream iss(trimmed);
  iss.imbue(std::locale::classic());
  int value{};
  iss >> value;
  if (iss.fail() || !detail::only_whitespace_remains(iss))
  {
    throw std::invalid_argument("safe_stoi: cannot parse '" + trimmed + "'");
  }
  return value;
}

/// \brief Locale-independent string to unsigned int conversion.
/// \param s The string to parse (null bytes and trailing whitespace tolerated).
/// \return The parsed unsigned int value.
/// \throws std::invalid_argument if the string contains non-numeric content
///         beyond optional trailing whitespace.
inline unsigned int
safe_stoui(const std::string & s)
{
  const std::string  trimmed = detail::null_truncate(s);
  std::istringstream iss(trimmed);
  iss.imbue(std::locale::classic());
  unsigned int value{};
  iss >> value;
  if (iss.fail() || !detail::only_whitespace_remains(iss))
  {
    throw std::invalid_argument("safe_stoui: cannot parse '" + trimmed + "'");
  }
  return value;
}

} // namespace BRAINSTools

#endif // LocaleSafeConversions_h
