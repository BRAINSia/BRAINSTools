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
/// Addresses GitHub issue #403.

#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

namespace BRAINSTools
{

/// \brief Locale-independent string to double conversion.
/// \param s The string to parse.
/// \return The parsed double value.
/// \throws std::invalid_argument if the string cannot be parsed.
inline double
safe_stod(const std::string & s)
{
  std::istringstream iss(s);
  iss.imbue(std::locale::classic());
  double value{};
  iss >> value;
  if (iss.fail() || !iss.eof())
  {
    throw std::invalid_argument("safe_stod: cannot parse '" + s + "'");
  }
  return value;
}

/// \brief Locale-independent string to float conversion.
/// \param s The string to parse.
/// \return The parsed float value.
/// \throws std::invalid_argument if the string cannot be parsed.
inline float
safe_stof(const std::string & s)
{
  std::istringstream iss(s);
  iss.imbue(std::locale::classic());
  float value{};
  iss >> value;
  if (iss.fail() || !iss.eof())
  {
    throw std::invalid_argument("safe_stof: cannot parse '" + s + "'");
  }
  return value;
}

/// \brief Locale-independent string to int conversion.
/// \param s The string to parse.
/// \return The parsed int value.
/// \throws std::invalid_argument if the string cannot be parsed.
inline int
safe_stoi(const std::string & s)
{
  std::istringstream iss(s);
  iss.imbue(std::locale::classic());
  int value{};
  iss >> value;
  if (iss.fail() || !iss.eof())
  {
    throw std::invalid_argument("safe_stoi: cannot parse '" + s + "'");
  }
  return value;
}

/// \brief Locale-independent string to unsigned int conversion.
/// \param s The string to parse.
/// \return The parsed unsigned int value.
/// \throws std::invalid_argument if the string cannot be parsed.
inline unsigned int
safe_stoui(const std::string & s)
{
  std::istringstream iss(s);
  iss.imbue(std::locale::classic());
  unsigned int value{};
  iss >> value;
  if (iss.fail() || !iss.eof())
  {
    throw std::invalid_argument("safe_stoui: cannot parse '" + s + "'");
  }
  return value;
}

} // namespace BRAINSTools

#endif // LocaleSafeConversions_h
