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
#ifndef OutputVolumeTemplate_h
#define OutputVolumeTemplate_h

#include <string>

#include "itkMacro.h"

/** Build an output filename from the single --outputVolumes pattern.
 *  "%s" is replaced by the volume type and "%d" by the list index, each at most
 *  once and in either order. "%d" accepts an optional zero-padded width such as
 *  "%03d". "%%" is a literal percent sign; any other use of '%' is rejected. */
inline std::string
ExpandOutputVolumeTemplate(const std::string & outputTemplate, const std::string & volumeType, unsigned int index)
{
  const auto invalid = [&outputTemplate]() {
    itkGenericExceptionMacro("Invalid --outputVolumes pattern \""
                             << outputTemplate << "\": only one %s, one %d (or %0Nd), and %% are allowed");
  };
  const auto isDigit = [](const char c) { return c >= '0' && c <= '9'; };

  std::string            result;
  bool                   typeUsed = false;
  bool                   indexUsed = false;
  std::string::size_type pos = 0;
  while (pos < outputTemplate.size())
  {
    const char c = outputTemplate[pos++];
    if (c != '%')
    {
      result += c;
      continue;
    }
    const char conversion = (pos < outputTemplate.size()) ? outputTemplate[pos++] : '\0';
    if (conversion == '%')
    {
      result += '%';
    }
    else if (conversion == 's' && !typeUsed)
    {
      result += volumeType;
      typeUsed = true;
    }
    else if (conversion == 'd' && !indexUsed)
    {
      result += std::to_string(index);
      indexUsed = true;
    }
    else if (conversion == '0' && !indexUsed)
    {
      // "%0Nd": N is one or two digits, the first not zero
      if (pos >= outputTemplate.size() || outputTemplate[pos] < '1' || outputTemplate[pos] > '9')
      {
        invalid();
      }
      std::string::size_type width = static_cast<std::string::size_type>(outputTemplate[pos++] - '0');
      if (pos < outputTemplate.size() && isDigit(outputTemplate[pos]))
      {
        width = width * 10 + static_cast<std::string::size_type>(outputTemplate[pos++] - '0');
      }
      if (pos >= outputTemplate.size() || outputTemplate[pos++] != 'd')
      {
        invalid();
      }
      const std::string digits = std::to_string(index);
      result.append(width > digits.size() ? width - digits.size() : 0, '0');
      result += digits;
      indexUsed = true;
    }
    else
    {
      invalid();
    }
  }
  return result;
}

#endif // OutputVolumeTemplate_h
