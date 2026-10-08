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

#include <cstdio>
#include <string>

#include "itkMacro.h"

/** Build an output filename from the single --outputVolumes pattern.
 *  "%s" is replaced by the volume type and "%d" by the list index, each at most
 *  once and in either order. "%d" accepts an optional zero-padded width such as
 *  "%03d". "%%" is a literal percent sign; any other use of '%' is rejected. */
inline std::string
ExpandOutputVolumeTemplate(const std::string & outputTemplate, const std::string & volumeType, unsigned int index)
{
  std::string result;
  bool        typeUsed = false;
  bool        indexUsed = false;
  for (std::string::size_type i = 0; i < outputTemplate.size(); ++i)
  {
    if (outputTemplate[i] != '%')
    {
      result += outputTemplate[i];
      continue;
    }
    ++i;
    const char next = (i < outputTemplate.size()) ? outputTemplate[i] : '\0';
    if (next == '%')
    {
      result += '%';
    }
    else if (next == 's' && !typeUsed)
    {
      result += volumeType;
      typeUsed = true;
    }
    else if (!indexUsed && (next == 'd' || (next == '0' && i + 2 < outputTemplate.size() &&
                                            outputTemplate[i + 1] >= '1' && outputTemplate[i + 1] <= '9')))
    {
      int width = 0;
      if (next == '0')
      {
        ++i;
        width = outputTemplate[i] - '0';
        if (i + 1 < outputTemplate.size() && outputTemplate[i + 1] >= '0' && outputTemplate[i + 1] <= '9')
        {
          ++i;
          width = width * 10 + (outputTemplate[i] - '0');
        }
        if (i + 1 >= outputTemplate.size() || outputTemplate[i + 1] != 'd')
        {
          itkGenericExceptionMacro("Invalid --outputVolumes pattern \"" << outputTemplate << "\"");
        }
        ++i;
      }
      char digits[128];
      std::snprintf(digits, sizeof(digits), "%0*u", width, index);
      result += digits;
      indexUsed = true;
    }
    else
    {
      itkGenericExceptionMacro("Invalid --outputVolumes pattern \""
                               << outputTemplate << "\": only one %s, one %d (or %0Nd), and %% are allowed");
    }
  }
  return result;
}

#endif // OutputVolumeTemplate_h
