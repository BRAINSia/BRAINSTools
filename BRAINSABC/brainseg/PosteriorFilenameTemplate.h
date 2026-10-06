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
#ifndef PosteriorFilenameTemplate_h
#define PosteriorFilenameTemplate_h

#include <string>

#include "itkMacro.h"

/** Build a posterior output filename from the --posteriorTemplate value.
 *  "%s" is replaced by the prior name and "%%" by a literal percent sign;
 *  any other use of '%' is rejected. */
inline std::string
ExpandPosteriorTemplate(const std::string & posteriorTemplate, const std::string & priorName)
{
  std::string result;
  bool        substituted = false;
  for (std::string::size_type i = 0; i < posteriorTemplate.size(); ++i)
  {
    const char c = posteriorTemplate[i];
    if (c != '%')
    {
      result += c;
      continue;
    }
    const char next = (i + 1 < posteriorTemplate.size()) ? posteriorTemplate[i + 1] : '\0';
    if (next == '%')
    {
      result += '%';
    }
    else if (next == 's' && !substituted)
    {
      result += priorName;
      substituted = true;
    }
    else
    {
      itkGenericExceptionMacro("Invalid --posteriorTemplate \"" << posteriorTemplate
                                                                << "\": only one %s and %% are allowed");
    }
    ++i;
  }
  return result;
}

#endif // PosteriorFilenameTemplate_h
