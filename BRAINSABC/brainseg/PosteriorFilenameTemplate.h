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
  std::string            result;
  bool                   substituted = false;
  std::string::size_type pos = 0;
  while (pos < posteriorTemplate.size())
  {
    const char c = posteriorTemplate[pos++];
    if (c != '%')
    {
      result += c;
      continue;
    }
    const char conversion = (pos < posteriorTemplate.size()) ? posteriorTemplate[pos++] : '\0';
    if (conversion == '%')
    {
      result += '%';
    }
    else if (conversion == 's' && !substituted)
    {
      result += priorName;
      substituted = true;
    }
    else
    {
      itkGenericExceptionMacro("Invalid --posteriorTemplate \"" << posteriorTemplate
                                                                << "\": only one %s and %% are allowed");
    }
  }
  return result;
}

#endif // PosteriorFilenameTemplate_h
