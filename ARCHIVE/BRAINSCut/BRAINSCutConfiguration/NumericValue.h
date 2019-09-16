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
#ifndef NumericValue_h
#define NumericValue_h
#include "ElementContainer.h"

template <typename T>
class NumericValue : public XMLContents<T>
{
public:
  using SuperClass = XMLContents<T>;
  int
  PrintSelf(std::ostream &, int indent) const override
  {
    // indent+=SuperClass::PrintSelf(os, indent);
    // os << this->PrintSpaces(indent) << "=== NumericValue ===" <<
    // this->m_Value << std::endl;
    return indent;
  }

  NumericValue(const std::string & name, T value)
    : XMLContents<T>(name)
    , m_Value(value)
  {}

  T
  GetValue(void) const override
  {
    return this->m_Value;
  }

  void
  SetValue(T val)
  {
    this->m_Value = val;
  }

protected:
  T m_Value;
};
#endif // NumericValue_h
