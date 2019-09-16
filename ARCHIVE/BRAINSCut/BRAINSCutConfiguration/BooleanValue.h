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
#ifndef BooleanValue_h
#define BooleanValue_h
#include "ElementContainer.h"
#include <iostream>
#include "itkMacro.h" //Needed for override

class BooleanValue : public XMLContents<bool>
{
public:
  using SuperClass = XMLContents<bool>;
  int
  PrintSelf(std::ostream & os, int indent) const override
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== BooleanValue ===!" << this->m_Value << "!" << std::endl;
    return indent + 2;
  }

  using ReturnType = bool;
  BooleanValue(const std::string & name, const ReturnType value)
    : XMLContents<ReturnType>(name)
    , m_Value(value)
  {}

  ReturnType
  GetValue(void) const override
  {
    return this->m_Value;
  }

  void
  SetValue(const ReturnType & s)
  {
    this->m_Value = s;
  }

  bool
  Verify() const override
  {
    return true;
  }

  void
  SetValue(const std::string & stringval)
  {
    std::string s = stringval;

    for (unsigned i = 0; i < s.size(); i++)
    {
      s[i] = ::tolower(s[i]);
    }

    if (s != "true" && s != "false")
    {
      std::string msg("Can't convert *");
      msg += stringval;
      msg += ") to boolean";
      throw BRAINSCutExceptionStringHandler(msg);
    }
    bool returnVal = (s == "true");
    SetValue(returnVal);
  }

private:
  bool m_Value;
};

#endif // BooleanValue_h
