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
#ifndef StringValue_h
#define StringValue_h
#include "ElementContainer.h"
#include <string>
#include <iostream>
#include "itkMacro.h" //Needed for override

class StringValue :
  public XMLContents<const std::string>
{
public:
  typedef XMLContents<const std::string> SuperClass;
  int PrintSelf(std::ostream & os, int indent) const override
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== StringValue ===!"
       << this->m_Value << "!" << std::endl;
    return indent + 2;
  }

  typedef const std::string ReturnType;
  StringValue(ReturnType & name, ReturnType & value) :
    XMLContents<ReturnType>(name),
    m_Value(value)
  {
  }

  ReturnType GetValue(void) const override
  {
    return this->m_Value;
  }

  void SetValue(ReturnType & s)
  {
    this->m_Value = s;
  }

  //
  // presumably, if you cared if the string wasn't empty
  // or you had particular values in mind, you'd derive from
  // StringValue;
  bool Verify() const override
  {
    bool returnvalue = true;

    if( this->m_Value == "" )
      {
      std::cerr << "Empty string value for " << this->GetName() << std::endl;
      returnvalue = false;
      }
    return returnvalue;
  }

private:
  std::string m_Value;
};

#endif // StringValue_h
