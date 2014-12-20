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
#ifndef IntValue_h
#define IntValue_h
#include "NumericValue.h"
#include "itkMacro.h" //Needed for ITK_OVERRIDE

class IntValue :
  public NumericValue<long>
{
public:
  typedef NumericValue<long> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const ITK_OVERRIDE
  {
    indent = SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== IntValue ===!" << this->m_Value
       << "!" << std::endl;
    return indent + 2;
  }

  typedef long ReturnType;
  IntValue(const std::string & name, ReturnType value) :
    NumericValue<ReturnType>(name, value)
  {
  }

  IntValue(const std::string & name, const std::string & stringval) :
    NumericValue<ReturnType>(name, 0)
  {
    this->SetValue(stringval);
  }

  void SetValue(const std::string & stringval);

  virtual bool Verify() const ITK_OVERRIDE;

private:
};
#endif // IntValue_h
