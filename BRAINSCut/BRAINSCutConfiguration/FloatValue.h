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
#ifndef FloatValue_h
#define FloatValue_h
#include "NumericValue.h"

class FloatValue :
  public NumericValue<double>
{
public:
  typedef NumericValue<double> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== FloatValue === !"
       <<  this->m_Value << "!" << std::endl;
    return indent + 2;
  }

  typedef double ReturnType;
  FloatValue(const std::string & name, ReturnType value) :
    NumericValue<ReturnType>(name, value)
  {
  }

  FloatValue(const std::string & name, const std::string & stringval) :
    NumericValue<ReturnType>(name, 0.0)
  {
    this->SetValue(stringval);
  }

  void SetValue(const std::string & stringval);

  virtual bool Verify() const;

private:
};
#endif // FloatValue_h
