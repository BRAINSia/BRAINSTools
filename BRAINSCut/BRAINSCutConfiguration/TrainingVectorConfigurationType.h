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
#ifndef TrainingVectorConfigurationType_h
#define TrainingVectorConfigurationType_h
#include "StringValue.h"
#include "IntValue.h"
#include "FloatValue.h"
#include "ElementParser.h"

class TrainingVectorConfigurationType : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== TrainingVectorConfigurationType ===" << std::endl;
    return indent + 2;
  }

  TrainingVectorConfigurationType() : ElementParser("TrainingVectorConfiguration")
  {
    this->Add(new FloatValue("MaskSmoothingValue", 0.0), "MaskSmoothingValue");
    this->Add(new IntValue("GradientProfileSize", 1), "GradientProfileSize");
    this->Add(new StringValue("TrainingVectorFilename", ""), "TrainingVectorFilename");
    this->Add(new StringValue("TestVectorFilename", ""), "TestVectorFilename");
    this->Add(new StringValue("TrainingModelFilename", ""), "TrainingModelFilename");
    this->Add(new StringValue("Normalization", ""), "Normalization");
  }
};

#endif // TrainingVectorConfigurationType_h
