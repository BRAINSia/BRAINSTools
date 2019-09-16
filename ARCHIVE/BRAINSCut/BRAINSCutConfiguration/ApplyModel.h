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
#ifndef ApplyModel_h
#define ApplyModel_h
#include "StringValue.h"
#include "IntValue.h"
#include "FloatValue.h"
#include "ElementParser.h"

class ApplyModelType : public ElementParser
{
public:
  using SuperClass = ElementParser;
  int
  PrintSelf(std::ostream & os, int indent) const override
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ApplyModelType ===" << std::endl;
    return indent + 2;
  }

  ApplyModelType()
    : ElementParser("ApplyModel")
  {
    this->Add(new FloatValue("MaskThresh", 0.5), "MaskThresh");
    this->Add(new FloatValue("GaussianSmoothingSigma", 0.5), "GaussianSmoothingSigma");
  }
};

#endif // ApplyModel_h
