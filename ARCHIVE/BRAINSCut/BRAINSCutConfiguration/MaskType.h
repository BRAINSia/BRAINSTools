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
#ifndef MaskType_h
#define MaskType_h
#include "StringValue.h"
#include "ElementParser.h"

class MaskType : public ElementParser
{
public:
  using SuperClass = ElementParser;
  int
  PrintSelf(std::ostream & os, int indent) const override
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== MaskType ===" << std::endl;
    return indent + 2;
  }

  MaskType()
    : ElementParser("Mask")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class MaskList : public ElementParser
{
public:
  MaskList()
    : ElementParser("MaskList")
  {}
};

#endif // MaskType
