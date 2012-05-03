#ifndef ProbabilityMapParser_h
#define ProbabilityMapParser_h
// #include "StringValue.h"
// #include "FloatValue.h"
#include "ElementParser.h"

class ProbabilityMapParser : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ProbabilityMapParser ===" << std::endl;
    return indent + 2;
  }

  ProbabilityMapParser() : ElementParser("ProbabilityMapParser")
  {
    this->Add(new StringValue("StructureID", ""), "StructureID");
    this->Add(new StringValue("GenerateVector", ""), "GenerateVector");
    this->Add(new FloatValue("Gaussian", 1.0), "Gaussian");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class ProbabilityMapList : public ElementParser
{
public:
  ProbabilityMapList() : ElementParser("ProbabilityMapList")
  {
  }
};

#endif
