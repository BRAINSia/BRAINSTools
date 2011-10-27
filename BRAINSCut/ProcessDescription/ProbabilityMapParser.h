#ifndef ProbabilityMapParser_h
#define ProbabilityMapParser_h
#include "StringValue.h"
#include "FloatValue.h"
#include "XMLElementParser.h"

class ProbabilityMapParser : public XMLElementParser
{
public:
  typedef XMLElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ProbabilityMapParser ===" << std::endl;
    return indent + 2;
  }

  ProbabilityMapParser() : XMLElementParser("ProbabilityMapParser")
  {
    this->Add(new StringValue("StructureID", ""), "StructureID");
    this->Add(new StringValue("GenerateVector", ""), "GenerateVector");
    this->Add(new FloatValue("Gaussian", 1.0), "Gaussian");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class ProbabilityMapList : public XMLElementParser
{
public:
  ProbabilityMapList() : XMLElementParser("ProbabilityMapList")
  {
    std::cout << __LINE__ << __FILE__ << std::endl;
  }
};

#endif
