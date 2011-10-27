#ifndef ProbabilityMap_h
#define ProbabilityMap_h
#include "StringValue.h"
#include "FloatValue.h"
#include "CompoundObjectBase.h"

class ProbabilityMap : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ProbabilityMap ===" << std::endl;
    return indent + 2;
  }

  ProbabilityMap() : CompoundObjectBase("ProbabilityMap")
  {
    this->Add(new StringValue("StructureID", ""), "StructureID");
    this->Add(new StringValue("GenerateVector", ""), "GenerateVector");
    this->Add(new FloatValue("Gaussian", 1.0), "Gaussian");
    this->Add(new StringValue("Filename", ""), "Filename");
    this->Add(new StringValue("rho", ""), "rho");
    this->Add(new StringValue("phi", ""), "phi");
    this->Add(new StringValue("theta", ""), "theta");
  }
};

class ProbabilityMapList : public CompoundObjectBase
{
public:
  ProbabilityMapList() : CompoundObjectBase("ProbabilityMapList")
  {
  }
};

#endif
