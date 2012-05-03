#ifndef SVMParams_h
#define SVMParams_h
#include "StringValue.h"
#include "IntValue.h"
#include "FloatValue.h"
#include "ElementParser.h"

class SVMParams : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== SVMParams ===" << std::endl;
    return indent + 2;
  }

  SVMParams() : ElementParser("AnnParams")
  {
    this->Add(new IntValue("MaximumVectorsPerEpoch",
                           2000), "MaximumVectorsPerEpoch");
    //    this->Add(new FloatValue("GaussianSize",5.5),"GaussianSize");
  }
};

#endif // SVMParams_h
