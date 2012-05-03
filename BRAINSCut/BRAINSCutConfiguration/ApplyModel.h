#ifndef ApplyModel_h
#define ApplyModel_h
#include "StringValue.h"
#include "IntValue.h"
#include "FloatValue.h"
#include "ElementParser.h"

class ApplyModelType : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ApplyModelType ===" << std::endl;
    return indent + 2;
  }

  ApplyModelType() : ElementParser("ApplyModel")
  {
    this->Add(new FloatValue("MaskThresh", 0.5), "MaskThresh");
    this->Add(new FloatValue("GaussianSmoothingSigma", 0.5), "GaussianSmoothingSigma");
  }
};

#endif // ApplyModel_h
