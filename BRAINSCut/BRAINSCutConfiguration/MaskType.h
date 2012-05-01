#ifndef MaskType_h
#define MaskType_h
#include "StringValue.h"
#include "ElementParser.h"

class MaskType : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== MaskType ===" << std::endl;
    return indent + 2;
  }

  MaskType() : ElementParser("Mask")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class MaskList : public ElementParser
{
public:
  MaskList() : ElementParser("MaskList")
  {
  }
};

#endif // MaskType
