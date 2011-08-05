#ifndef MaskType_h
#define MaskType_h
#include "StringValue.h"
#include "CompoundObjectBase.h"

class MaskType : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== MaskType ===" << std::endl;
    return indent + 2;
  }

  MaskType() : CompoundObjectBase("Mask")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class MaskList : public CompoundObjectBase
{
public:
  MaskList() : CompoundObjectBase("MaskList")
  {
  }
};

#endif // MaskType
