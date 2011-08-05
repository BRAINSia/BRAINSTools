#if 0
#ifndef LandmarkType_h
#define LandmarkType_h
#include "CompoundObjectBase.h"
#include "StringValue.h"
class LandmarkType : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== LandmarkType ===" << std::endl;
    return indent + 2;
  }

  LandmarkType() : CompoundObjectBase("Landmark")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class LandmarkList : public CompoundObjectBase
{
public:
  LandmarkList() : CompoundObjectBase("LandmarkList")
  {
  }
};

#endif // xoLandmarkType_h
#endif
