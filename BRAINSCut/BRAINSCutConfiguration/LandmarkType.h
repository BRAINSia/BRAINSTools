#if 0
#ifndef LandmarkType_h
#define LandmarkType_h
#include "ElementParser.h"
#include "StringValue.h"
class LandmarkType : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== LandmarkType ===" << std::endl;
    return indent + 2;
  }

  LandmarkType() : ElementParser("Landmark")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class LandmarkList : public ElementParser
{
public:
  LandmarkList() : ElementParser("LandmarkList")
  {
  }
};

#endif // xoLandmarkType_h
#endif
