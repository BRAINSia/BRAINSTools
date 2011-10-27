#if 0
#ifndef LandmarkType_h
#define LandmarkType_h
#include "XMLElementParser.h"
#include "StringValue.h"
class LandmarkType : public XMLElementParser
{
public:
  typedef XMLElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== LandmarkType ===" << std::endl;
    return indent + 2;
  }

  LandmarkType() : XMLElementParser("Landmark")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class LandmarkList : public XMLElementParser
{
public:
  LandmarkList() : XMLElementParser("LandmarkList")
  {
  }
};

#endif // xoLandmarkType_h
#endif
