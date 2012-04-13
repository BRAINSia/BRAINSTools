#ifndef SpatialLocationType_h
#define SpatialLocationType_h
// #include "StringValue.h"
#include "XMLElementParser.h"

class SpatialLocationType : public XMLElementParser
{
public:
  typedef XMLElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== SpatialLocationType ===" << std::endl;
    return indent + 2;
  }

  SpatialLocationType() : XMLElementParser("SpatialLocation")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class SpatialLocationList : public XMLElementParser
{
public:
  SpatialLocationList() : XMLElementParser("SpatialLocationList")
  {
  }
};

#endif // SpatialLocationType
