#ifndef SpatialLocationType_h
#define SpatialLocationType_h
// #include "StringValue.h"
#include "ElementParser.h"

class SpatialLocationType : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== SpatialLocationType ===" << std::endl;
    return indent + 2;
  }

  SpatialLocationType() : ElementParser("SpatialLocation")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class SpatialLocationList : public ElementParser
{
public:
  SpatialLocationList() : ElementParser("SpatialLocationList")
  {
  }
};

#endif // SpatialLocationType
