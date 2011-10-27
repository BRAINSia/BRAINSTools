#ifndef MaskType_h
#define MaskType_h
#include "StringValue.h"
#include "XMLElementParser.h"

class MaskType : public XMLElementParser
{
public:
  typedef XMLElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== MaskType ===" << std::endl;
    return indent + 2;
  }

  MaskType() : XMLElementParser("Mask")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class MaskList : public XMLElementParser
{
public:
  MaskList() : XMLElementParser("MaskList")
  {
  }
};

#endif // MaskType
