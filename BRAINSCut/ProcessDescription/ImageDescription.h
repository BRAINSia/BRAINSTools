#ifndef ImageDescription_h
#define ImageDescription_h

#include "StringValue.h"
#include "XMLElementParser.h"

class ImageDescription : public XMLElementParser
{
public:
  typedef XMLElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ImageDescription ===" << std::endl;
    return indent + 2;
  }

  ImageDescription() : XMLElementParser("Image")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class ImageList : public XMLElementParser
{
public:
  ImageList() : XMLElementParser("ImageList")
  {
  }
};

#endif // ImageDescription_h
