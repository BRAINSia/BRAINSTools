#ifndef ImageDescription_h
#define ImageDescription_h

#include "StringValue.h"
#include "ElementParser.h"

class ImageDescription : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ImageDescription ===" << std::endl;
    return indent + 2;
  }

  ImageDescription() : ElementParser("Image")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class ImageList : public ElementParser
{
public:
  ImageList() : ElementParser("ImageList")
  {
  }
};

#endif // ImageDescription_h
