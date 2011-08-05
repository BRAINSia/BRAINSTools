#ifndef ImageDescription_h
#define ImageDescription_h

#include "StringValue.h"
#include "CompoundObjectBase.h"

class ImageDescription : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ImageDescription ===" << std::endl;
    return indent + 2;
  }

  ImageDescription() : CompoundObjectBase("Image")
  {
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("Filename", ""), "Filename");
  }
};

class ImageList : public CompoundObjectBase
{
public:
  ImageList() : CompoundObjectBase("ImageList")
  {
  }
};

#endif // ImageDescription_h
