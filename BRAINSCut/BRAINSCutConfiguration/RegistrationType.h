#ifndef RegistrationType_h
#define RegistrationType_h
#include "StringValue.h"
#include "ElementParser.h"

class RegistrationType : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== RegistrationType ===" << std::endl;
    return indent + 2;
  }

  RegistrationType() : ElementParser("Registration")
  {
    this->Add(new StringValue("SubjToAtlasRegistrationFilename", ""),
              "SubjToAtlasRegistrationFilename");
    this->Add(new StringValue("AtlasToSubjRegistrationFilename", ""),
              "AtlasToSubjRegistrationFilename");
    this->Add(new StringValue("ID", ""),
              "ID");
  }
};

class RegistrationList : public ElementParser
{
public:
  RegistrationList() : ElementParser("RegistrationList")
  {
  }
};
#endif
