#ifndef RegistrationConfigurationParser_h
#define RegistrationConfigurationParser_h
#include "XMLElementParser.h"
// #include "StringValue.h"
// #include "IntValue.h"

class RegistrationConfigurationParser : public XMLElementParser
{
public:
  typedef XMLElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== RegistrationConfiguration ==="
       << std::endl;
    return indent + 2;
  }

  RegistrationConfigurationParser() : XMLElementParser("RegistrationConfiguration")
  {
    this->Add(new StringValue("ImageTypeToUse", ""), "ImageTypeToUse");
    this->Add(new StringValue("ID", ""), "ID");
    this->Add(new IntValue("BRAINSROIAutoDilateSize", 1), "BRAINSROIAutoDilateSize");
  }
};

#endif // RegistrationConfigurationParser_h
