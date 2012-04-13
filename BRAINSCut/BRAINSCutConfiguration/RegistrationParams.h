#ifndef RegistrationParams_h
#define RegistrationParams_h
#include "CompoundObjectBase.h"
#include "StringValue.h"

class RegistrationParams : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== RegistrationParams ==="
       << std::endl;
    return indent + 2;
  }

  RegistrationParams() : CompoundObjectBase("RegistrationParams")
  {
    this->Add(new StringValue("ImageTypeToUse", ""), "ImageTypeToUse");
    this->Add(new StringValue("ID", ""), "ID");
  }
};

#endif // RegistrationParams_h
