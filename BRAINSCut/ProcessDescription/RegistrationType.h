#ifndef RegistrationType_h
#define RegistrationType_h
#include "StringValue.h"
#include "CompoundObjectBase.h"

class RegistrationType : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== RegistrationType ===" << std::endl;
    return indent + 2;
  }

  RegistrationType() : CompoundObjectBase("Registration")
  {
    this->Add(new StringValue("SubjToAtlasRegistrationFilename", ""),
              "SubjToAtlasRegistrationFilename");
    this->Add(new StringValue("AtlasToSubjRegistrationFilename", ""),
              "AtlasToSubjRegistrationFilename");
    this->Add(new StringValue("AtlasBinaryFilename", "NA"),
              "AtlasBinaryFilename");
    this->Add(new StringValue("SubjectBinaryFilename", "NA"),
              "SubjectBinaryFilename");
    //    this->Add(new StringValue("LandmarkType",""),
    //             "LandmarkType");
    this->Add(new StringValue("ID", ""),
              "ID");
  }
};

class RegistrationList : public CompoundObjectBase
{
public:
  RegistrationList() : CompoundObjectBase("RegistrationList")
  {
  }
};
#endif
