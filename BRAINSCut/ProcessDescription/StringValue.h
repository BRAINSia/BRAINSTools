#ifndef StringValue_h
#define StringValue_h
#include "ProcessObjectBase.h"
#include <string>
#include <iostream>

class StringValue :
  public ValueObjectBase<const std::string>
{
public:
  typedef ValueObjectBase<const std::string> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== StringValue ===!"
       << this->m_Value << "!" << std::endl;
    return indent + 2;
  }

  typedef const std::string ReturnType;
  StringValue(ReturnType & name, ReturnType & value) :
    ValueObjectBase<ReturnType>(name),
    m_Value(value)
  {
  }

  ReturnType GetValue(void) const
  {
    return this->m_Value;
  }

  void SetValue(ReturnType & s)
  {
    this->m_Value = s;
  }

  //
  // presumably, if you cared if the string wasn't empty
  // or you had particular values in mind, you'd derive from
  // StringValue;
  virtual bool Verify() const
  {
    bool returnvalue = true;

    if( this->m_Value == "" )
      {
      std::cerr << "Empty string value for " << this->GetName() << std::endl;
      returnvalue = false;
      }
    return returnvalue;
  }

private:
  std::string m_Value;
};

#endif // StringValue_h
