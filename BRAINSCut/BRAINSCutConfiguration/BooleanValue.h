#ifndef BooleanValue_h
#define BooleanValue_h
#include "ElementContainer.h"
#include <iostream>

class BooleanValue :
  public XMLContents<const bool>
{
public:
  typedef XMLContents<const bool> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== BooleanValue ===!"
       << this->m_Value << "!" << std::endl;
    return indent + 2;
  }

  typedef const bool ReturnType;
  BooleanValue(const std::string & name, ReturnType value) :
    XMLContents<ReturnType>(name),
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

  virtual bool Verify() const
  {
    return true;
  }

  void SetValue(const std::string & stringval)
  {
    std::string s = stringval;

    for( unsigned i = 0; i < s.size(); i++ )
      {
      s[i] = ::tolower(s[i]);
      }

    if( s != "true" &&
        s != "false" )
      {
      std::string msg("Can't convert *");
      msg += stringval;
      msg += ") to boolean";
      throw BRAINSCutExceptionStringHandler(msg);
      }
    bool returnVal = ( s == "true" );
    SetValue( returnVal );
  }

private:
  bool m_Value;
};

#endif // BooleanValue_h
