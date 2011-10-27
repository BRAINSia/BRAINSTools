#ifndef IntValue_h
#define IntValue_h
#include "NumericValue.h"

class IntValue :
  public NumericValue<long>
{
public:
  typedef NumericValue<long> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent = SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== IntValue ===!" << this->m_Value
       << "!" << std::endl;
    return indent + 2;
  }

  typedef long ReturnType;
  IntValue(const std::string & name, ReturnType value) :
    NumericValue<ReturnType>(name, value)
  {
  }

  IntValue(const std::string & name, const std::string & stringval) :
    NumericValue<ReturnType>(name, 0)
  {
    this->SetValue(stringval);
  }

  void SetValue(const std::string & stringval);

  virtual bool Verify() const;

private:
};
#endif // IntValue_h
