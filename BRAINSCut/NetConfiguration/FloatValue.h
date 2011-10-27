#ifndef FloatValue_h
#define FloatValue_h
#include "NumericValue.h"

class FloatValue :
  public NumericValue<double>
{
public:
  typedef NumericValue<double> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== FloatValue === !"
       <<  this->m_Value << "!" << std::endl;
    return indent + 2;
  }

  typedef double ReturnType;
  FloatValue(const std::string & name, ReturnType value) :
    NumericValue<ReturnType>(name, value)
  {
  }

  FloatValue(const std::string & name, const std::string & stringval) :
    NumericValue<ReturnType>(name, 0.0)
  {
    this->SetValue(stringval);
  }

  void SetValue(const std::string & stringval);

  virtual bool Verify() const;

private:
};
#endif // FloatValue_h
