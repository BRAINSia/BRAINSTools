#ifndef NumericValue_h
#define NumericValue_h
#include "ProcessObjectBase.h"

template <typename T>
class NumericValue :
  public ValueObjectBase<T>
{
public:
  typedef ValueObjectBase<T> SuperClass;
  virtual int PrintSelf(std::ostream &, int indent) const
  {
    // indent+=SuperClass::PrintSelf(os, indent);
    // os << this->PrintSpaces(indent) << "=== NumericValue ===" <<
    // this->m_Value << std::endl;
    return indent;
  }

  NumericValue(const std::string & name, T value) :
    ValueObjectBase<T>(name),
    m_Value(value)
  {
  }

  T GetValue(void) const
  {
    return this->m_Value;
  }

  void SetValue(T val)
  {
    this->m_Value = val;
  }

protected:
  T m_Value;
};
#endif // NumericValue_h
