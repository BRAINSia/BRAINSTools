#ifndef DoubleToString_h
#define DoubleToString_h

#include "double-conversion.h"
#include "itkMacro.h"

class DoubleToString
{
public:
  DoubleToString();
  std::string operator()(double val);

private:
  DoubleToString & operator=(const DoubleToString &); // not defined

  const double_conversion::DoubleToStringConverter & m_DoubleToStringConverter;
};

#endif // DoubleToString_h
