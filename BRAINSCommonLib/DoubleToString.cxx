#include "DoubleToString.h"

DoubleToString
::DoubleToString() :
  m_DoubleToStringConverter(double_conversion::DoubleToStringConverter::EcmaScriptConverter() )
{
}

std::string
DoubleToString
::operator()(double val)
{
  char                             buf[256];
  double_conversion::StringBuilder builder(buf, sizeof(buf) );

  builder.Reset();
  if( !m_DoubleToStringConverter.ToShortest(val, &builder) )
    {
    itkGenericExceptionMacro(<< "Conversion failed for " << val);
    }
  return std::string(builder.Finalize() );
}
