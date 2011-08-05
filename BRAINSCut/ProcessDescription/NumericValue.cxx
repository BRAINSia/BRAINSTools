#include "FloatValue.h"
#include "IntValue.h"
#include <stdlib.h>

//
// for each type define Verify, and SetValue (from string)

// eventually set up domain/range checking?
bool
IntValue::Verify() const
{
  return true;
}

bool
FloatValue::Verify() const
{
  return true;
}

void
IntValue::SetValue(const std::string & stringval)
{
  char *test;
  long  val = strtol(stringval.c_str(), &test, 10);

  if( test == stringval.c_str() )
    {
    std::string msg("Can't convert *");
    msg += stringval;
    msg += ") to integer";
    throw ProcessObjectException(msg);
    }
  this->m_Value = val;
}

void
FloatValue::SetValue(const std::string & stringval)
{
  char * test;
  double val = strtod(stringval.c_str(), &test);

  if( test == stringval.c_str() )
    {
    std::string msg("Can't convert *");
    msg += stringval;
    msg += ") to float";
    throw ProcessObjectException(msg);
    }
  this->NumericValue<double>::SetValue(val);
}
