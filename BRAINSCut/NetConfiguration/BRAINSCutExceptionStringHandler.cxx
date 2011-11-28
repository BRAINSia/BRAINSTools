#include "BRAINSCutExceptionStringHandler.h"

BRAINSCutExceptionStringHandler
::BRAINSCutExceptionStringHandler(const std::string & errorString)
{
  this->m_ErrorString = "****ERROR**** [BRAINSCutExceptionStringHandler]:: ";
  this->m_ErrorString += errorString;
  this->m_ErrorString += "\n";
}

BRAINSCutExceptionStringHandler
::BRAINSCutExceptionStringHandler(const char *errorString)
{
  BRAINSCutExceptionStringHandler( std::string( errorString ) );
}

const std::string &
BRAINSCutExceptionStringHandler
::Error() const
{
  return m_ErrorString;
}

std::ostream & operator<<(std::ostream& stream, BRAINSCutExceptionStringHandler ob)
{
  stream << ob.Error();
  return stream;
}
