#include "BRAINSCutExceptionStringHandler.h"

static std::string LocalFormatErrorStringWrapper(const std::string & errorString)
{
  std::string buildErrorString("****ERROR**** [BRAINSCutExceptionStringHandler]:: ");

  buildErrorString += errorString;
  buildErrorString += "\n";
  return buildErrorString
}

BRAINSCutExceptionStringHandler
::BRAINSCutExceptionStringHandler(const std::string & errorString)
{
  this->m_ErrorString = LocalFormatErrorStringWrapper(errorString);
}

BRAINSCutExceptionStringHandler
::BRAINSCutExceptionStringHandler(const char *errorString)
{
  this->m_ErrorString = LocalFormatErrorStringWrapper(errorString);
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
