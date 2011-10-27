#ifndef BRAINSCutExceptionHandler_h
#define BRAINSCutExceptionHandler_h
#include <string>

class BRAINSCutExceptionStringHandler
{
public:
  BRAINSCutExceptionStringHandler(const std::string & errorString)
  {
    this->m_ErrorString = "****ERROR**** [BRAINSCutExceptionStringHandler]:: ";
    this->m_ErrorString += errorString;
    this->m_ErrorString += "\n";
  }

  BRAINSCutExceptionStringHandler(const char *errorString)
  {
    BRAINSCutExceptionStringHandler( std::string( errorString ) );
  }

  const std::string & Error() const
  {
    return m_ErrorString;
  }

private:
  std::string m_ErrorString;
};
#endif
