#include "Log.h"

#include "muException.h"

namespace mu
{
Log *
Log
::GetInstance()
{
  // Allow only one instance
  static Log instance;

  return &instance;
}

Log
::Log()
{
  m_EchoFlag = true;
  m_OutputFileName = "";
}

Log
::~Log()
{
  this->CloseFile();
}

Log
::Log(const Log & l)
{
  m_EchoFlag = l.m_EchoFlag;
  m_OutputFileName = l.m_OutputFileName;
}

void
Log
::CloseFile()
{
  if( m_Output.is_open() )
    {
    m_Output.close();
    }
}

void
Log
::SetOutputFileName(const char *s)
{
  if( m_Output.is_open() )
    {
    m_Output.close();
    }

  m_Output.open(s);

  if( m_Output.fail() )
    {
    muExceptionMacro(
      << "[Log::SetOutputFileName] Failed to open " << s);
    }
}

void
Log
::SetOutputFileName(const std::string & s)
{
  this->SetOutputFileName( s.c_str() );
}

void
Log
::WriteString(const char *s)
{
  if( s == NULL )
    {
    std::cout << "[Log::WriteString] NULL argument" << std::endl << std::flush;
    return;
    }

  if( m_Output.good() )
    {
    m_Output << s;
    m_Output.flush();
    }

  if( m_EchoFlag )
    {
    std::cout << s;
    (std::cout).flush();
    }
}

void
Log
::WriteString(const std::string & s)
{
  this->WriteString( s.c_str() );
}
} // namespace mu
