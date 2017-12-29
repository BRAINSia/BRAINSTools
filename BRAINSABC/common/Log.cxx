/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "Log.h"

#include "muException.h"

#include "itkMacro.h" //Needed for nullptr

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
::Log() :
  m_OutputFileName("")
{
  m_EchoFlag = true;
}

Log
::~Log()
{
  this->CloseFile();
}

Log
::Log(const Log & l) :
  m_EchoFlag( l.m_EchoFlag ),
  m_OutputFileName( l.m_OutputFileName )
{
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
    muExceptionMacro(<< "[Log::SetOutputFileName] Failed to open " << s);
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
  if( s == nullptr )
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
