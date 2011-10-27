#include "StreamToRead.h"
#include <iostream>
#include <fstream>
#include "BRAINSCutExceptionStringHandler.h"

StreamToRead::ReturnType
StreamToRead::GetValue(void)
{
  if( this->m_Filename == "" )
    {
    std::string msg("Missing File Name For ");
    msg += this->GetName();
    throw BRAINSCutExceptionStringHandler(msg);
    }
  if( this->m_F != 0 )
    {
    if( this->m_F->is_open() )
      {
      this->m_F->close();
      }
    delete this->m_F;
    }
  this->m_F = new StreamType(this->m_Filename.c_str(), std::ifstream::in);
  if( m_F->fail() )
    {
    std::string msg("Can't open ");
    msg += this->m_Filename;
    throw BRAINSCutExceptionStringHandler(msg);
    }
  return this->m_F;
}

StreamToRead::
~StreamToRead()
{
  this->Close();
}

void
StreamToRead::Close()
{
  if( m_F != 0 )
    {
    m_F->close();
    delete m_F;
    }
  m_F = 0;
}
