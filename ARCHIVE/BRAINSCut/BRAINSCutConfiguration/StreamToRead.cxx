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
#include "StreamToRead.h"
#include <iostream>
#include <fstream>
#include "BRAINSCutExceptionStringHandler.h"
#include "itkMacro.h" //Needed for nullptr

int
StreamToRead ::PrintSelf( std::ostream & os, int indent ) const
{
  indent += Superclass::PrintSelf( os, indent );
  os << this->PrintSpaces( indent ) << "=== StreamToRead ===" << std::endl;
  return indent + 2;
}

StreamToRead ::StreamToRead( const std::string & name, const std::string & filename )
  : FileToRead< std::fstream * >( name, filename )
  , m_F( nullptr )
{
  if ( this->m_Filename == "" )
  {
    std::string msg( "Missing File Name For " );
    msg += this->GetName();
    throw BRAINSCutExceptionStringHandler( msg );
  }
  if ( this->m_F != nullptr )
  {
    if ( this->m_F->is_open() )
    {
      this->m_F->close();
    }
    delete this->m_F;
  }
  this->m_F = new FileStreamType( this->m_Filename.c_str(), std::ifstream::in );
  if ( m_F->fail() )
  {
    std::string msg( "Can't open " );
    msg += this->m_Filename;
    throw BRAINSCutExceptionStringHandler( msg );
  }
}

StreamToRead::OutputType
StreamToRead ::GetValue() const
{
  // this->m_F.clear();              // forget we hit the end of file
  // this->m_F.seekg(0, ios::beg);   // move to the start of the file
  return this->m_F;
}

StreamToRead::~StreamToRead()
{
  this->Close();
}

void
StreamToRead::Close()
{
  if ( m_F != nullptr )
  {
    m_F->close();
    delete m_F;
  }
  m_F = nullptr;
}
