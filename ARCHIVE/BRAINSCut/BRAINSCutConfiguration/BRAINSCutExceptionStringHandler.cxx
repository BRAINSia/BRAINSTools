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
#include "BRAINSCutExceptionStringHandler.h"

static std::string
LocalFormatErrorStringWrapper( const std::string & errorString )
{
  std::string buildErrorString( "****ERROR**** [BRAINSCutExceptionStringHandler]:: " );

  buildErrorString += errorString;
  buildErrorString += "\n";
  return buildErrorString;
}

BRAINSCutExceptionStringHandler ::BRAINSCutExceptionStringHandler( const std::string & errorString )
{
  this->m_ErrorString = LocalFormatErrorStringWrapper( errorString );
}

BRAINSCutExceptionStringHandler ::BRAINSCutExceptionStringHandler( const char * errorString )
{
  this->m_ErrorString = LocalFormatErrorStringWrapper( errorString );
}

const std::string &
BRAINSCutExceptionStringHandler ::Error() const
{
  return m_ErrorString;
}

std::ostream &
operator<<( std::ostream & stream, BRAINSCutExceptionStringHandler ob )
{
  stream << ob.Error();
  return stream;
}
