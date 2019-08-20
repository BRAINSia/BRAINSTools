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
#ifndef BRAINSCutExceptionHandler_h
#define BRAINSCutExceptionHandler_h
#include <string>
#include <ostream>

class BRAINSCutExceptionStringHandler
{
public:
  BRAINSCutExceptionStringHandler( const std::string & errorString );
  BRAINSCutExceptionStringHandler( const char * errorString );
  const std::string &
  Error() const;

  friend std::ostream &
  operator<<( std::ostream & stream, BRAINSCutExceptionStringHandler ob );

private:
  std::string m_ErrorString;
};

#endif
