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
#ifndef __muException_h
#define __muException_h

#include <iostream>
#include <sstream>
#include <string>

#include "Log.h"

namespace mu
{
/** \class Exception
 */
class Exception : public std::exception
{
public:

  Exception()
  throw ( ) :
    m_Message("")
  {
  }

  ~Exception()
  throw ( )
  {
  }

  void SetMessage(const char *s)
  {
    m_Message = s;
  }

  void Print(std::ostream & os) const
  {
    os << m_Message << std::endl;
  }

  const char * what() const
  throw ( )
  {
    return m_Message.c_str();
  }

protected:

  std::string m_Message;
};
} // namespace mu

inline std::ostream & operator<<(std::ostream & os, mu::Exception & e)
{
  ( &e )->Print(os);
  return os;
}

#define muExceptionMacro(x)                                             \
    {                                                                     \
    muLogMacro( << "mu::Exception, in " << __FILE__ << " line " << __LINE__; \
                std::cerr << "\n" x << "\n");                           \
    std::stringstream oss;                                              \
    oss << "mu::Exception, in " << __FILE__ << " line " << __LINE__;    \
    oss << "\n" x << std::ends;                                         \
    mu::Exception e;                                                    \
    e.SetMessage( oss.str().c_str() );                                  \
    throw e;                                                            \
    }

#endif
