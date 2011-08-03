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
  throw ( )
  {
    m_Message = "";
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
