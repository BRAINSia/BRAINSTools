#ifndef StreamToRead_h
#define StreamToRead_h
#include "FileToRead.h"
#include <iostream>
#include <fstream>

class StreamToRead :
  public FileToRead<std::fstream *>
{
public:
  typedef FileToRead<std::fstream *> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== StreamToRead ===" << std::endl;
    return indent + 2;
  }

  typedef std::fstream StreamType;
  typedef StreamType * ReturnType;
  StreamToRead(const std::string & name, const std::string & filename) : FileToRead<std::fstream *>(name, filename),
    m_F(0)
  {
  }

  StreamToRead() : m_F(0)
  {
  }

  virtual ~StreamToRead();
  virtual ReturnType GetValue(void);

  virtual void Close();

private:
  ReturnType m_F;
};
#endif // StreamToRead_h
