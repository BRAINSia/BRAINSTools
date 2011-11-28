#ifndef StreamToRead_h
#define StreamToRead_h
#include "FileToRead.h"
#include <iostream>
#include <fstream>

typedef std::fstream FileStreamType;

class StreamToRead :
  public FileToRead<FileStreamType *>
{
public:
  typedef StreamToRead               Self;
  typedef FileToRead<std::fstream *> Superclass;
  typedef Superclass::OutputType     OutputType;

  virtual OutputType GetValue() const;

  virtual void Close();

  virtual int PrintSelf(std::ostream & os, int indent) const;

  StreamToRead(const std::string & name, const std::string & filename);
  virtual ~StreamToRead();
protected:
  StreamToRead();         // purposefully not implemented : m_F(0) { }
  void operator=(Self &); // purposefully not implemented

private:
  OutputType m_F;
};
#endif // StreamToRead_h
