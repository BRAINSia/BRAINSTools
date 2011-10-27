#ifndef FileToRead_h
#define FileToRead_h
#include "FileSystemDescriptor.h"
#include <iostream>

template <typename OutputType>
class FileToRead :
  public FileSystemDescriptor<OutputType>
{
public:
  typedef FileSystemDescriptor<OutputType> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== FileToRead ===" << std::endl;
    return indent + 2;
  }

  FileToRead(const std::string & name, const std::string & filename) : FileSystemDescriptor<OutputType>(name,
                                                                                                        filename)
  {
  }

  FileToRead()
  {
  }

  virtual bool Verify()
  {
    bool returnvalue = true;

    if( this->m_Filename == "" )
      {
      std::cerr << "No filename specified." << std::endl;
      returnvalue = false;
      }
    if( !this->Exists() )
      {
      std::cerr << "File does not exists" << this->m_Filename << std::endl;
      returnvalue = false;
      }
    if( !this->IsReadable() )
      {
      std::cerr << "File is not readable " << this->m_Filename << std::endl;
      returnvalue = false;
      }
    return returnvalue;
  }
};

#endif // FileToRead_h
