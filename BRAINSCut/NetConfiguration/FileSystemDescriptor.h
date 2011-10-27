#ifndef FileSystemDescriptor_H
#define FileSystemDescriptor_H
#include "XMLElementContainer.h"
#include <itksys/SystemTools.hxx>
#include <fstream>
template <typename OutputType>
class FileSystemDescriptor :
  public XMLContents<OutputType>
{
public:
  typedef XMLContents<OutputType> SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== FileSystemDescriptor ==="
       << this->m_Filename << std::endl;
    return indent + 2;
  }

  FileSystemDescriptor(const std::string & name, const std::string & filename) :
    XMLContents<OutputType>(name)
  {
    this->SetFileName(filename);
  }

  FileSystemDescriptor()
  {
  }

  const std::string & GetFilename() const
  {
    return m_Filename;
  }

  void SetFileName(const std::string & s)
  {
    m_Filename = s;
  }

  bool Exists()
  {
    return itksys::SystemTools::FileExists( m_Filename.c_str() );
  }

  bool IsReadable()
  {
    std::fstream f(this->m_Filename.c_str(), std::fstream::in);
    bool         rval = f.is_open();

    if( rval == true )
      {
      f.close();
      }
    return rval;
  }

  bool IsDirectory()
  {
    return itksys::SystemTools::FileIsDirectory( m_Filename.c_str() );
  }

  virtual void Close() = 0;

protected:
  std::string m_Filename;
};

#endif // FileSystemDescriptor_H
