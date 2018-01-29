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
#ifndef FileSystemDescriptor_H
#define FileSystemDescriptor_H
#include "ElementContainer.h"
#include <itksys/SystemTools.hxx>
#include <fstream>
template <typename TOutputType>
class FileSystemDescriptor :
  public XMLContents<TOutputType>
{
public:
  typedef XMLContents<TOutputType> SuperClass;
  typedef TOutputType              OutputType;
  int PrintSelf(std::ostream & os, int indent) const override
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== FileSystemDescriptor ==="
       << this->m_Filename << std::endl;
    return indent + 2;
  }

  FileSystemDescriptor(const std::string & name, const std::string & filename) :
    XMLContents<TOutputType>(name)
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

  bool Exists() const
  {
    return itksys::SystemTools::FileExists( m_Filename.c_str() );
  }

  bool IsReadable() const
  {
    std::fstream f(this->m_Filename.c_str(), std::fstream::in);
    const bool   rval = f.is_open();

    if( rval == true )
      {
      f.close();
      }
    return rval;
  }

  bool IsDirectory() const
  {
    return itksys::SystemTools::FileIsDirectory( m_Filename.c_str() );
  }

  virtual void Close() = 0;

  bool Verify() constexpr override  = 0;

  OutputType GetValue() constexpr override  = 0;

protected:
  std::string m_Filename;
};

#endif // FileSystemDescriptor_H
