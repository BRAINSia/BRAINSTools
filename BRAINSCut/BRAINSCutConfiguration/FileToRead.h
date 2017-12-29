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
  int PrintSelf(std::ostream & os, int indent) const override
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

  bool Verify() const override
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
