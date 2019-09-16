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
#ifndef XMLParser_h
#define XMLParser_h
#include <string>
#include <expat.h>
#include "itkMacro.h"

class XMLParser
{
public:
  XMLParser(const std::string & filename)
    : m_Filename(filename)
    , m_UserData(nullptr)
    , m_Buffer(nullptr)
  {}

  void
  SetUserData(void * userData)
  {
    m_UserData = userData;
  }

  void *
  GetUserData()
  {
    return m_UserData;
  }

  bool
  Parse();

  virtual void
  StartElement(void * userData, const XML_Char * name, const XML_Char ** atts) = 0;

  virtual void
  EndElement(void * userData, const XML_Char * name) = 0;

  virtual ~XMLParser()
  {
    XML_ParserFree(this->m_Parser);
    delete[] this->m_Buffer;
  }

private:
  std::string m_Filename;
  void *      m_UserData;
  char *      m_Buffer;
  XML_Parser  m_Parser;
};
#endif // XMLParser_h
