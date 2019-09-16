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
#include "XMLParser.h"
#include <itksys/SystemTools.hxx>
#include <iostream>
#include <fstream>
#include <ElementContainer.h>
#include "BRAINSCutExceptionStringHandler.h"
static void
SE(void * parser, const XML_Char * name, const XML_Char ** atts)
{
  XMLParser * Parser = static_cast<XMLParser *>(parser);

  Parser->StartElement(Parser->GetUserData(), name, atts);
}

static void
EE(void * parser, const XML_Char * name)
{
  XMLParser * Parser = static_cast<XMLParser *>(parser);

  Parser->EndElement(Parser->GetUserData(), name);
}

bool
XMLParser::Parse()
{
  this->m_Parser = XML_ParserCreate(nullptr);

  XML_SetElementHandler(this->m_Parser, &SE, &EE);
  XML_SetUserData(this->m_Parser, this);

  std::ifstream inputstream;

  inputstream.open(this->m_Filename.c_str(), std::ios::binary | std::ios::in);
  if (inputstream.fail())
  {
    std::string message = "Can't open ";
    message += this->m_Filename;
    message += '\n';
    BRAINSCutExceptionStringHandler exception(message);
    throw exception;
  }

  // Default stream parser just reads a block at a time.
  std::streamsize filesize = itksys::SystemTools::FileLength(this->m_Filename.c_str());

  this->m_Buffer = new char[filesize];

  inputstream.read(this->m_Buffer, filesize);

  if (static_cast<std::streamsize>(inputstream.gcount()) != filesize)
  {
    BRAINSCutExceptionStringHandler exception("File Read Error");
    throw exception;
  }
  return XML_Parse(this->m_Parser, this->m_Buffer, inputstream.gcount(), true);
}
