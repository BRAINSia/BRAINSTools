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
#ifndef XMLConfigurationFileParser_h
#define XMLConfigurationFileParser_h

#include <BRAINSCutConfiguration.h>
#include "XMLParser.h"
#include <list>
#include <vector>
#include <iostream>

class StringMap :
  public std::map<std::string, std::string>
{
public:
  typedef std::map<std::string, std::string> Superclass;
  typedef Superclass::const_iterator         const_iterator;
  std::string Get(const char *eleType,
                  const char *key)
  {
    StringMap::const_iterator it = this->find( std::string(key) );

    if( it == this->end() )
      {
      std::string msg("Missing Name attribute ");
      msg +=  key;
      msg += " in ";
      msg += eleType;
      throw BRAINSCutExceptionStringHandler(msg);
      }
    return it->second;
  }

  std::string GetIfExist(const char *eleType,
                         const char *key)
  {
    StringMap::const_iterator it = this->find( std::string(key) );

    if( it == this->end() )
      {
      std::cout << " *** Note ***" << std::endl
                << " *** Note ***" << std::endl
                << " Missing Name of attribute, mark as NA "
                << eleType << " " << key
                << std::endl << std::endl;
      return "NA";
      }
    return it->second;
  }
};

class XMLConfigurationFileParser : public XMLParser
{
public:
  XMLConfigurationFileParser(const std::string & filename ) : XMLParser(filename)
  {
    myConfiguration = new BRAINSCutConfiguration();
  }

  virtual void StartElement(void *userData, const XML_Char *name, const XML_Char * *atts) ITK_OVERRIDE;

  virtual void EndElement(void *userData, const XML_Char *name) ITK_OVERRIDE;

  BRAINSCutConfiguration * GetConfiguration();

  void ValidateDataSets();

private:
  // void ReadXML();

  BRAINSCutConfiguration * myConfiguration;
};

#endif // XMLConfigurationFileParser_h
