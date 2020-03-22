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
#include "AtlasDefinition.h"
#include <expat.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <itksys/SystemTools.hxx>
#include "itkNumberToString.h"

namespace // anon
{
extern "C"
{
  namespace AtlasXMLParser
  {
  void
  XMLstart(void * data, const char * el, const char **)
  {
    // std::cerr << "Start, El = ("
    //         << el
    //         << ")";
    // std::cerr << " Attributes = (";
    // for(unsigned i = 0; attr[i] != 0; i+= 2)
    //   {
    //   if(i > 0) std::cerr << ",";
    //   std::cerr << attr[i] << "=" << attr[i+1];
    //   }
    // std::cerr << ")" << std::endl;
    auto * _this = reinterpret_cast<AtlasDefinition *>(data);

    _this->XMLStart(el);
  }

  void
  XMLend(void * data, const char * el)
  {
    // std::cerr << "End, El = ("
    //         << el
    //         << ")" << std::endl;
    auto * _this = reinterpret_cast<AtlasDefinition *>(data);

    _this->XMLEnd(el);
  }

  void
  XMLcharhandler(void * data, const char * txt, int txtlen)
  {
    // std::cerr << "Char data = (";
    // for(unsigned i = 0; i < txtlen; i++)
    //   {
    //   std::cerr << txt[i];
    //   }
    // std::cerr << ")" << std::endl;
    auto buf = new char[txtlen + 1];
    int  i = 0;

    for (i = 0; i < txtlen; i++)
    {
      buf[i] = txt[i];
    }
    buf[i] = '\0';
    auto * _this = reinterpret_cast<AtlasDefinition *>(data);
    _this->XMLChar(buf);
    delete[] buf;
  }
  } // namespace AtlasXMLParser
}
} // namespace

AtlasDefinition::AtlasDefinition()

{
  this->m_TissueTypes.resize(0);
}

const AtlasDefinition::TissueTypeVector &
AtlasDefinition::TissueTypes()
{
  if (this->m_TissueTypes.empty())
  {
    for (PriorMapType::const_iterator it = this->m_PriorMap.begin(); it != this->m_PriorMap.end(); ++it)
    {
      this->m_TissueTypes.push_back(it->first);
    }
  }
  return this->m_TissueTypes;
}

void
AtlasDefinition::XMLStart(const char * el)
{
  std::string El(el);

  this->m_XMLElementStack.push_back(El);
  if (El == "Prior")
  {
    this->m_LastPriorBounds.clear();
  }
}

double
AtlasDefinition ::StrToD(const char * str, const char * message) const
{
  char * last = nullptr;
  double rval = strtod(str, &last);

  if (str == static_cast<const char *>(last))
  {
    std::cerr << message << ' ' << this->m_LastXMLString << std::endl;
    throw;
  }
  return rval;
}

long
AtlasDefinition ::StrToL(const char * str, const char * message) const
{
  char * last = nullptr;
  long   rval = strtol(str, &last, 10);

  if (str == static_cast<const char *>(last))
  {
    std::cerr << message << ' ' << this->m_LastXMLString << std::endl;
    throw;
  }
  return rval;
}

void
AtlasDefinition::XMLEnd(const char * el)
{
  std::string  El(el);
  const char * start = this->m_LastXMLString.c_str();

  // pop the current element name off the stack.
  this->m_XMLElementStack.pop_back();

  if (El == "Atlas")
  {
    return; // final element
  }
  const std::string ContainingElement = this->m_XMLElementStack.back();
  if (El == "filename")
  {
    this->m_LastFilename = this->m_LastXMLString;
  }
  else if (El == "type")
  {
    if (ContainingElement == "Prior")
    {
      this->m_LastPriorType = this->m_LastXMLString;
    }
    else
    {
      this->m_LastType = this->m_LastXMLString;
    }
  }
  else if (El == "AtlasImage")
  {
    this->m_TemplateVolumes[this->m_LastType] = this->m_LastFilename;
  }
  else if (El == "Weight")
  {
    this->m_LastWeight = this->StrToD(start, "Bad Weight given");
  }
  else if (El == "lower")
  {
    this->m_LastLower = this->StrToD(start, "Bad lower bound");
  }
  else if (El == "upper")
  {
    this->m_LastUpper = this->StrToD(start, "Bad upper bound");
  }
  else if (El == "GaussianClusterCount")
  {
    this->m_LastGaussianClusterCount = this->StrToL(start, "Bad GaussianClusterCount");
  }
  else if (El == "LabelCode")
  {
    this->m_LastLabelCode = this->StrToL(start, "Bad LabelCode");
  }
  else if (El == "UseForBias")
  {
    this->m_LastUseForBias = (this->StrToL(start, "Bad UseForBias") != 0);
  }
  else if (El == "IsForegroundPrior")
  {
    this->m_LastIsForegroundPrior = (this->StrToL(start, "Bad IsForegroundPrior") != 0);
  }
  else if (El == "BrainMask")
  {
    this->m_TemplateBrainMask = this->m_LastFilename;
  }
  else if (El == "HeadRegion")
  {
    this->m_TemplateHeadRegion = this->m_LastFilename;
  }
  else if (El == "Prior")
  {
    Prior curPrior;
    curPrior.SetFilename(this->m_LastFilename);
    curPrior.SetWeight(this->m_LastWeight);
    curPrior.SetGaussianClusterCount(this->m_LastGaussianClusterCount);
    curPrior.SetLabelCode(this->m_LastLabelCode);
    curPrior.SetUseForBias(this->m_LastUseForBias);
    curPrior.SetIsForegroundPrior(this->m_LastIsForegroundPrior);
    curPrior.SetBoundsList(this->m_LastPriorBounds);
    this->m_PriorMap[this->m_LastPriorType] = curPrior;
  }
  else if (El == "bounds")
  {
    BoundsType & curBounds = this->m_LastPriorBounds[this->m_LastType];
    curBounds.SetLower(this->m_LastLower);
    curBounds.SetUpper(this->m_LastUpper);
  }
  else
  {
    std::cerr << "Unhandled XML Element type " << El << std::endl;
    throw;
  }
}

void
AtlasDefinition::XMLChar(const char * buf)
{
  this->m_LastXMLString = buf;
}

void
AtlasDefinition::InitFromXML(const std::string & XMLFilename)
{
  std::ifstream xmlFile(XMLFilename.c_str(), std::ifstream::in);

  if (!xmlFile.good())
  {
    std::cout << "ERROR:  XML file " << XMLFilename << " can not be read properly " << std::flush << std::endl;
    throw;
  }
  std::streamsize fSize = itksys::SystemTools::FileLength(XMLFilename.c_str());

  XML_Parser parser = XML_ParserCreate(nullptr);
  XML_SetUserData(parser, static_cast<void *>(this));
  XML_SetElementHandler(parser, AtlasXMLParser::XMLstart, AtlasXMLParser::XMLend);
  XML_SetCharacterDataHandler(parser, AtlasXMLParser::XMLcharhandler);

  auto filebuf = new char[fSize];
  if (filebuf == nullptr)
  {
    std::cout << "ERROR:  memory char[" << fSize << "] can not be allocated properly " << std::flush << std::endl;
    throw;
  }

  xmlFile.read(filebuf, fSize);
  if (static_cast<std::streamsize>(xmlFile.gcount()) != fSize)
  {
    std::cout << "ERROR:  file not read proplerly " << XMLFilename << std::flush << std::endl;
    delete[] filebuf;
    throw;
  }
  xmlFile.close();

  int parserReturn(1);
  try
  {
    parserReturn = XML_Parse(parser, filebuf, fSize, 1);
  }
  catch (...)
  {
    parserReturn = 0;
  }
  if (parserReturn == 0)
  {
    delete[] filebuf;
    std::cerr << "XML File parsing error" << std::endl;
    throw;
  }
  delete[] filebuf;
}

void
AtlasDefinition::DebugPrint()
{
  itk::NumberToString<double> doubleConvert;
  std::cout << "<Atlas>" << std::endl;

  for (auto & elem : m_TemplateVolumes)
  {
    std::cout << "  <AtlasImage>" << std::endl
              << "    <type>" << elem.first << "</type>" << std::endl
              << "    <filename>" << elem.second << "</filename>" << std::endl
              << "  </AtlasImage>" << std::endl;
  }

  std::cout << "  <BrainMask>" << std::endl
            << "    <filename>" << m_TemplateBrainMask << "</filename>" << std::endl
            << "  </BrainMask>" << std::endl;

  std::cout << "  <HeadRegion>" << std::endl
            << "    <filename>" << m_TemplateHeadRegion << "</filename>" << std::endl
            << "  </HeadRegion>" << std::endl;
  for (auto & elem : m_PriorMap)
  {
    std::cout << "  <Prior>" << std::endl
              << "    <type>" << elem.first << "</type>" << std::endl
              << "    <filename>" << elem.second.GetFilename() << "</filename>" << std::endl
              << "    <Weight>" << doubleConvert(elem.second.GetWeight()) << "</Weight>" << std::endl
              << "    <GaussianClusterCount>" << elem.second.GetGaussianClusterCount() << "</GaussianClusterCount>"
              << std::endl
              << "    <LabelCode>" << elem.second.GetLabelCode() << "</LabelCode>" << std::endl
              << "    <UseForBias>" << elem.second.GetUseForBias() << "</UseForBias>" << std::endl
              << "    <IsForegroundPrior>" << elem.second.GetIsForegroundPrior() << "</IsForegroundPrior>" << std::endl;
    const BoundsMapType & curBounds = elem.second.GetBoundsList();
    for (const auto & curBound : curBounds)
    {
      std::cout << "    <bounds>" << std::endl
                << "      <type>" << curBound.first << "</type>" << std::endl
                << "      <lower>" << doubleConvert(curBound.second.GetLower()) << "</lower>" << std::endl
                << "      <upper>" << doubleConvert(curBound.second.GetUpper()) << "<upper>" << std::endl
                << "    </bounds>" << std::endl;
    }
    std::cout << "  </Prior>" << std::endl;
  }
  std::cout << "</Atlas>" << std::endl;
}
