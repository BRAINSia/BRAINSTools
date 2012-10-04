#include "AtlasDefinition.h"
#include <expat.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <itksys/SystemTools.hxx>

namespace // anon
{
extern "C"
{
namespace AtlasXMLParser
{
void
XMLstart(void *data, const char *el, const char * *)
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
  AtlasDefinition *_this = reinterpret_cast<AtlasDefinition *>( data );

  _this->XMLStart(el);
}

void
XMLend(void *data, const char *el)
{
  // std::cerr << "End, El = ("
  //         << el
  //         << ")" << std::endl;
  AtlasDefinition *_this = reinterpret_cast<AtlasDefinition *>( data );

  _this->XMLEnd(el);
}

void
XMLcharhandler(void *data, const char *txt, int txtlen)
{
  // std::cerr << "Char data = (";
  // for(unsigned i = 0; i < txtlen; i++)
  //   {
  //   std::cerr << txt[i];
  //   }
  // std::cerr << ")" << std::endl;
  char *buf = new char[txtlen + 1];
  int   i;

  for( i = 0; i < txtlen; i++ )
    {
    buf[i] = txt[i];
    }
  buf[i] = '\0';
  AtlasDefinition *_this = reinterpret_cast<AtlasDefinition *>( data );
  _this->XMLChar(buf);
  delete[] buf;
}
}
}
}
#if 0
const char *AtlasDefinition::tissueTypes[] =
  {
  "WM",
  "GM",
  "BGM",
  "CSF",
  "VB",
  "NOTCSF",
  "NOTGM",
  "NOTVB",
  "NOTWM",
  "AIR",
  0
  };
#endif

AtlasDefinition::AtlasDefinition() :
  m_LastWeight(0.0),
  m_LastLower(0.0),
  m_LastUpper(0.0),
  m_LastGaussianClusterCount(0),
  m_LastLabelCode(0),
  m_LastUseForBias(0),
  m_LastIsForegroundPrior(0)
{
  this->m_TissueTypes.resize(0);
}

const AtlasDefinition::TissueTypeVector &
AtlasDefinition::TissueTypes()
{
  if( this->m_TissueTypes.size() == 0 )
    {
    for( PriorMapType::const_iterator it = this->m_PriorMap.begin();
         it != this->m_PriorMap.end();
         ++it )
      {
      this->m_TissueTypes.push_back(it->first);
      }
    }
  return this->m_TissueTypes;
}

void
AtlasDefinition::XMLStart(const char *el)
{
  std::string El(el);

  this->m_XMLElementStack.push_back(El);
  if( El == "Prior" )
    {
    this->m_LastPriorBounds.clear();
    }
}

void
AtlasDefinition::XMLEnd(const char *el)
{
  std::string El(el);
  const char *start = this->m_LastXMLString.c_str();
  char *      last;

  // pop the current element name off the stack.
  this->m_XMLElementStack.pop_back();

  if( El == "Atlas" )
    {
    return; // final element
    }
  const std::string ContainingElement = this->m_XMLElementStack.back();
  if( El == "filename" )
    {
    this->m_LastFilename = this->m_LastXMLString;
    }
  else if( El == "type" )
    {
    if( ContainingElement == "Prior" )
      {
      this->m_LastPriorType = this->m_LastXMLString;
      }
    else
      {
      this->m_LastType = this->m_LastXMLString;
      }
    }
  else if( El == "AtlasImage" )
    {
    this->m_TemplateVolumes[this->m_LastType] = this->m_LastFilename;
    }
  else if( El == "Weight" )
    {
    double w = strtod(start, &last);
    if( start == (const char *)last )
      {
      std::cerr << "Bad Weight given " << this->m_LastXMLString
                << std::endl;
      throw;
      }
    this->m_LastWeight = w;
    }
  else if( El == "lower" )
    {
    double l = strtod(start, &last);
    if( start == (const char *)last )
      {
      std::cerr << "Bad lower bound " << this->m_LastXMLString
                << std::endl;
      throw;
      }
    this->m_LastLower = l;
    }
  else if( El == "upper" )
    {
    double u = strtod(start, &last);
    if( start == (const char *)last )
      {
      std::cerr << "Bad upper bound " << this->m_LastXMLString
                << std::endl;
      throw;
      }
    this->m_LastUpper = u;
    }
  else if( El == "GaussianClusterCount" )
    {
    int GCC = strtol(start, &last, 10);
    if( start == (const char *)last )
      {
      std::cerr << "Bad GaussianClusterCount given " << this->m_LastXMLString
                << std::endl;
      throw;
      }
    this->m_LastGaussianClusterCount = GCC;
    }
  else if( El == "LabelCode" )
    {
    int GCC = strtol(start, &last, 10);
    if( start == (const char *)last )
      {
      std::cerr << "Bad LabelCode given " << this->m_LastXMLString
                << std::endl;
      throw;
      }
    this->m_LastLabelCode = GCC;
    }
  else if( El == "UseForBias" )
    {
    int UFB = strtol(start, &last, 10);
    if( start == (const char *)last )
      {
      std::cerr << "Bad UseForBias given " << this->m_LastXMLString
                << std::endl;
      throw;
      }
    this->m_LastUseForBias = UFB;
    }
  else if( El == "IsForegroundPrior" )
    {
    int UFB = strtol(start, &last, 10);
    if( start == (const char *)last )
      {
      std::cerr << "Bad IsForegroundPrior given " << this->m_LastXMLString
                << std::endl;
      throw;
      }
    this->m_LastIsForegroundPrior = UFB;
    }
  else if( El == "BrainMask" )
    {
    this->m_TemplateBrainMask = this->m_LastFilename;
    }
  else if( El == "HeadRegion" )
    {
    this->m_TemplateHeadRegion = this->m_LastFilename;
    }
  else if( El == "Prior" )
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
  else if( El == "bounds" )
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
AtlasDefinition::XMLChar(const char *buf)
{
  this->m_LastXMLString = buf;
}

void
AtlasDefinition::InitFromXML(const std::string & XMLFilename)
{
  std::ifstream xmlFile(XMLFilename.c_str(),
                        std::ifstream::in);

  if( !xmlFile.good() )
    {
    throw;
    }
  std::streamsize fSize =
    itksys::SystemTools::FileLength(XMLFilename.c_str() );

  XML_Parser parser = XML_ParserCreate(0);
  XML_SetUserData( parser, static_cast<void *>( this ) );
  XML_SetElementHandler(parser, AtlasXMLParser::XMLstart, AtlasXMLParser::XMLend);
  XML_SetCharacterDataHandler(parser, AtlasXMLParser::XMLcharhandler);

  char *filebuf = new char[fSize];
  if( filebuf == NULL )
    {
    throw;
    }

  xmlFile.read(filebuf, fSize);
  if( static_cast<std::streamsize>(xmlFile.gcount() ) != fSize )
    {
    delete [] filebuf;
    throw;
    }
  xmlFile.close();

  if( XML_Parse(parser, filebuf, fSize, 1) == 0 )
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
  std::cout << "<Atlas>" << std::endl;

  for( TemplateMap::iterator it
         = m_TemplateVolumes.begin();
       it != m_TemplateVolumes.end();
       ++it )
    {
    std::cout << "  <AtlasImage>" << std::endl
              << "    <type>" << it->first << "</type>" << std::endl
              << "    <filename>" << it->second << "</filename>" << std::endl
              << "  </AtlasImage>" << std::endl;
    }

  std::cout << "  <BrainMask>" << std::endl
            << "    <filename>"
            << m_TemplateBrainMask
            << "</filename>" << std::endl
            << "  </BrainMask>" << std::endl;

  std::cout << "  <HeadRegion>" << std::endl
            << "    <filename>"
            << m_TemplateHeadRegion
            << "</filename>" << std::endl
            << "  </HeadRegion>" << std::endl;
  for( PriorMapType::iterator it
         = m_PriorMap.begin();
       it != m_PriorMap.end();
       ++it )
    {
    std::cout << "  <Prior>" << std::endl
              << "    <type>"
              << it->first
              << "</type>" << std::endl
              << "    <filename>"
              << it->second.GetFilename()
              << "</filename>" << std::endl
              << "    <Weight>"
              << it->second.GetWeight()
              << "</Weight>" << std::endl
              << "    <GaussianClusterCount>"
              << it->second.GetGaussianClusterCount()
              << "</GaussianClusterCount>" << std::endl
              << "    <LabelCode>"
              << it->second.GetLabelCode()
              << "</LabelCode>" << std::endl
              << "    <UseForBias>"
              << it->second.GetUseForBias()
              << "</UseForBias>" << std::endl
              << "    <IsForegroundPrior>"
              << it->second.GetIsForegroundPrior()
              << "</IsForegroundPrior>" << std::endl;
    const BoundsMapType & curBounds = it->second.GetBoundsList();
    for( BoundsMapType::const_iterator it2
           = curBounds.begin();
         it2 != curBounds.end();
         ++it2 )
      {
      std::cout << "    <bounds>" << std::endl
                << "      <type>"
                << it2->first
                << "</type>" << std::endl
                << "      <lower>"
                << it2->second.GetLower()
                << "</lower>" << std::endl
                << "      <upper>"
                << it2->second.GetUpper()
                << "<upper>" << std::endl
                << "    </bounds>" << std::endl;
      }
    std::cout << "  </Prior>" << std::endl;
    }
  std::cout << "</Atlas>" << std::endl;
}
