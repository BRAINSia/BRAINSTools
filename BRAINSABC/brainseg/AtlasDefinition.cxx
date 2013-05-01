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
  if( this->m_TissueTypes.empty() )
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

double
AtlasDefinition
::StrToD(const char *str, const char *message) const
{
  char * last;
  double rval = strtod(str, &last);

  if( str == static_cast<const char *>(last) )
    {
    std::cerr << message << ' ' << this->m_LastXMLString << std::endl;
    throw;
    }
  return rval;
}

long
AtlasDefinition
::StrToL(const char *str, const char *message) const
{
  char *last;
  long  rval = strtol(str, &last, 10);

  if( str == static_cast<const char *>(last) )
    {
    std::cerr << message << ' ' << this->m_LastXMLString << std::endl;
    throw;
    }
  return rval;
}

void
AtlasDefinition::XMLEnd(const char *el)
{
  std::string El(el);
  const char *start = this->m_LastXMLString.c_str();

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
    this->m_LastWeight = this->StrToD(start, "Bad Weight given");
    }
  else if( El == "lower" )
    {
    this->m_LastLower = this->StrToD(start, "Bad lower bound");
    }
  else if( El == "upper" )
    {
    this->m_LastUpper = this->StrToD(start, "Bad upper bound");
    }
  else if( El == "GaussianClusterCount" )
    {
    this->m_LastGaussianClusterCount = this->StrToL(start, "Bad GaussianClusterCount");
    }
  else if( El == "LabelCode" )
    {
    this->m_LastLabelCode = this->StrToL(start, "Bad LabelCode");
    }
  else if( El == "UseForBias" )
    {
    this->m_LastUseForBias = this->StrToL(start, "Bad UseForBias");
    }
  else if( El == "IsForegroundPrior" )
    {
    this->m_LastIsForegroundPrior = this->StrToL(start, "Bad IsForegroundPrior");
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

  int parserReturn(1);
  try
    {
    parserReturn = XML_Parse(parser, filebuf, fSize, 1);
    }
  catch( ... )
    {
    parserReturn = 0;
    }
  if( !parserReturn )
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
              << doubleConvert(it->second.GetWeight() )
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
                << doubleConvert(it2->second.GetLower() )
                << "</lower>" << std::endl
                << "      <upper>"
                << doubleConvert(it2->second.GetUpper() )
                << "<upper>" << std::endl
                << "    </bounds>" << std::endl;
      }
    std::cout << "  </Prior>" << std::endl;
    }
  std::cout << "</Atlas>" << std::endl;
}
