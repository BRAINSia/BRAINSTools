#ifndef XMLParser_h
#define XMLParser_h
#include <string>
#include <expat.h>

class XMLParser
{
public:
  XMLParser(const std::string & filename) : m_Filename(filename),
    m_UserData(0),
    m_Buffer(0)
  {
  }

  void SetUserData(void *userData)
  {
    m_UserData = userData;
  }

  void * GetUserData()
  {
    return m_UserData;
  }

  bool Parse();

  virtual void StartElement(void *userData, const XML_Char *name, const XML_Char * *atts) = 0;

  virtual void EndElement(void *userData, const XML_Char *name) = 0;

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
