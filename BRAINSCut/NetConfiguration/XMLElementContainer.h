#ifndef XMLElementContainer_H
#define XMLElementContainer_H
#include <string>
#include <iostream>

class XMLElementContainer
{
public:
  XMLElementContainer(const std::string & name)
  {
    this->SetName(name);
  }

  XMLElementContainer()
  {
  }

  virtual ~XMLElementContainer()
  {
  }

  virtual bool Verify() const = 0;

  std::string PrintSpaces(const int howmany) const
  {
    std::string spaces("");

    for( int i = 0; i < howmany; i++ )
      {
      spaces = spaces + " ";
      }
    return spaces;
  }

  virtual int PrintSelf(std::ostream & os, int indent) const = 0;

  const std::string & GetName() const
  {
    return m_Name;
  }

  void SetName(const std::string & s)
  {
    m_Name = s;
  }

private:
  std::string m_Name;
};

template <typename TOutputType>
class XMLContents :
  public XMLElementContainer
{
public:
  typedef XMLElementContainer SuperClass;
  typedef TOutputType         OutputType;
  virtual int PrintSelf(std::ostream &, int indent) const
  {
    // SuperClass::PrintSelf(os);
    // os << this->PrintSpaces(indent) << "=== XMLContents ===" <<
    // std::endl;
    // return indent+2;
    return indent;
  }

  XMLContents(const std::string & s) :
    XMLElementContainer(s)
  {
  }

  XMLContents()
  {
  }

  virtual ~XMLContents()
  {
  }

  virtual OutputType GetValue(void) const = 0;
};

#endif // XMLElementContainer_H
