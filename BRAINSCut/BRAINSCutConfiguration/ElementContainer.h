#ifndef ElementContainer_H
#define ElementContainer_H
#include <string>
#include <iostream>

class ElementContainer
{
public:
  ElementContainer(const std::string & name)
  {
    this->SetName(name);
  }

  ElementContainer()
  {
  }

  virtual ~ElementContainer()
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
  public ElementContainer
{
public:
  typedef ElementContainer SuperClass;
  typedef TOutputType      OutputType;
  virtual int PrintSelf(std::ostream &, int indent) const
  {
    // SuperClass::PrintSelf(os);
    // os << this->PrintSpaces(indent) << "=== XMLContents ===" <<
    // std::endl;
    // return indent+2;
    return indent;
  }

  XMLContents(const std::string & s) :
    ElementContainer(s)
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

#endif // ElementContainer_H
