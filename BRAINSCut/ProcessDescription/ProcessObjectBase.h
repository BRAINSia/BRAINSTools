#ifndef ProcessObjectBase_H
#define ProcessObjectBase_H
#include <string>
#include <iostream>

class ProcessObjectException
{
public:
  ProcessObjectException(const char *errorString)
  {
    this->m_ErrorString = errorString;
  }

  ProcessObjectException(const std::string & errorString)
  {
    this->m_ErrorString = errorString;
  }

  const std::string & Error() const
  {
    return m_ErrorString;
  }

private:
  std::string m_ErrorString;
};

class ProcessObjectBase
{
public:
  ProcessObjectBase(const std::string & name)
  {
    this->SetName(name);
  }

  ProcessObjectBase()
  {
  }

  virtual ~ProcessObjectBase()
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

#if 0
    {
    os << "===ProcessObjectBase===!" << std::endl;
    }
#endif
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

template <typename OutputType>
class ValueObjectBase :
  public ProcessObjectBase
{
public:
  typedef ProcessObjectBase SuperClass;
  virtual int PrintSelf(std::ostream &, int indent) const
  {
    // SuperClass::PrintSelf(os);
    // os << this->PrintSpaces(indent) << "=== ValueObjectBase ===" <<
    // std::endl;
    // return indent+2;
    return indent;
  }

  ValueObjectBase(const std::string & s) :
    ProcessObjectBase(s)
  {
  }

  ValueObjectBase()
  {
  }

  virtual ~ValueObjectBase()
  {
  }

  virtual OutputType GetValue(void) const = 0;
};

#endif // ProcessObjectBase_H
