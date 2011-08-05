#ifndef CompoundObjectBase_h
#define CompoundObjectBase_h
#include "ProcessObjectBase.h"
#include "StringValue.h"
#include <map>
#include <vector>
#include <ctype.h>

class CompoundObjectBase :
  public ProcessObjectBase
{
private:
  typedef ProcessObjectBase SuperClass;
  void tolower(std::string & s) const
  {
    for( unsigned i = 0; i < s.size(); i++ )
      {
      s[i] = ::tolower(s[i]);
      }
  }

public:
  typedef std::map<std::string, ProcessObjectBase *> MapType;
  typedef MapType::iterator                          iterator;
  typedef MapType::const_iterator                    const_iterator;
  typedef std::vector<std::string>                   StringVectorType;
  CompoundObjectBase(const char *name) :
    ProcessObjectBase(name)
  {
  }

  CompoundObjectBase()
  {
  }

  virtual ~CompoundObjectBase()
  {
    for( iterator it = m_Map.begin();
         it != m_Map.end(); ++it )
      {
      delete it->second;
      }
  }

  virtual bool Verify() const
  {
    static bool             ErrorPrinted = false;
    MapType::const_iterator it;

    for( it = m_Map.begin(); it != m_Map.end(); ++it )
      {
      if( it->second->Verify() == false )
        {
        if( ErrorPrinted == false )
          {
          it->second->PrintSelf(std::cerr, 0);
          ErrorPrinted = true;
          }
        return false;
        }
      }
    return true;
  }

  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    // SuperClass::PrintSelf(os);
    // os << this->PrintSpaces(indent) << "=== CompoundObjectBase ===" <<
    // std::endl;
    MapType::const_iterator errdump;

    for( errdump = m_Map.begin(); errdump != m_Map.end(); ++errdump )
      {
      os << this->PrintSpaces(indent) << "--" << errdump->first << "--"
         << std::endl;
      errdump->second->PrintSelf(std::cerr, indent + 2);
      std::cerr << this->PrintSpaces(indent)
                << "--------------------------------" << std::endl;
      }
    return indent + 2;
  }

  void Add(ProcessObjectBase *toAdd, const char *name)
  {
    std::string Name(name);

    this->tolower(Name);
    if( m_Map.find(Name) != m_Map.end() )
      {
      std::string msg("Duplicate object ");
      msg += name;
      toAdd->PrintSelf(std::cout, 2);
      throw ProcessObjectException(msg);
      }
    m_Map[Name] = toAdd;
  }

  void Add(ProcessObjectBase *toAdd,
           const std::string & name)
  {
    this->Add( toAdd, name.c_str() );
  }

  template <typename T>
  const T * Get(const char *name) const
  {
    std::string lowercaseName(name);

    this->tolower(lowercaseName);
    MapType::const_iterator it = m_Map.find(lowercaseName);
    if( it == m_Map.end() )
      {
      //  std::cout << "ERROR:  map does not contain element for " <<
      // lowercaseName << std::endl;
      return 0;
      }

    return dynamic_cast<const T *>( it->second );     // returns zero if it can't
    // cast to the required type.
  }

  template <typename T>
  T * Get(const char *name)
  {
    std::string lowercaseName(name);

    this->tolower(lowercaseName);
    MapType::iterator it = m_Map.find(lowercaseName);
    if( it == m_Map.end() )
      {
      //  std::cout << "ERROR:  map does not contain element for " <<
      // lowercaseName << std::endl;
      return 0;
      }

    return dynamic_cast<T *>( it->second );     // returns zero if it can't
                                                // cast to the required type.
  }

  template <typename T, typename ValType>
  void SetAttribute(const char *name, const ValType & val)
  {
    T *t = this->Get<T>(name);

    if( t == 0 )
      {
      std::cout << "ERROR:  Can not set name " << name << " to value "
                << val << std::endl;
      }
    t->SetValue(val);
  }

  template <typename T>
  const typename T::ReturnType GetAttribute(const char *name) const
  {
    const T *t = this->Get<T>(name);

    if( t == 0 )
      {
      //      return std::string();
      std::string err("Can't find ");
      err += name;
      throw ProcessObjectException(err);
      }
    return t->GetValue();
  }

  template <typename T>
  typename T::ReturnType GetAttribute(const char *name)
  {
    T *t = this->Get<T>(name);

    if( t == 0 )
      {
      //      return std::string();
      std::string err("Can't find ");
      err += name;
      throw ProcessObjectException(err);
      }
    return t->GetValue();
  }

  template <typename T>
  typename T::ReturnType GetAttributeIfExist(const char *name)
  {
    T *t = this->Get<T>(name);

    if( t == 0 )
      {
      std::cout << " Can't find "
                << name
                << ", returning null."
                << std::endl;

      return 0;
      }
    return t->GetValue();
  }

  template <typename T>
  const typename T::ReturnType GetAttributeIfExist(const char *name) const
  {
    T *t = this->Get<T>(name);

    if( t == 0 )
      {
      std::cout << " Can't find "
                << name
                << ", returning null."
                << std::endl;

      return 0;
      }
    return t->GetValue();
  }

  //
  // given a particular attribute name, collect that
  // attributes. Only works on homogenous compound objects,
  // i.e. lists of a particular type
  template <class T>
  // type
  StringVectorType CollectAttValues(const char *attributeName) const
  {
    StringVectorType rval( this->size() );
    const_iterator   it;
    unsigned         i = 0;

    for( it = this->begin(); it != this->end(); ++it )
      {
      T *current = dynamic_cast<T *>( it->second );
      if( current == 0 )           // skip everything not of the
        {                          // requested type
        continue;
        }
      const std::string val
        ( current->GetAttribute<StringValue>(
          attributeName) );
      if( val != "" )
        {
        rval[i] = val;
        i++;
        }
      }
    return rval;
  }

  //
  // given an attribute with a particular value, return the
  // list element matching that value
  template <class T>
  // list element type
  const T * GetMatching(const char *attName, const char *attValue) const
  {
    const_iterator it;

    for( it = this->begin(); it != this->end(); ++it )
      {
      const T *current = dynamic_cast<const T *>( it->second );
      if( current == 0 )
        {
        continue;
        }
      if( current->GetAttribute<StringValue>(attName) ==
          attValue )
        {
        return current;
        }
      }
    return 0;
  }

  //
  // given an attribute with a particular value, return the
  // list element matching that value
  template <class T>
  // list element type
  T * GetMatching(const char *attName, const char *attValue)
  {
    iterator it;

    for( it = this->begin(); it != this->end(); ++it )
      {
      T *current = dynamic_cast<T *>( it->second );
      if( current == 0 )
        {
        continue;
        }
      if( current->GetAttribute<StringValue>(attName) ==
          attValue )
        {
        return current;
        }
      }
    return 0;
  }

  const_iterator begin() const
  {
    return m_Map.begin();
  }

  // cppcheck-suppress functionConst
  iterator begin()
  {
    return m_Map.begin();
  }

  const_iterator end() const
  {
    return m_Map.end();
  }

  // cppcheck-suppress functionConst
  iterator end()
  {
    return m_Map.end();
  }

  unsigned size() const
  {
    return m_Map.size();
  }

private:
  MapType m_Map;
};

#endif // CompoundObjectBase_h
