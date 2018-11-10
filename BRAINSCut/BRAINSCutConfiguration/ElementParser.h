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
#ifndef ElementParser_h
#define ElementParser_h
#include "ElementContainer.h"
#include "BRAINSCutExceptionStringHandler.h"

#include "StringValue.h"
#include "BooleanValue.h"
#include "IntValue.h"
#include "FloatValue.h"

#include <map>
#include <vector>
#include <cctype>

#include "itkMacro.h" //Needed for nullptr

class ElementParser :
  public ElementContainer
{
private:
  using SuperClass = ElementContainer;
  void toLower(std::string & s) const
  {
    for( unsigned i = 0; i < s.size(); i++ )
      {
      s[i] = ::tolower(s[i]);
      }
  }

public:
  using MapType = std::map<std::string, ElementContainer *>;
  using iterator = MapType::iterator;
  typedef MapType::const_iterator                   const_iterator;
  using StringVectorType = std::vector<std::string>;

  ElementParser(const char *name) :
    ElementContainer(name)
  {
  }

  ElementParser()
  {
  }

  ~ElementParser() override
  {
    for( iterator it = m_Map.begin();
         it != m_Map.end(); ++it )
      {
      delete it->second;
      }
  }

  bool Verify() const override
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

  int PrintSelf(std::ostream & os, int indent) const override
  {
    // SuperClass::PrintSelf(os);
    // os << this->PrintSpaces(indent) << "=== ElementParser ===" <<
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

// AddElementByName
  void Add(ElementContainer *toAdd, const char *name)
  {
    std::string Name(name);

    this->toLower(Name);
    if( m_Map.find(Name) != m_Map.end() )
      {
      std::string msg("Duplicate object ");
      msg += name;
      toAdd->PrintSelf(std::cout, 2);
      throw BRAINSCutExceptionStringHandler(msg);
      }
    m_Map[Name] = toAdd;
  }

  void Add(ElementContainer *toAdd,
           const std::string & name)
  {
    this->Add( toAdd, name.c_str() );
  }

  template <typename T>
  const T * Get(const char *name) const
  {
    std::string lowercaseName(name);

    this->toLower(lowercaseName);
    MapType::const_iterator it = m_Map.find(lowercaseName);
    if( it == m_Map.end() )
      {
      std::cout << "ERROR:  map does not contain element for "
                << lowercaseName << std::endl;
      return nullptr;
      }

    return dynamic_cast<const T *>( it->second );     // returns zero if it can't
    // cast to the required type.
  }

  template <typename T>
  T * Get(const char *name)
  {
    std::string lowercaseName(name);

    this->toLower(lowercaseName);
    MapType::iterator it = m_Map.find(lowercaseName);
    if( it == m_Map.end() )
      {
      std::cout << "ERROR:  map does not contain element for "
                << lowercaseName << std::endl;
      return nullptr;
      }

    return dynamic_cast<T *>( it->second );     // returns zero if it can't
                                                // cast to the required type.
  }

  template <typename T, typename ValType>
  void SetAttribute(const char *name, const ValType & val)
  {
    T *t = this->Get<T>(name);

    if( t == nullptr )
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

    if( t == nullptr )
      {
      //      return std::string();
      std::string err("Can't find ");
      err += name;
      throw BRAINSCutExceptionStringHandler(err);
      }
    return t->GetValue();
  }

  template <typename T>
  typename T::ReturnType GetAttribute(const char *name)
  {
    T *t = this->Get<T>(name);

    if( t == nullptr )
      {
      //      return std::string();
      std::string err("Can't find ");
      err += name;
      throw BRAINSCutExceptionStringHandler(err);
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
  template <typename T>
  // type
  StringVectorType CollectAttValues(const char *attributeName) const
  {
    StringVectorType rval( this->size() );
    const_iterator   it;
    unsigned         i = 0;

    for( it = this->begin(); it != this->end(); ++it )
      {
      T *current = dynamic_cast<T *>( it->second );
      if( current == nullptr )           // skip everything not of the
        {                          // requested type
        continue;
        }
      const std::string val
        ( current->template GetAttribute<StringValue>(attributeName) );
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
  template <typename T>
  // list element type
  const T * GetMatching(const char *attName, const char *attValue) const
  {
    const_iterator it;

    for( it = this->begin(); it != this->end(); ++it )
      {
      const T *current = dynamic_cast<const T *>( it->second );
      if( current == nullptr )
        {
        continue;
        }
      if( current->template GetAttribute<StringValue>(attName) ==
          attValue )
        {
        return current;
        }
      }
    return nullptr;
  }

  //
  // given an attribute with a particular value, return the
  // list element matching that value
  template <typename T>
  // list element type
  T * GetMatching(const char *attName, const char *attValue)
  {
    iterator it;

    for( it = this->begin(); it != this->end(); ++it )
      {
      T *current = dynamic_cast<T *>( it->second );
      if( current == nullptr )
        {
        continue;
        }
      if( current->template GetAttribute<StringValue>(attName) ==
          attValue )
        {
        return current;
        }
      }
    return nullptr;
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

#endif // ElementParser_h
