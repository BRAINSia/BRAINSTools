#ifndef __AtlasDefinition_h
#define __AtlasDefinition_h

#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include "DoubleToString.h"
/** \class AtlasDefiniton */
class AtlasDefinition
{
public:
  //  static const char * tissueTypes[];
  typedef std::vector<std::string>           TissueTypeVector;
  typedef std::map<std::string, std::string> TemplateMap;
  typedef std::vector<float>                 FloatVector;
  typedef std::vector<int>                   IntVector;
  /** \class BoundsType */
  class BoundsType
  {
public:
    BoundsType() : m_Low(0), m_High(0)
    {
    }

    void SetLower(const double v)
    {
      m_Low = v;
    }

    double GetLower() const
    {
      return m_Low;
    }

    void SetUpper(const double v)
    {
      m_High = v;
    }

    double GetUpper() const
    {
      return m_High;
    }

    void Print(void) const
    {
      DoubleToString ds;
      std::cout << "RANGE:  [" << ds(m_Low) << "," << ds(m_High) << "]" << std::endl;
    }

private:
    double m_Low;
    double m_High;
  };
  typedef std::map<std::string, BoundsType> BoundsMapType;

  AtlasDefinition();
  void InitFromXML(const std::string & XMLFilename);

  void DebugPrint();

  const TemplateMap & GetTemplateVolumes() const
  {
    return m_TemplateVolumes;
  }

  const std::string & GetTemplateBrainMask() const
  {
    return m_TemplateBrainMask;
  }

  const std::string & GetTemplateHeadRegion() const
  {
    return m_TemplateHeadRegion;
  }

  std::string GetPriorFilename(const std::string & tissueType) const
  {
    if( m_PriorMap.find(tissueType) == m_PriorMap.end() )
      {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
      }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if( mit == m_PriorMap.end() )
      {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
      }
    return mit->second.GetFilename();
  }

  double GetWeight(const std::string & tissueType) const
  {
    if( m_PriorMap.find(tissueType) == m_PriorMap.end() )
      {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
      }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if( mit == m_PriorMap.end() )
      {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
      }
    return mit->second.GetWeight();
  }

  int GetGaussianClusterCount(const std::string & tissueType) const
  {
    if( m_PriorMap.find(tissueType) == m_PriorMap.end() )
      {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
      }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if( mit == m_PriorMap.end() )
      {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
      }
    return mit->second.GetGaussianClusterCount();
  }

  int GetLabelCode(const std::string & tissueType) const
  {
    // HACK:  All the get functions need review to remove duplicate code.
    if( m_PriorMap.find(tissueType) == m_PriorMap.end() )
      {
      std::cout << "MISSING LABEL CODE IN ATLAS:  " << tissueType << std::endl;
      throw;
      }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if( mit == m_PriorMap.end() )
      {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
      }
    return mit->second.GetLabelCode();
  }

  int GetUseForBias(const std::string & tissueType) const
  {
    if( m_PriorMap.find(tissueType) == m_PriorMap.end() )
      {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
      }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if( mit == m_PriorMap.end() )
      {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
      }
    return mit->second.GetUseForBias();
  }

  bool GetIsForegroundPrior(const std::string & tissueType) const
  {
    if( m_PriorMap.find(tissueType) == m_PriorMap.end() )
      {
      std::cout << "MISSING IsForegroungPrior IN ATLAS:  " << tissueType << std::endl;
      throw;
      }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if( mit == m_PriorMap.end() )
      {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
      }
    return mit->second.GetIsForegroundPrior();
  }

  const BoundsType & GetBounds(const std::string & tissueType,
                               const std::string & Modality) const
  {
    if( m_PriorMap.find(tissueType) == m_PriorMap.end() )
      {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
      }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if( mit == m_PriorMap.end() )
      {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
      }
    return mit->second.GetBounds(Modality);
  }

  double GetLow(const std::string & tissueType,
                const std::string & Modality) const
  {
    return this->GetBounds(tissueType, Modality).GetLower();
  }

  double GetHigh(const std::string & tissueType,
                 const std::string & Modality) const
  {
    return this->GetBounds(tissueType, Modality).GetUpper();
  }

  const TissueTypeVector & TissueTypes();

  void XMLStart(const char *el);

  void XMLEnd(const char *el);

  void XMLChar(const char *buf);

private:
  /** convert to float or throw exception on failure */
  double StrToD(const char *str, const char *message) const;

  /** convert to long or throw exception on failure */
  long   StrToL(const char *str, const char *message) const;

  TemplateMap      m_TemplateVolumes;
  std::string      m_TemplateBrainMask;
  std::string      m_TemplateHeadRegion;
  TissueTypeVector m_TissueTypes;
  class Prior
  {
public:
    Prior() :
      m_Weight(0.0),
      m_GaussianClusterCount(0),
      m_LabelCode(0),
      m_UseForBias(false),
      m_IsForegroundPrior(false)
    {
    }

    const std::string & GetFilename() const
    {
      return m_Filename;
    }

    void SetFilename(const std::string & s)
    {
      m_Filename = s;
    }

    double GetWeight() const
    {
      return m_Weight;
    }

    void SetWeight(double x)
    {
      m_Weight = x;
    }

    int GetGaussianClusterCount() const
    {
      return m_GaussianClusterCount;
    }

    void SetGaussianClusterCount(int i)
    {
      m_GaussianClusterCount = i;
    }

    int GetLabelCode() const
    {
      return m_LabelCode;
    }

    void SetLabelCode(const int i)
    {
      m_LabelCode = i;
    }

    int GetUseForBias() const
    {
      return m_UseForBias;
    }

    void SetUseForBias(const bool i)
    {
      m_UseForBias = i;
    }

    int GetIsForegroundPrior() const
    {
      return m_IsForegroundPrior;
    }

    void SetIsForegroundPrior(const bool i)
    {
      m_IsForegroundPrior = i;
    }

    const BoundsType & GetBounds(const std::string & Modality) const
    {
      BoundsMapType::const_iterator bit = m_BoundsMap.find(Modality);

      if( bit == m_BoundsMap.end() )
        {
        std::cout << "MISSING MODALIITY TYPE IN ATLAS:  " << Modality << std::endl;
        throw;
        }

      return bit->second;
    }

    void SetBounds(const std::string & Modality, const BoundsType & b)
    {
      this->m_BoundsMap[Modality] = b;
    }

    void SetBounds(const std::string & Modality, double lower, double upper)
    {
      BoundsType b;

      b.SetLower(lower);
      b.SetUpper(upper);
      this->SetBounds(Modality, b);
    }

    void SetBoundsList(BoundsMapType & boundsMap)
    {
      this->m_BoundsMap = boundsMap;
    }

    const BoundsMapType & GetBoundsList() const
    {
      return this->m_BoundsMap;
    }

private:
    std::string   m_Filename;
    double        m_Weight;
    int           m_GaussianClusterCount;
    int           m_LabelCode;
    bool          m_UseForBias;
    bool          m_IsForegroundPrior;
    BoundsMapType m_BoundsMap;
  };
  typedef std::map<std::string, Prior> PriorMapType;
  PriorMapType m_PriorMap;
  //
  // XML Parsing variables
  std::string            m_LastXMLString;
  std::string            m_LastPriorType;
  std::string            m_LastFilename;
  std::string            m_LastType;
  std::list<std::string> m_XMLElementStack;

  double m_LastWeight;
  double m_LastLower;
  double m_LastUpper;

  int           m_LastGaussianClusterCount;
  int           m_LastLabelCode;
  bool          m_LastUseForBias;
  bool          m_LastIsForegroundPrior;
  BoundsMapType m_LastPriorBounds;
};

#endif // AtlasDefinition_h
