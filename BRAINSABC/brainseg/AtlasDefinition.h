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
#ifndef __AtlasDefinition_h
#define __AtlasDefinition_h

#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include "itkNumberToString.h"
/** \class AtlasDefiniton */
class AtlasDefinition
{
public:
  //  static const char * tissueTypes[];
  using TissueTypeVector = std::vector<std::string>;
  using TemplateMap = std::map<std::string, std::string>;
  using FloatVector = std::vector<float>;
  using IntVector = std::vector<int>;
  /** \class BoundsType */
  class BoundsType
  {
  public:
    BoundsType() = default;

    void
    SetLower(const double v)
    {
      m_Low = v;
    }

    double
    GetLower() const
    {
      return m_Low;
    }

    void
    SetUpper(const double v)
    {
      m_High = v;
    }

    double
    GetUpper() const
    {
      return m_High;
    }

    void
    Print() const
    {
      itk::NumberToString<double> ds;
      std::cout << "RANGE:  [" << ds(m_Low) << "," << ds(m_High) << "]" << std::endl;
    }

  private:
    double m_Low{ 0 };
    double m_High{ 0 };
  };
  using BoundsMapType = std::map<std::string, BoundsType>;

  AtlasDefinition();
  void
  InitFromXML(const std::string & XMLFilename);

  void
  DebugPrint();

  const TemplateMap &
  GetTemplateVolumes() const
  {
    return m_TemplateVolumes;
  }

  const std::string &
  GetTemplateBrainMask() const
  {
    return m_TemplateBrainMask;
  }

  const std::string &
  GetTemplateHeadRegion() const
  {
    return m_TemplateHeadRegion;
  }

  std::string
  GetPriorFilename(const std::string & tissueType) const
  {
    if (m_PriorMap.find(tissueType) == m_PriorMap.end())
    {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
    }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if (mit == m_PriorMap.end())
    {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
    }
    return mit->second.GetFilename();
  }

  double
  GetWeight(const std::string & tissueType) const
  {
    if (m_PriorMap.find(tissueType) == m_PriorMap.end())
    {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
    }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if (mit == m_PriorMap.end())
    {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
    }
    return mit->second.GetWeight();
  }

  int
  GetGaussianClusterCount(const std::string & tissueType) const
  {
    if (m_PriorMap.find(tissueType) == m_PriorMap.end())
    {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
    }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if (mit == m_PriorMap.end())
    {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
    }
    return mit->second.GetGaussianClusterCount();
  }

  int
  GetLabelCode(const std::string & tissueType) const
  {
    // HACK:  All the get functions need review to remove duplicate code.
    if (m_PriorMap.find(tissueType) == m_PriorMap.end())
    {
      std::cout << "MISSING LABEL CODE IN ATLAS:  " << tissueType << std::endl;
      throw;
    }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if (mit == m_PriorMap.end())
    {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
    }
    return mit->second.GetLabelCode();
  }

  int
  GetUseForBias(const std::string & tissueType) const
  {
    if (m_PriorMap.find(tissueType) == m_PriorMap.end())
    {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
    }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if (mit == m_PriorMap.end())
    {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
    }
    return mit->second.GetUseForBias();
  }

  bool
  GetIsForegroundPrior(const std::string & tissueType) const
  {
    if (m_PriorMap.find(tissueType) == m_PriorMap.end())
    {
      std::cout << "MISSING IsForegroungPrior IN ATLAS:  " << tissueType << std::endl;
      throw;
    }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if (mit == m_PriorMap.end())
    {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
    }
    return mit->second.GetIsForegroundPrior();
  }

  const BoundsType &
  GetBounds(const std::string & tissueType, const std::string & Modality) const
  {
    if (m_PriorMap.find(tissueType) == m_PriorMap.end())
    {
      std::cout << "MISSING TISSUE TYPE IN ATLAS:  " << tissueType << std::endl;
      throw;
    }
    PriorMapType::const_iterator mit = m_PriorMap.find(tissueType);
    if (mit == m_PriorMap.end())
    {
      // HACK:  Should throw and exception here with line number and file
      std::cout << "ERROR:  Invalid tissueType requested in GetPriorFilename" << std::endl;
      throw;
    }
    return mit->second.GetBounds(Modality);
  }

  double
  GetLow(const std::string & tissueType, const std::string & Modality) const
  {
    return this->GetBounds(tissueType, Modality).GetLower();
  }

  double
  GetHigh(const std::string & tissueType, const std::string & Modality) const
  {
    return this->GetBounds(tissueType, Modality).GetUpper();
  }

  const TissueTypeVector &
  TissueTypes();

  void
  XMLStart(const char * el);

  void
  XMLEnd(const char * el);

  void
  XMLChar(const char * buf);

private:
  /** convert to float or throw exception on failure */
  double
  StrToD(const char * str, const char * message) const;

  /** convert to long or throw exception on failure */
  long
  StrToL(const char * str, const char * message) const;

  TemplateMap      m_TemplateVolumes;
  std::string      m_TemplateBrainMask;
  std::string      m_TemplateHeadRegion;
  TissueTypeVector m_TissueTypes;
  class Prior
  {
  public:
    Prior() = default;

    const std::string &
    GetFilename() const
    {
      return m_Filename;
    }

    void
    SetFilename(const std::string & s)
    {
      m_Filename = s;
    }

    double
    GetWeight() const
    {
      return m_Weight;
    }

    void
    SetWeight(double x)
    {
      m_Weight = x;
    }

    int
    GetGaussianClusterCount() const
    {
      return m_GaussianClusterCount;
    }

    void
    SetGaussianClusterCount(int i)
    {
      m_GaussianClusterCount = i;
    }

    int
    GetLabelCode() const
    {
      return m_LabelCode;
    }

    void
    SetLabelCode(const int i)
    {
      m_LabelCode = i;
    }

    int
    GetUseForBias() const
    {
      return m_UseForBias;
    }

    void
    SetUseForBias(const bool i)
    {
      m_UseForBias = i;
    }

    int
    GetIsForegroundPrior() const
    {
      return m_IsForegroundPrior;
    }

    void
    SetIsForegroundPrior(const bool i)
    {
      m_IsForegroundPrior = i;
    }

    const BoundsType &
    GetBounds(const std::string & Modality) const
    {
      BoundsMapType::const_iterator bit = m_BoundsMap.find(Modality);

      if (bit == m_BoundsMap.end())
      {
        std::cout << "MISSING MODALIITY TYPE IN ATLAS:  " << Modality << std::endl;
        throw;
      }

      return bit->second;
    }

    void
    SetBounds(const std::string & Modality, const BoundsType & b)
    {
      this->m_BoundsMap[Modality] = b;
    }

    void
    SetBounds(const std::string & Modality, double lower, double upper)
    {
      BoundsType b;

      b.SetLower(lower);
      b.SetUpper(upper);
      this->SetBounds(Modality, b);
    }

    void
    SetBoundsList(BoundsMapType & boundsMap)
    {
      this->m_BoundsMap = boundsMap;
    }

    const BoundsMapType &
    GetBoundsList() const
    {
      return this->m_BoundsMap;
    }

  private:
    std::string   m_Filename;
    double        m_Weight{ 0.0 };
    int           m_GaussianClusterCount{ 0 };
    int           m_LabelCode{ 0 };
    bool          m_UseForBias{ false };
    bool          m_IsForegroundPrior{ false };
    BoundsMapType m_BoundsMap;
  };
  using PriorMapType = std::map<std::string, Prior>;
  PriorMapType m_PriorMap;
  //
  // XML Parsing variables
  std::string            m_LastXMLString;
  std::string            m_LastPriorType;
  std::string            m_LastFilename;
  std::string            m_LastType;
  std::list<std::string> m_XMLElementStack;

  double m_LastWeight{ 0.0 };
  double m_LastLower{ 0.0 };
  double m_LastUpper{ 0.0 };

  int           m_LastGaussianClusterCount{ 0 };
  int           m_LastLabelCode{ 0 };
  bool          m_LastUseForBias{ false };
  bool          m_LastIsForegroundPrior{ false };
  BoundsMapType m_LastPriorBounds;
};

#endif // AtlasDefinition_h
