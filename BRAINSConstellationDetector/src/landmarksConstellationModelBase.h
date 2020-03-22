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
#ifndef _landmarksConstellationModelBase_h
#define _landmarksConstellationModelBase_h

class landmarksConstellationModelBase
{
protected:
  using ValMapType = std::map<std::string, float>;
  using ValMapConstIterator = ValMapType::const_iterator;
  using ValMapIterator = ValMapType::iterator;

public:
  landmarksConstellationModelBase() = default;

  virtual ~landmarksConstellationModelBase() = default;

  virtual unsigned int
  GetNumDataSets() const
  {
    return m_NumDataSets;
  }

  virtual unsigned int
  GetSearchboxDims() const
  {
    return m_SearchboxDims;
  }

  virtual float
  GetResolutionUnits() const
  {
    return m_ResolutionUnits;
  }

  virtual const ValMapType &
  GetRadii() const
  {
    return this->m_Radius;
  }

  virtual float
  GetRadius(const std::string & PointName) const
  {
    ValMapConstIterator it(m_Radius.find(PointName));

    if (it == m_Radius.end())
    {
      throw;
    }
    return it->second;
  }

  virtual float
  GetHeight(const std::string & PointName) const
  {
    ValMapConstIterator it(m_Height.find(PointName));

    if (it == m_Height.end())
    {
      throw;
    }
    return it->second;
  }

  virtual void
  SetRadius(const std::string & PointName, float x)
  {
    m_Radius[PointName] = x;
  }

  virtual void
  SetHeight(const std::string & PointName, float x)
  {
    m_Height[PointName] = x;
  }

  virtual float
  GetInitialRotationAngle() const
  {
    return m_InitialRotationAngle;
  }

  virtual float
  GetInitialRotationStep() const
  {
    return m_InitialRotationStep;
  }

  virtual unsigned int
  GetNumRotationSteps() const
  {
    return m_NumRotationSteps;
  }

  virtual void
  SetNumDataSets(unsigned int x)
  {
    m_NumDataSets = x;
  }

  virtual void
  SetSearchboxDims(unsigned int x)
  {
    m_SearchboxDims = x;
  }

  virtual void
  SetResolutionUnits(float x)
  {
    m_ResolutionUnits = x;
  }

  virtual void
  SetInitialRotationAngle(float x)
  {
    m_InitialRotationAngle = x;
  }

  virtual void
  SetInitialRotationStep(float x)
  {
    m_InitialRotationStep = x;
  }

  virtual void
  SetNumRotationSteps(unsigned int x)
  {
    m_NumRotationSteps = x;
  }

protected:
  unsigned int m_NumDataSets{ 0 };
  unsigned int m_SearchboxDims{ 0 };
  float        m_ResolutionUnits{ 1.0 };
  ValMapType   m_Height;
  ValMapType   m_Radius;
  float        m_InitialRotationAngle{ 0.0 };
  float        m_InitialRotationStep{ 0.0 };
  unsigned int m_NumRotationSteps{ 0 };
};

#endif // _landmarksConstellationModelBase_h
