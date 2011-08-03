#ifndef _landmarksConstellationModelBase_h
#define  _landmarksConstellationModelBase_h

class landmarksConstellationModelBase
{
protected:
  typedef std::map<std::string, float> ValMapType;
  typedef ValMapType::const_iterator   ValMapConstIterator;
  typedef ValMapType::iterator         ValMapIterator;
public:
  landmarksConstellationModelBase() : m_NumDataSets(0),
    m_SearchboxDims(0),
    m_ResolutionUnits(1.0),
    m_InitialRotationAngle(0.0),
    m_InitialRotationStep(0.0),
    m_NumRotationSteps(0)
  {
  }

  virtual ~landmarksConstellationModelBase()
  {
  }

  virtual unsigned int GetNumDataSets() const
  {
    return m_NumDataSets;
  }

  virtual unsigned int GetSearchboxDims() const
  {
    return m_SearchboxDims;
  }

  virtual float GetResolutionUnits() const
  {
    return m_ResolutionUnits;
  }

  virtual const ValMapType & GetRadii() const
  {
    return this->m_Radius;
  }

  virtual float GetRadius(const std::string & PointName) const
  {
    ValMapConstIterator it( m_Radius.find(PointName) );

    if( it == m_Radius.end() )
      {
      throw;
      }
    return it->second;
  }

  virtual float GetHeight(const std::string & PointName) const
  {
    ValMapConstIterator it( m_Height.find(PointName) );

    if( it == m_Height.end() )
      {
      throw;
      }
    return it->second;
  }

  virtual void SetRadius(const std::string & PointName, float x)
  {
    m_Radius[PointName] = x;
  }

  virtual void SetHeight(const std::string & PointName, float x)
  {
    m_Height[PointName] = x;
  }

  virtual float GetInitialRotationAngle() const
  {
    return m_InitialRotationAngle;
  }

  virtual float GetInitialRotationStep() const
  {
    return m_InitialRotationStep;
  }

  virtual unsigned int GetNumRotationSteps() const
  {
    return m_NumRotationSteps;
  }

  virtual void SetNumDataSets(unsigned int x)
  {
    m_NumDataSets = x;
  }

  virtual void SetSearchboxDims(unsigned int x)
  {
    m_SearchboxDims = x;
  }

  virtual void SetResolutionUnits(float x)
  {
    m_ResolutionUnits = x;
  }

  virtual void SetInitialRotationAngle(float x)
  {
    m_InitialRotationAngle = x;
  }

  virtual void SetInitialRotationStep(float x)
  {
    m_InitialRotationStep = x;
  }

  virtual void SetNumRotationSteps(unsigned int x)
  {
    m_NumRotationSteps = x;
  }

protected:
  unsigned int m_NumDataSets;
  unsigned int m_SearchboxDims;
  float        m_ResolutionUnits;
  ValMapType   m_Height;
  ValMapType   m_Radius;
  float        m_InitialRotationAngle;
  float        m_InitialRotationStep;
  unsigned int m_NumRotationSteps;
};

#endif // _landmarksConstellationModelBase_h
