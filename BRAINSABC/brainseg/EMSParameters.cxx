#include "EMSParameters.h"

#include "itksys/SystemTools.hxx"
#include "DoubleToString.h"

EMSParameters
::EMSParameters() :
  m_Suffix(""),
  m_AtlasDirectory(""),
  m_AtlasOrientation("RAI"),
  m_OutputDirectory(""),
  m_OutputFormat("Meta")

{
  m_DoAtlasWarp = true;

  m_Images.clear();
  m_ImageOrientations.clear();

  m_FilterMethod = "CurvatureFlow";
  m_FilterIterations = 1;
  m_FilterTimeStep = 0.01;

  m_MaxBiasDegree = 4;

  m_AtlasWarpGridX = 5;
  m_AtlasWarpGridY = 5;
  m_AtlasWarpGridZ = 5;

  m_Prior1 = 1.0;
  m_Prior2 = 1.0;
  m_Prior3 = 1.0;
  m_Prior4 = 1.0;

  m_AtlasLinearMapType = "BSpline";
  m_ImageLinearMapType = "Rigid";
}

EMSParameters
::~EMSParameters()
{
}

void
EMSParameters
::AddImage(std::string s, std::string orient)
{
  m_Images.push_back(s);
  m_ImageOrientations.push_back(orient);
}

void
EMSParameters
::ClearImages()
{
  m_Images.clear();
  m_ImageOrientations.clear();
}

bool
EMSParameters
::CheckValues() const
{
  if( m_Suffix.length() == 0 )
    {
    std::cerr << "Invalid Suffix length" << std::endl;
    return false;
    }

  if( m_AtlasDirectory.length() == 0 )
    {
    std::cerr << "Invalid Atlas Directory" << std::endl;
    return false;
    }

  if( m_OutputDirectory.length() == 0 )
    {
    std::cerr << "Invalid OutputDirectory" << std::endl;
    return false;
    }

  bool validFormat = false;
  if( itksys::SystemTools::Strucmp(m_OutputFormat.c_str(), "Analyze") == 0 )
    {
    validFormat = true;
    }
  if( itksys::SystemTools::Strucmp(m_OutputFormat.c_str(), "GIPL") == 0 )
    {
    validFormat = true;
    }
  if( itksys::SystemTools::Strucmp(m_OutputFormat.c_str(), "Nrrd") == 0 )
    {
    validFormat = true;
    }
  if( itksys::SystemTools::Strucmp(m_OutputFormat.c_str(), "Meta") == 0 )
    {
    validFormat = true;
    }
  if( itksys::SystemTools::Strucmp(m_OutputFormat.c_str(), "NIFTI") == 0 )
    {
    validFormat = true;
    }

  if( !validFormat )
    {
    std::cerr << "No valid image format defined." << std::endl;
    return false;
    }

  if( m_Images.empty() )
    {
    std::cerr << "Images not defined." << std::endl;
    return false;
    }

  return true;
}

void
EMSParameters
::PrintSelf(std::ostream & os, itk::Indent) const
{
  DoubleToString doubleConvert;

  os << "Suffix = " << m_Suffix << std::endl;
  os << "Atlas directory = " << m_AtlasDirectory << std::endl;
  os << "Atlas orientation = " << m_AtlasOrientation << std::endl;
  os << "Output directory = " << m_OutputDirectory << std::endl;
  os << "Output format = " << m_OutputFormat << std::endl;
  os << "Images:" << std::endl;
  for( unsigned int k = 0; k < m_Images.size(); k++ )
    {
    os << "  " << m_Images[k] << " --- " << m_ImageOrientations[k] << std::endl;
    }
  os << "Filter iterations = " << m_FilterIterations << std::endl;
  os << "Filter time step = " << doubleConvert(m_FilterTimeStep) << std::endl;
  os << "Max bias degree = " << m_MaxBiasDegree << std::endl;
  os << "Prior 1 = " << m_Prior1 << std::endl;
  os << "Prior 2 = " << m_Prior2 << std::endl;
  os << "Prior 3 = " << m_Prior3 << std::endl;
  os << "Prior 4 = " << m_Prior4 << std::endl;
  if( m_DoAtlasWarp )
    {
    os << "Atlas warping, grid = "
       << m_AtlasWarpGridX << "x"
       << m_AtlasWarpGridY << "x"
       << m_AtlasWarpGridZ << std::endl;
    }
  else
    {
    os << "No atlas warping..." << std::endl;
    }
}
