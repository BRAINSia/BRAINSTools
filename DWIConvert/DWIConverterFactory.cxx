//
// Created by Hui Xie on 12/19/16.
//

#include <utility>

#include "DWIConverterFactory.h"
#include <set>

DWIConverterFactory::DWIConverterFactory(std::string  DicomDirectory,
                                         const bool   UseBMatrixGradientDirections,
                                         const double smallGradientThreshold)
  : m_DicomDirectory(std::move(DicomDirectory))
  , m_UseBMatrixGradientDirections(UseBMatrixGradientDirections)
  , m_SmallGradientThreshold(smallGradientThreshold)
{}

bool
DWIConverterFactory::isNIIorNrrd(const std::string & filename)
{
  constexpr size_t NUMEXT = 4;
  const char *     extensions[NUMEXT] = { ".nii", ".nii.gz", ".nhdr", ".nrrd" };
  for (auto & extension : extensions)
  {
    if (filename.find(extension) != std::string::npos)
    {
      return true; // Return true if one of the valid extensions found
    }
  }
  return false;
}

std::unique_ptr<DWIConverter>
DWIConverterFactory::New()
{

  // Directory of DICOM slices?
  if (itksys::SystemTools::FileIsDirectory(m_DicomDirectory.c_str()))
  {
    DWIDICOMConverterBase::InputNamesGeneratorType::Pointer inputNames =
      DWIDICOMConverterBase::InputNamesGeneratorType::New();
    inputNames->SetUseSeriesDetails(true);
    inputNames->SetLoadSequences(true);
    inputNames->SetLoadPrivateTags(true);
    inputNames->SetInputDirectory(m_DicomDirectory);
    // Series details split a DWI by Sequence Name; keep every group of the first Series Instance UID.
    const auto            seriesUID = [](const std::string & id) { return id.substr(0, id.find('_')); };
    const auto            seriesIdentifiers = inputNames->GetSeriesUIDs();
    std::set<std::string> merged;
    for (const auto & id : seriesIdentifiers)
    {
      if (seriesUID(id) == seriesUID(seriesIdentifiers.front()) && merged.insert(id).second)
      {
        const auto & files = inputNames->GetFileNames(id);
        m_InputFileNames.insert(m_InputFileNames.end(), files.begin(), files.end());
      }
    }
  }
  // single file multiSlice Volume?
  else if (itksys::SystemTools::FileExists(m_DicomDirectory.c_str()))
  {
    m_InputFileNames.push_back(m_DicomDirectory);
  }

  std::unique_ptr<DWIConverter> converter;
  // nothing at all found?
  if (m_InputFileNames.empty())
  {
    std::cerr << "Error: no DICOMfiles found in inputDirectory: " << m_DicomDirectory << std::endl;
    return nullptr;
  }

  // there is a logic error below --huixie
  else if (m_InputFileNames.size() == 1 && isNIIorNrrd(m_InputFileNames[0])) // FSL Reader or NRRD Reader
  {
    itkGenericExceptionMacro(<< "INVALID PATH, create FSLDWIConverter in main program" << std::endl);
  }
  else // Assume multi file dicom file reading
  {

    // below code has logic error, it has empty hole in the tail of m_Headers when HasPixelData == false;
    /*
     *
     * m_Headers.resize(m_InputFileNames.size());
    int                                 headerCount = 0;
    for( unsigned i = 0; i < m_Headers.size(); ++i )
    {
      brains::DCMTKFileReader *curReader = new brains::DCMTKFileReader;
      curReader->SetFileName(m_InputFileNames[i]);
      try
      {
        curReader->LoadFile();
      }
      catch (const ... )
      {
        std::cerr << "Error reading slice" << m_InputFileNames[i] << std::endl;
        delete curReader;
        curReader = nullptr;
      }
      // check for pixel data.
      if(curReader)
      {
        if(!curReader->HasPixelData() )
        {
          delete curReader;
        }
        else
        {
          this->m_Headers[headerCount] = curReader;
          ++headerCount;
        }
      }
    }*/

    // modified by HuiXie
    m_Headers.clear();
    int headerCount = 0;
    for (auto & m_InputFileName : m_InputFileNames)
    {
      auto curReader = std::make_unique<brains::DCMTKFileReader>();
      curReader->SetFileName(m_InputFileName);
      try
      {
        curReader->LoadFile();
      }
      catch (...)
      {
        std::cerr << "Error reading slice" << m_InputFileName << std::endl;
        continue;
      }
      if (curReader->HasPixelData())
      {
        m_Headers.push_back(std::move(curReader));
        ++headerCount;
      }
    }
    // end of modified by HuiXie


    // no headers found, nothing to do.
    if (headerCount == 0)
    {
      std::cerr << "No pixel data in series " << m_DicomDirectory << std::endl;
      return nullptr;
    }
    DWIDICOMConverterBase::OrderSlicesVolumeMajor(m_Headers);
    m_InputFileNames.resize(0);
    //
    // clean the filename by traversing the header vector
    for (unsigned int i = 0; i < this->m_Headers.size(); ++i)
    {
      m_InputFileNames.push_back(m_Headers[i]->GetFileName());
    }
    try
    {
      m_Headers[0]->GetElementLO(0x0008, 0x0070, this->m_Vendor);
      strupper(this->m_Vendor);
    }
    catch (const itk::ExceptionObject & excp)
    {
      std::cerr << "Can't get vendor name from DICOM file" << excp << std::endl;
      return nullptr;
    }

    if (StringContains(this->m_Vendor, "PHILIPS"))
    {
      converter =
        std::make_unique<PhilipsDWIConverter>(std::move(m_Headers), m_InputFileNames, m_UseBMatrixGradientDirections);
    }
    else if (StringContains(this->m_Vendor, "SIEMENS"))
    {
      converter = std::make_unique<SiemensDWIConverter>(
        std::move(m_Headers), m_InputFileNames, m_UseBMatrixGradientDirections, m_SmallGradientThreshold);
    }
    else if (StringContains(this->m_Vendor, "GE"))
    {
      converter =
        std::make_unique<GEDWIConverter>(std::move(m_Headers), m_InputFileNames, m_UseBMatrixGradientDirections);
    }
    else if (StringContains(this->m_Vendor, "HITACHI"))
    {
      converter =
        std::make_unique<HitachiDWIConverter>(std::move(m_Headers), m_InputFileNames, m_UseBMatrixGradientDirections);
    }
    else
    {
      // generic converter can't do anything except load a DICOM
      // directory
      converter = std::make_unique<GenericDWIConverter>(m_InputFileNames);
      this->m_Vendor = "GENERIC";
    }
  }
  return converter;
}
std::string
DWIConverterFactory::GetVendor()
{
  return m_Vendor;
}
