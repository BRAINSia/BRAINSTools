//
// Created by Hui Xie on 12/19/16.
//

#include "DWIConverterFactory.h"

DWIConverterFactory::DWIConverterFactory(const std::string DicomDirectory,
                                         const bool UseBMatrixGradientDirections,
                                         const bool FSLFileFormatHorizontalBy3Rows,
                                         const double smallGradientThreshold)
        : m_DicomDirectory(DicomDirectory)
        , m_UseBMatrixGradientDirections(UseBMatrixGradientDirections)
        , m_FSLFileFormatHorizontalBy3Rows(FSLFileFormatHorizontalBy3Rows)
        , m_SmallGradientThreshold(smallGradientThreshold)
{
}

DWIConverterFactory::~DWIConverterFactory()
{
  for( std::vector<itk::DCMTKFileReader *>::iterator it = this->m_Headers.begin();
       it != this->m_Headers.end(); ++it )
  {
    delete (*it);
  }
}

bool DWIConverterFactory::isNIIorNrrd( const std::string & filename )
{
  const size_t NUMEXT=4;
  const char * extensions [NUMEXT] = { ".nii", ".nii.gz", ".nhdr", ".nrrd"};
  for( size_t i=0; i< NUMEXT; ++i)
  {
    if( filename.find(extensions[i]) != std::string::npos )
    {
      return true; //Return true if one of the valid extensions found
    }
  }
  return false;
}

DWIConverter* DWIConverterFactory::New()
{

  // Directory of DICOM slices?
  if(itksys::SystemTools::FileIsDirectory(m_DicomDirectory.c_str()))
  {
    DWIDICOMConverterBase::InputNamesGeneratorType::Pointer inputNames =
            DWIDICOMConverterBase::InputNamesGeneratorType::New();
    inputNames->SetUseSeriesDetails( true);
    inputNames->SetLoadSequences( true );
    inputNames->SetLoadPrivateTags( true );
    inputNames->SetInputDirectory(m_DicomDirectory);
    m_InputFileNames = inputNames->GetInputFileNames();
  }
    // single file multiSlice Volume?
  else if(itksys::SystemTools::FileExists(m_DicomDirectory.c_str()))
  {
    m_InputFileNames.push_back(m_DicomDirectory);
  }

  DWIConverter *converter(ITK_NULLPTR);
  // nothing at all found?
  if(m_InputFileNames.size() < 1)
  {
    std::cerr << "Error: no DICOMfiles found in inputDirectory: " << m_DicomDirectory
              << std::endl;
    return ITK_NULLPTR;
  }
  else if( m_InputFileNames.size() == 1 &&  isNIIorNrrd( m_InputFileNames[0])) // FSL Reader or NRRD Reader
  {
    itkGenericExceptionMacro(<< "INVALID PATH, create FSLDWIConverter in main program" << std::endl);
    converter = new FSLDWIConverter(m_InputFileNames,"","", this->m_FSLFileFormatHorizontalBy3Rows);
  }
  else  // Assume multi file dicom file reading
  {
    m_Headers.resize(m_InputFileNames.size());
    int                                 headerCount = 0;
    for( unsigned i = 0; i < m_Headers.size(); ++i )
    {
      itk::DCMTKFileReader *curReader = new itk::DCMTKFileReader;
      curReader->SetFileName(m_InputFileNames[i]);
      try
      {
        curReader->LoadFile();
      }
      catch( ... )
      {
        std::cerr << "Error reading slice" << m_InputFileNames[i] << std::endl;
        delete curReader;
        curReader = ITK_NULLPTR;
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
          headerCount++;
        }
      }
    }
    // no headers found, nothing to do.
    if( headerCount == 0 )
    {
      std::cerr << "No pixel data in series " << m_DicomDirectory << std::endl;
      return ITK_NULLPTR;
    }
    m_InputFileNames.resize(0);
    //
    // clean the filename by traversing the header vector
    for(unsigned int i = 0; i < this->m_Headers.size(); ++i)
    {
      m_InputFileNames.push_back(m_Headers[i]->GetFileName());
    }
    try
    {
      m_Headers[0]->GetElementLO(0x0008, 0x0070, this->m_Vendor);
      strupper(this->m_Vendor);
    }
    catch(itk::ExceptionObject &excp)
    {
      std::cerr << "Can't get vendor name from DICOM file" << excp << std::endl;
      return ITK_NULLPTR;
    }

    if(StringContains(this->m_Vendor,"PHILIPS"))
    {
      converter = new PhilipsDWIConverter(m_Headers,m_InputFileNames,
                                          m_UseBMatrixGradientDirections, m_FSLFileFormatHorizontalBy3Rows);
    }
    else if(StringContains(this->m_Vendor,"SIEMENS"))
    {
      converter = new SiemensDWIConverter(m_Headers,m_InputFileNames,
                                          m_UseBMatrixGradientDirections,
                                          m_SmallGradientThreshold, m_FSLFileFormatHorizontalBy3Rows);
    }
    else if(StringContains(this->m_Vendor,"GE"))
    {
      converter = new GEDWIConverter(m_Headers,m_InputFileNames,
                                     m_UseBMatrixGradientDirections, m_FSLFileFormatHorizontalBy3Rows);
    }
    else if(StringContains(this->m_Vendor,"HITACHI"))
    {
      converter = new HitachiDWIConverter(m_Headers,m_InputFileNames,
                                          m_UseBMatrixGradientDirections, m_FSLFileFormatHorizontalBy3Rows);
    }
    else
    {
      // generic converter can't do anything except load a DICOM
      // directory
      converter = new GenericDWIConverter(m_InputFileNames, m_FSLFileFormatHorizontalBy3Rows);
      this->m_Vendor = "GENERIC";
    }
  }
  return converter;
}
std::string DWIConverterFactory::GetVendor() { return m_Vendor; }