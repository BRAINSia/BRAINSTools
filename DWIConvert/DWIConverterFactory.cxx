//
// Created by Hui Xie on 12/19/16.
//

#include "DWIConverterFactory.h"

DWIConverterFactory::DWIConverterFactory(const std::string DicomDirectory,
                                         const bool UseBMatrixGradientDirections,
                                         const double smallGradientThreshold)
        : m_DicomDirectory(DicomDirectory)
        , m_UseBMatrixGradientDirections(UseBMatrixGradientDirections)
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

  DWIConverter *converter(nullptr);
  // nothing at all found?
  if(m_InputFileNames.size() < 1)
  {
    std::cerr << "Error: no DICOMfiles found in inputDirectory: " << m_DicomDirectory
              << std::endl;
    return nullptr;
  }

  // there is a logic error below --huixie
  else if( m_InputFileNames.size() == 1 &&  isNIIorNrrd( m_InputFileNames[0])) // FSL Reader or NRRD Reader
  {
    itkGenericExceptionMacro(<< "INVALID PATH, create FSLDWIConverter in main program" << std::endl);
    converter = new FSLDWIConverter(m_InputFileNames,"","");
  }
  else  // Assume multi file dicom file reading
  {

    // below code has logic error, it has empty hole in the tail of m_Headers when HasPixelData == false;
    /*
     *
     * m_Headers.resize(m_InputFileNames.size());
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
          headerCount++;
        }
      }
    }*/

    //modified by HuiXie
    m_Headers.clear();
    int  headerCount = 0;
    for( unsigned i = 0; i < m_InputFileNames.size(); ++i )
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
          m_Headers.push_back(curReader);
          headerCount++;
        }
      }
    }
    //end of modified by HuiXie


    // no headers found, nothing to do.
    if( headerCount == 0 )
    {
      std::cerr << "No pixel data in series " << m_DicomDirectory << std::endl;
      return nullptr;
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
      return nullptr;
    }

    if(StringContains(this->m_Vendor,"PHILIPS"))
    {
      converter = new PhilipsDWIConverter(m_Headers,m_InputFileNames,
                                          m_UseBMatrixGradientDirections);
    }
    else if(StringContains(this->m_Vendor,"SIEMENS"))
    {
      converter = new SiemensDWIConverter(m_Headers,m_InputFileNames,
                                          m_UseBMatrixGradientDirections,
                                          m_SmallGradientThreshold);
    }
    else if(StringContains(this->m_Vendor,"GE"))
    {
      converter = new GEDWIConverter(m_Headers,m_InputFileNames,
                                     m_UseBMatrixGradientDirections);
    }
    else if(StringContains(this->m_Vendor,"HITACHI"))
    {
      converter = new HitachiDWIConverter(m_Headers,m_InputFileNames,
                                          m_UseBMatrixGradientDirections);
    }
    else
    {
      // generic converter can't do anything except load a DICOM
      // directory
      converter = new GenericDWIConverter(m_InputFileNames);
      this->m_Vendor = "GENERIC";
    }
  }
  return converter;
}
std::string DWIConverterFactory::GetVendor() { return m_Vendor; }
