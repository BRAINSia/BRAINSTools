#ifndef __DWIConverterFactory_h
#define __DWIConverterFactory_h

#include "itkImageSeriesReader.h"
#include "itkDCMTKFileReader.h"
#include "itksys/SystemTools.hxx"
#include "DWIConverter.h"
#include "PhilipsDWIConverter.h"
#include "GEDWIConverter.h"
#include "SiemensDWIConverter.h"
#include "HitachiDWIConverter.h"
#include "GenericDWIConverter.h"

/** 'factory' object. Does initial scanning of the DICOM directory and
* then it can instantiate the correct DWIConverter class
*/
class DWIConverterFactory
{
public:
  DWIConverterFactory(const std::string DicomDirectory,
                      bool UseBMatrixGradientDirections,
                      double smallGradientThreshold) : m_DicomDirectory(DicomDirectory),
                                                       m_UseBMatrixGradientDirections(UseBMatrixGradientDirections),
                                                       m_SmallGradientThreshold(smallGradientThreshold)
    {
    }
  ~DWIConverterFactory()
    {
      for( std::vector<itk::DCMTKFileReader *>::iterator it = this->m_Headers.begin();
           it != this->m_Headers.end(); ++it )
        {
        delete (*it);
        }
    }
  DWIConverter *New()
    {
      // Directory of DICOM slices?
      if(itksys::SystemTools::FileIsDirectory(m_DicomDirectory.c_str()))
        {
        DWIConverter::InputNamesGeneratorType::Pointer inputNames =
          DWIConverter::InputNamesGeneratorType::New();
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
      // nothing at all found?
      if(m_InputFileNames.size() < 1)
        {
        std::cerr << "Error: no DICOMfiles found in inputDirectory: " << m_DicomDirectory
                  << std::endl;
        return 0;
        }

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
          curReader = 0;
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
        return 0;
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
        return 0;
        }
      DWIConverter *converter(0);
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
        converter = new GenericDWIConverter(m_Headers,m_InputFileNames,
                                            m_UseBMatrixGradientDirections);
        }
      return converter;
    }
  std::string GetVendor() { return m_Vendor; }
private:
  std::string m_DicomDirectory;
  std::string m_Vendor;
  bool        m_UseBMatrixGradientDirections;
  double      m_SmallGradientThreshold;

  DWIConverter::DCMTKFileVector m_Headers;
  DWIConverter::FileNamesContainer m_InputFileNames;

};

#endif // __DWIConverterFactory_h
