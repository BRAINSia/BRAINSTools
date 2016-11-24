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
#include "FSLDWIConverter.h"

/** 'factory' object. Does initial scanning of the DICOM directory and
* then it can instantiate the correct DWIConverter class
*/
class DWIConverterFactory
{
public:
  DWIConverterFactory(const std::string DicomDirectory,
                      bool UseBMatrixGradientDirections,
                      double smallGradientThreshold,
                      bool useIdentityMeaseurementFrame)
    : m_DicomDirectory(DicomDirectory)
    , m_UseBMatrixGradientDirections(UseBMatrixGradientDirections)
    , m_SmallGradientThreshold(smallGradientThreshold)
    , m_useIdentityMeaseurementFrame(useIdentityMeaseurementFrame)
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

  static bool isNIIorNIFTI( const std::string & filename )
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

  DWIConverter *New()
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
      else if( m_InputFileNames.size() == 1 &&  isNIIorNIFTI( m_InputFileNames[0])) // FSL Reader or NRRD Reader
      {
        converter = new FSLDWIConverter(m_InputFileNames);
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
        converter->SetUseIdentityMeaseurementFrame(this->m_useIdentityMeaseurementFrame);
      }
      return converter;
    }
  std::string GetVendor() { return m_Vendor; }
private:
  std::string m_DicomDirectory;
  std::string m_Vendor;
  bool        m_UseBMatrixGradientDirections;
  double      m_SmallGradientThreshold;
  bool        m_useIdentityMeaseurementFrame;

  DWIDICOMConverterBase::DCMTKFileVector m_Headers;
  DWIConverter::FileNamesContainer m_InputFileNames;

};

#endif // __DWIConverterFactory_h
