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
#include "NRRDDWIConverter.h"

/** 'factory' object. Does initial scanning of the DICOM directory and
* then it can instantiate the correct DWIConverter class
*/
class DWIConverterFactory
{
public:
  DWIConverterFactory(const std::string DicomDirectory,
                        const bool UseBMatrixGradientDirections,
                        const bool FSLFileFormatHorizontalBy3Rows,
                        const double smallGradientThreshold);

  ~DWIConverterFactory();

  static bool isNIIorNrrd( const std::string & filename );
  DWIConverter* New();
  std::string GetVendor();

private:
  std::string m_DicomDirectory;
  std::string m_Vendor;
  bool        m_UseBMatrixGradientDirections;
  bool        m_FSLFileFormatHorizontalBy3Rows;
  double      m_SmallGradientThreshold;

  DWIDICOMConverterBase::DCMTKFileVector m_Headers;
  DWIConverter::FileNamesContainer m_InputFileNames;

};

#endif // __DWIConverterFactory_h
