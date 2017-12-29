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
#ifndef __PhilipsDWIConverter_h
#define __PhilipsDWIConverter_h

#include "DWIDICOMConverterBase.h"
#include "itkExtractImageFilter.h"

/** specific converter for Philips scanners */
class PhilipsDWIConverter : public DWIDICOMConverterBase
{
public:
  PhilipsDWIConverter(DWIDICOMConverterBase::DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      const bool useBMatrixGradientDirections,
                      const bool FSLFileFormatHorizontalBy3Rows) ;
  ~PhilipsDWIConverter() override;

  void LoadDicomDirectory() override;
  void ExtractDWIData() override;
protected:
  void AddFlagsToDictionary() override;
  /** # of trailing images to ignore */
  unsigned int        m_NTrailingImagesToIgnore;
};

#endif // __PhilipsDWIConverter_h
