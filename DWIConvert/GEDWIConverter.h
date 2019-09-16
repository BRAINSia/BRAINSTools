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
#ifndef __GEDWIConverter_h
#define __GEDWIConverter_h
#include "DWIDICOMConverterBase.h"

/** Specific converter for GE Scanners */
class GEDWIConverter : public DWIDICOMConverterBase
{
public:
  GEDWIConverter(DWIDICOMConverterBase::DCMTKFileVector & allHeaders,
                 DWIConverter::FileNamesContainer &       inputFileNames,
                 const bool                               useBMatrixGradientDirections);

  ~GEDWIConverter() override;
  void
  LoadDicomDirectory() override;

  void
  ExtractDWIData() override;

protected:
  void
  AddFlagsToDictionary() override;
};

#endif // __GEDWIConverter_h
