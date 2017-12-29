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

#ifndef __HitachiDWIConverter_h
#define __HitachiDWIConverter_h
#include "DWIConverter.h"
#include "DWIDICOMConverterBase.h"

/** specific converter for Hitachi scanners */
class HitachiDWIConverter : public DWIDICOMConverterBase
{
public:
  HitachiDWIConverter(DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      const bool useBMatrixGradientDirections,
                      const bool FSLFileFormatHorizontalBy3Rows);

  virtual ~HitachiDWIConverter();
  /* load dicom directory -- no postprocessing necessary after letting
   * superclass do its thing.
   */
  virtual void LoadDicomDirectory() override;
  /** extract gradient vectors.
   *  Hitachi apparently supports the Supplement 49 definition
   *  for Diffusion data.-- see page 94 of the Supplement 49 document:
   *  ftp://medical.nema.org/medical/dicom/final/sup49_ft.pdf
   */
  void ExtractDWIData() override;

protected:
  virtual void AddFlagsToDictionary() override;

};

#endif // __HitachiDWIConverter_h
