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
#ifndef __GenericDWIConverter_h
#define __GenericDWIConverter_h
#include "DWIConverter.h"

class GenericDWIConverter : public DWIConverter
{
public:
  GenericDWIConverter( DWIConverter::FileNamesContainer &inputFileNames , const bool FSLFileFormatHorizontalBy3Rows);

  virtual void LoadFromDisk() override;


  virtual ~GenericDWIConverter();

protected:
  virtual void ExtractDWIData() override;
  virtual void AddFlagsToDictionary() override;

  /**
 * @brief Return common fields.  Does nothing for FSL
 * @return empty map
 */
  virtual CommonDicomFieldMapType GetCommonDicomFieldsMap() const override;

};

#endif // __GenericDWIConverter_h
