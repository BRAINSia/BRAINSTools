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
#ifndef __FSLDWIConverter_h
#define __FSLDWIConverter_h
#include "DWIConverter.h"
#include "StringContains.h"

/** specific converter for FSL nifti formatted files*/
class FSLDWIConverter : public DWIConverter
{
public:

  FSLDWIConverter( const DWIConverter::FileNamesContainer & inputFileNames,
  const std::string inputBValues, const std::string inputBVectors, const bool FSLFileFormatHorizontalBy3Rows);

  virtual ~FSLDWIConverter() {}

  virtual void AddFlagsToDictionary() ITK_OVERRIDE;

  /**
   * @brief FSL datasets are always in  normal sequential volume arrangement.
   */
   virtual void LoadFromDisk() ITK_OVERRIDE;

   /**
    * @brief  find the bvalues and gradient vectors
    */
  void ExtractDWIData() ITK_OVERRIDE;

  /**
   * @brief Return common fields.  Does nothing for FSL
   * @return empty map
   */
  virtual CommonDicomFieldMapType GetCommonDicomFieldsMap() const ITK_OVERRIDE
  {
    return CommonDicomFieldMapType();
  }

private:
  std::string m_inputBValues;
  std::string m_inputBVectors;
};

#endif // __FSLDWIConverter_h
