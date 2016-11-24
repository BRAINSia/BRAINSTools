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
  FSLDWIConverter(DWIConverter::DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      bool useBMatrixGradientDirections,
                      double smallGradientThreshold) : DWIConverter(allHeaders,inputFileNames,
                                                                    useBMatrixGradientDirections),
                                                       m_SmallGradientThreshold(smallGradientThreshold)
    {
    }
  virtual ~FSLDWIConverter() {}

  /**
   * @brief FSL datasets are always in  normal sequential volume arrangement.
   */
  virtual void LoadDicomDirectory() ITK_OVERRIDE
    {
    }

   /**
    * @brief  find the bvalues and gradient vectors */
    */
  void ExtractDWIData() ITK_OVERRIDE
    {
    }
private:
};

#endif // __FSLDWIConverter_h
