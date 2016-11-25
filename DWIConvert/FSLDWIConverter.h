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
  const std::string inputBValues,
  const std::string inputBVectors, const bool FSLFileFormatHorizontalBy3Rows)
    : DWIConverter( inputFileNames, FSLFileFormatHorizontalBy3Rows )
    , m_inputBValues(inputBValues)
    , m_inputBVectors(inputBVectors)
  {
  }
  virtual ~FSLDWIConverter() {}

  virtual void AddFlagsToDictionary() ITK_OVERRIDE
  {
    //TODO:  Move the QFORM/SFORM codes here
  }

  /**
   * @brief FSL datasets are always in  normal sequential volume arrangement.
   */
   virtual void LoadFromDisk() ITK_OVERRIDE
    {
      //HACK: TODO:
      const bool allowLossyConversion=false;

      const std::string fslNIFTIFile = m_InputFileNames[0];

      Volume4DType::Pointer inputVol;

      // string to use as template if no bval or bvec filename is given.
      ReadVolume<Volume4DType>(inputVol, fslNIFTIFile, allowLossyConversion);
      this->m_SlicesPerVolume = inputVol->GetLargestPossibleRegion().GetSize()[2];
      this->m_NVolume = inputVol->GetLargestPossibleRegion().GetSize()[3];
      this->m_NSlice = this->m_SlicesPerVolume * this->m_NVolume;
      this->m_Volume = FourDToThreeDImage(inputVol);

    }

   /**
    * @brief  find the bvalues and gradient vectors
    */
  void ExtractDWIData() ITK_OVERRIDE
    {
      const std::string fslNIFTIFile = m_InputFileNames[0];
      this->ReadGradientInformation(m_inputBValues,m_inputBVectors,fslNIFTIFile);
    }
private:
  std::string m_inputBValues;
  std::string m_inputBVectors;
};

#endif // __FSLDWIConverter_h
