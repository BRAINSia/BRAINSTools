//
// Created by Johnson, Hans J on 11/25/16.
//

#ifndef BRAINSTOOLS_NRRDDWICONVERTER_H
#define BRAINSTOOLS_NRRDDWICONVERTER_H

#include "DWIConverter.h"
#include "StringContains.h"

/** specific converter for FSL nifti formatted files*/
class NRRDDWIConverter : public DWIConverter
{
public:



  NRRDDWIConverter( const DWIConverter::FileNamesContainer & inputFileNames );

  ~NRRDDWIConverter() override {}

  void AddFlagsToDictionary() override;

  /**
   * @brief FSL datasets are always in  normal sequential volume arrangement.
   */
  void LoadFromDisk() override;

  /**
   * @brief  find the bvalues and gradient vectors
   */
  void ExtractDWIData() override;

  /**
   * @brief Return common fields.  Does nothing for FSL
   * @return empty map
   */
  CommonDicomFieldMapType GetCommonDicomFieldsMap() const override;

private:
  Volume4DType::Pointer CreateVolume(VectorVolumeType::Pointer & inputVol);
  std::string m_inputBValues;
  std::string m_inputBVectors;
};


#endif //BRAINSTOOLS_NRRDDWICONVERTER_H
