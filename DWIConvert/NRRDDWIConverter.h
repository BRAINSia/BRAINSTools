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

  typedef itk::VectorImage<PixelValueType, 3> VectorVolumeType;

  NRRDDWIConverter( const DWIConverter::FileNamesContainer & inputFileNames, const bool FSLFileFormatHorizontalBy3Rows );

  virtual ~NRRDDWIConverter() {}

  virtual void AddFlagsToDictionary() ITK_OVERRIDE;

  /**
   * @brief FSL datasets are always in  normal sequential volume arrangement.
   */
  virtual void LoadFromDisk() ITK_OVERRIDE;

  /**
   * @brief  find the bvalues and gradient vectors
   */
  virtual void ExtractDWIData() ITK_OVERRIDE;

  /**
   * @brief Return common fields.  Does nothing for FSL
   * @return empty map
   */
  virtual CommonDicomFieldMapType GetCommonDicomFieldsMap() const ITK_OVERRIDE
  {
    return CommonDicomFieldMapType();
  }

private:
  Volume4DType::Pointer CreateVolume(VectorVolumeType::Pointer & inputVol);
  std::string m_inputBValues;
  std::string m_inputBVectors;
};


#endif //BRAINSTOOLS_NRRDDWICONVERTER_H
