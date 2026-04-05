//
// Created by Johnson, Hans J on 11/24/16.
//
#ifndef BRAINSTOOLS_DWIDICOMCONVERTERBASE_H
#define BRAINSTOOLS_DWIDICOMCONVERTERBASE_H

#include <vector>
#include <iostream>
#include "DWIConverter.h"
#include "itkDCMTKSeriesFileNames.h"
#include "itkMacro.h"
#include "itkDCMTKImageIO.h"
#include "itkImage.h"
#include "itkDCMTKFileReader.h"
#include "itkNumberToString.h"
#include "DWIConvertUtils.h"

class DWIDICOMConverterBase : public DWIConverter
{
public:
  using InputNamesGeneratorType = itk::DCMTKSeriesFileNames;
  using DCMTKFileVector = std::vector<itk::DCMTKFileReader *>;

  DWIDICOMConverterBase(DCMTKFileVector            allHeaders,
                        const FileNamesContainer & inputFileNames,
                        const bool                 useBMatrixGradientDirections);

  /**
   * @brief Return common fields.  Does nothing for FSL
   * @return empty map
   */
  CommonDicomFieldMapType
  GetCommonDicomFieldsMap() const override;

  void
  LoadFromDisk() override;

  virtual void
  LoadDicomDirectory();
  double
  readThicknessFromDicom() const;
  int
  getDicomSpacing(double * const spacing) const;

protected:
  enum VRType
  {
    DCM_CS,
    DCM_LO,
    DCM_SH,
    DCM_DS,
  };

  /**
   * @brief Create a dictionary of dicom extracted information about the scans
   *   using dcmtk dcm2xml tool, identify the information desired to be kept
   *   <element tag="0018,1314" vr="DS" vm="1" len="2" name="FlipAngle">90</element>
   * @param dcm_primary_name "0018" in example above
   * @param dcm_seconary_name "1314" in example above
   * @param dcm_human_readable_name "FlipAngle" in example above
   * @param vr "DCM_DS" for enumeration as indicated by vr in example above
   */
  void
  _addToStringDictionary(const std::string & dcm_primary_name,
                         const std::string & dcm_seconary_name,
                         const std::string & dcm_human_readable_name,
                         const enum VRType   vr);

  /** the SliceOrderIS flag can be computed (as above) but if it's
   *  invariant, the derived classes can just set the flag. This method
   *  fixes up the VolumeDirectionCos after the flag is set.
   */
  void
  SetDirectionsFromSliceOrder();


  /* given a sequence of dicom files where all the slices for location
   * 0 are folled by all the slices for location 1, etc. This method
   * transforms it into a sequence of volumes
   */
  void
  DeInterleaveVolume();
  /* determine if slice order is inferior to superior */
  void
  DetermineSliceOrderIS();

  /** force use of the BMatrix to compute gradients in Siemens data instead of
   *  the reported gradients. which are in many cases bogus.
   */
  const bool m_UseBMatrixGradientDirections;
  /** one file reader per DICOM file in dataset */
  const DCMTKFileVector m_Headers;

  /** matrix with just spacing information, used a couple places */
  /** the current dataset is represented in a single file */
  bool m_MultiSliceVolume;
  /** slice order is inferior/superior? */
  bool m_SliceOrderIS;


  /** track if images is interleaved */
  bool m_IsInterleaved;

  /**
   * @brief Try to extract DWI gradient data using the DICOM Supplement 49
   *   standard tags (shared by Hitachi scanners and other vendors that follow
   *   the standard): b-value from (0018,9087) and gradient direction from
   *   (0018,9089) inside the SharedFunctionalGroupsSequence (5200,9229).
   *
   *   This method fills m_BValues and m_DiffusionVectors and returns true on
   *   success.  Returns false (without modifying the member vectors) if the
   *   tags are absent or cannot be read, allowing the caller to fall back to
   *   vendor-specific tag parsing.
   *
   *   Addresses issue #294: vendor converters can call this as a fallback so
   *   that scanners emitting only the standard Supplement 49 tags are
   *   supported without a dedicated converter class.
   */
  bool
  TryExtractSupp49DWIData();
};

#endif // BRAINSTOOLS_DWIDICOMCONVERTERBASE_H
