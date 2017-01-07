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
#ifndef __DWIConverter_h
#define __DWIConverter_h

#include <vector>
#include <algorithm>

#include <iostream>
#include <string>
#include <sstream>
#include <ctype.h>

#include "itkMacro.h"
#include "itkIntTypes.h"
#include "itkDCMTKSeriesFileNames.h"
#include "itkDCMTKImageIO.h"
#include "itkRawImageIO.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itksys/Directory.hxx"
#include "itksys/SystemTools.hxx"
#include "itksys/Base64.h"


#include "itkDCMTKFileReader.h"
#include "itkMatrix.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkRawImageIO.h"
#include "itkImage.h"
#include "itkDCMTKFileReader.h"
#include "itkDCMTKImageIO.h"
#include "itkNumberToString.h"

#include "dcmtk/oflog/helpers/loglog.h"

#include "DWIConvertUtils.h"
#include "DWIMetaDataDictionaryValidator.h"
#include "StringContains.h"

/** the DWIConverter is a base class for all scanner-specific
 *  converters.  It handles the tasks that are required for all
 *  scanners. In particular it loads the DICOM directory, and fills
 *  out various data fields needed by the DWIConvert program in order
 *  to write out NRRD and other files.
 */
class DWIConverter
{
public:
  /* The internal default format for DWIConverter is an itk::VectorImage<PixelValueType,3> */
  typedef std::vector< std::string > FileNamesContainer;

  typedef Vector3DType::SpacingType            SpacingType;
  //typedef itk::ImageFileReader<Volume3DUnwrappedType>   SingleFileReaderType;

  typedef itk::Vector<double, 3>                        PointType;

  typedef std::map<std::string,std::string> CommonDicomFieldMapType;
  DWIConverter( const FileNamesContainer &inputFileNames, const bool FSLFileFormatHorizontalBy3Rows );
  virtual ~DWIConverter();

  virtual void LoadFromDisk() = 0 ;



  /** extract dwi data -- vendor specific so must happen in subclass
   *  implementing this method.
   */
  virtual void ExtractDWIData() = 0;
  virtual CommonDicomFieldMapType GetCommonDicomFieldsMap() const =0;

  /** access methods for image data */
  const DWIMetaDataDictionaryValidator::GradientTableType &GetDiffusionVectors() const;
  /**
    * @brief NRRD file format stores a single BValue, and sets all the BVectors to scaled version that represents
    *                  the magnitude of the BValue offset.
    */
  void ConvertToSingleBValueScaledDiffusionVectors();
  /**
   * @brief FSL Format requires unit gradient directions and separate
   *        BValues for each gradient direction
   */
  void ConvertToMutipleBValuesUnitScaledBVectors();

 /**
  * @brief FSL orientation to allow for convenient display of images and
  *        conformance with conventions used by dcm2niix & fslview
  * @param toFSL prefers FSL's internal data format layout, if false, prefer Dicom natural data layout
  *         FSL [1 0 0; 0 -1 0; 0 0 1]    Dicom [1 0 0; 0 1 0; 0 0 1]
  * @return Returns a 4D image pointer properly formatted
  */
 Volume4DType::Pointer OrientForFSLConventions (const bool toFSL=true );

 const std::vector<double> &GetBValues() const;
 void SetBValues( const std::vector<double> & inBValues );
 double GetMaxBValue() const;

 Vector3DType::Pointer GetDiffusionVolume() const ;

 SpacingType GetSpacing() const;

 Vector3DType::PointType GetOrigin() const;
 void SetOrigin(Vector3DType::PointType origin);

 RotationMatrixType   GetLPSDirCos() const;

 RotationMatrixType GetMeasurementFrame() const;

 std::string GetNRRDSpaceDefinition() const;

 unsigned short GetRows() const;

 unsigned short GetCols() const;
 unsigned short GetSlices() const;


    /**
   * @brief Force overwriting the gradient directions by inserting values read from specified file
   * @param gradientVectorFile The file with gradients specified for overwriting
   */
  void ReadGradientInformation(const std::string& inputBValues, const std::string &inputBVectors, const std::string &inputVolumeNameTemplate);

  /**
   * @brief ConvertBVectorsToIdentityMeasurementFrame, Convert the values of the gradients to
   * use an identity measurement frame. This is required by FSL outputs.
   */
  void ConvertBVectorsToIdentityMeasurementFrame();

  std::string MakeFileComment(
            const std::string& version,
            bool useBMatrixGradientDirections,
            bool useIdentityMeaseurementFrame,
            double smallGradientThreshold,
            const std::string conversionMode) const;

  void ManualWriteNRRDFile(
            const std::string& outputVolumeHeaderName,
            const std::string commentstring) const;


/** the DICOM datasets are read as 3D volumes, but they need to be
 *  written as 4D volumes for image types other than NRRD.
 */
  void WriteFSLFormattedFileSet(const std::string& outputVolumeHeaderName,
                             const std::string outputBValues, const std::string outputBVectors, Volume4DType::Pointer img4D) const;

  void WriteFSLFormattedFileSet(const std::string& outputVolumeHeaderName,
                                                const std::string outputBValues, const std::string outputBVectors) const;


  /**
   * @brief Choose if we are going to allow for lossy conversion by typecasting
   * to the only internally supported format of short int
   * @param allowLossyConvertsion (true = automatically convert to short int)
   */
  void SetAllowLossyConversion(const bool newValue);




protected:
  double ComputeMaxBvalue(const std::vector<double> &bValues) const;
  size_t has_valid_nifti_extension( std::string outputVolumeHeaderName ) const;

  /** add vendor-specific flags; */
  virtual void AddFlagsToDictionary() = 0;

  /** the names of all the filenames, needed to use
   *  itk::ImageSeriesReader
   */
  const FileNamesContainer  m_InputFileNames;
  bool m_allowLossyConversion; // Allow type-cast conversion from float to short storage format


  /** double conversion instance, for optimal printing of numbers as  text */
  itk::NumberToString<double> m_DoubleConvert;
  bool       m_FSLFileFormatHorizontalBy3Rows; // Format of FSL files on disk

  //the default data model for all of DWIConvert
  Vector3DType::Pointer   m_Vector3DVolume;

  /** measurement from for gradients if different than patient
   *  reference frame.
   */
  RotationMatrixType   m_MeasurementFrame;
  /** list of B Values for each volume */
  std::vector<double>  m_BValues;
  /** list of gradient vectors */
  DWIMetaDataDictionaryValidator::GradientTableType  m_DiffusionVectors;
  // A map of common dicom fields to be propagated to image
  std::map<std::string,std::string> m_CommonDicomFieldsMap;

// this is always "left-posterior-superior" in all cases that we currently support
    const std::string           m_NRRDSpaceDefinition;

};

#endif // __DWIConverter_h
