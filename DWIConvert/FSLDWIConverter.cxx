//
// Created by Johnson, Hans J on 11/25/16.
//

#include "FSLDWIConverter.h"

FSLDWIConverter::FSLDWIConverter( const DWIConverter::FileNamesContainer & inputFileNames,
  const std::string inputBValues,
  const std::string inputBVectors, const bool FSLFileFormatHorizontalBy3Rows)
  : DWIConverter( inputFileNames, FSLFileFormatHorizontalBy3Rows )
  , m_inputBValues(inputBValues)
  , m_inputBVectors(inputBVectors)
{
}

void
FSLDWIConverter::AddFlagsToDictionary()
{
}

void
FSLDWIConverter::LoadFromDisk()
{
  const bool allowLossyConversion=false;

  const std::string fslNIFTIFile = m_InputFileNames[0];

  Volume4DType::Pointer inputVol;

  // string to use as template if no bval or bvec filename is given.
  ReadVolume<Volume4DType>(inputVol, fslNIFTIFile, allowLossyConversion);
  // Reorient from FSL standard format to ITK/Dicom standard format
  inputVol= DicomToFSLOrientationImageConverter(inputVol);
  this->m_SlicesPerVolume = inputVol->GetLargestPossibleRegion().GetSize()[2];
  this->m_NVolume = inputVol->GetLargestPossibleRegion().GetSize()[3];
  this->m_NSlice = this->m_SlicesPerVolume * this->m_NVolume;
  this->m_Volume = FourDToThreeDImage(inputVol);
}

void
FSLDWIConverter::ExtractDWIData()
{
  const std::string fslNIFTIFile = m_InputFileNames[0];
  this->ReadGradientInformation(m_inputBValues,m_inputBVectors,fslNIFTIFile);
  // Reorient from FSL standard format to ITK/Dicom standard format
  this->m_DiffusionVectors = DicomToFSLOrientationGradientTableConverter(m_DiffusionVectors);
}
