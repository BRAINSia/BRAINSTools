//
// Created by Johnson, Hans J on 11/25/16.
//

#include "FSLDWIConverter.h"

FSLDWIConverter::FSLDWIConverter(const DWIConverter::FileNamesContainer & inputFileNames,
                                 const std::string                        inputBValues,
                                 const std::string                        inputBVectors)
  : DWIConverter(inputFileNames)
  , m_inputBValues(inputBValues)
  , m_inputBVectors(inputBVectors)
{}

FSLDWIConverter::CommonDicomFieldMapType
FSLDWIConverter::GetCommonDicomFieldsMap() const
{
  return CommonDicomFieldMapType();
}

void
FSLDWIConverter::AddFlagsToDictionary()
{}

void
FSLDWIConverter::LoadFromDisk()
{
  const std::string fslNIFTIFile = m_InputFileNames[0];

  Volume4DType::Pointer inputVol;

  // string to use as template if no bval or bvec filename is given.
  if (ReadScalarVolume<Volume4DType>(inputVol, fslNIFTIFile, this->m_allowLossyConversion) != EXIT_SUCCESS)
    throw;
  // Reorient from FSL standard format to ITK/Dicom standard format
  this->m_SlicesPerVolume = inputVol->GetLargestPossibleRegion().GetSize()[2];
  this->m_NVolume = inputVol->GetLargestPossibleRegion().GetSize()[3];
  this->m_NSlice = this->m_SlicesPerVolume * this->m_NVolume;
  this->m_Volume = FourDToThreeDImage(inputVol);
}

void
FSLDWIConverter::ExtractDWIData()
{
  const std::string fslNIFTIFile = m_InputFileNames[0];
  this->ReadGradientInformation(m_inputBValues, m_inputBVectors, fslNIFTIFile);
  this->OrientForFSLConventions(false); // Orient for Dicom data layout
}
