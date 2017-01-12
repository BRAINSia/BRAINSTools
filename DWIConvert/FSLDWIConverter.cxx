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

FSLDWIConverter::CommonDicomFieldMapType FSLDWIConverter::GetCommonDicomFieldsMap() const
{
  return CommonDicomFieldMapType();
}

void
FSLDWIConverter::AddFlagsToDictionary()
{
}

void
FSLDWIConverter::LoadFromDisk()
{
  const std::string fslNIFTIFile = m_InputFileNames[0];

  ScalarImage4DType::Pointer inputVol;

  // string to use as template if no bval or bvec filename is given.
  ReadScalarVolume<ScalarImage4DType>(inputVol, fslNIFTIFile, this->m_allowLossyConversion);
  // Reorient from FSL standard format to ITK/Dicom standard format

  //this->m_vectorImage3D = convertScalarImage4DToVectorImage3D(inputVol);
  // Reorient from FSL standard format to ITK/Dicom standard format
  this->m_SlicesPerVolume = inputVol->GetLargestPossibleRegion().GetSize()[2];
  this->m_NVolume = inputVol->GetLargestPossibleRegion().GetSize()[3];
  this->m_NSlice = this->m_SlicesPerVolume * this->m_NVolume;
  this->m_scalarImage3D = FourDToThreeDUnwrappedImage(inputVol);

}

void
FSLDWIConverter::ExtractDWIData()
{
  const std::string fslNIFTIFile = m_InputFileNames[0];
  this->ReadGradientInformation(m_inputBValues,m_inputBVectors,fslNIFTIFile);
  this->OrientForFSLConventions(false); //Orient for Dicom data layout
}
