//
// Created by Johnson, Hans J on 11/25/16.
//

#include "NRRDDWIConverter.h"


NRRDDWIConverter::NRRDDWIConverter( const DWIConverter::FileNamesContainer & inputFileNames,
  const bool FSLFileFormatHorizontalBy3Rows)
  : DWIConverter( inputFileNames, FSLFileFormatHorizontalBy3Rows )
{
}

void
NRRDDWIConverter::AddFlagsToDictionary()
{
  //TODO:  Move the QFORM/SFORM codes here
}


DWIConverter::Volume4DType::Pointer
NRRDDWIConverter::CreateVolume(VectorVolumeType::Pointer & vector3DVolume)
{
  VectorVolumeType::SizeType      inputSize = vector3DVolume->GetLargestPossibleRegion().GetSize();
  VectorVolumeType::SpacingType   inputSpacing = vector3DVolume->GetSpacing();
  VectorVolumeType::PointType     inputOrigin = vector3DVolume->GetOrigin();
  VectorVolumeType::DirectionType inputDirection = vector3DVolume->GetDirection();

  Volume4DType::Pointer fourDVolume = Volume4DType::New();
  Volume4DType::SizeType      volSize;
  Volume4DType::SpacingType   volSpacing;
  Volume4DType::PointType     volOrigin;
  Volume4DType::DirectionType volDirection;

  for( unsigned int i = 0; i < 3; ++i )
  {
    volSize[i] = inputSize[i];
    volSpacing[i] = inputSpacing[i];
    volOrigin[i] = inputOrigin[i];
    for( unsigned int j = 0; j < 3; ++j )
    {
      volDirection[i][j] = inputDirection[i][j];
    }
    volDirection[3][i] = 0.0;
    volDirection[i][3] = 0.0;
  }
  volDirection[3][3] = 1.0;
  volSpacing[3] = 1.0;
  volOrigin[3] = 0.0;
  volSize[3] = vector3DVolume->GetNumberOfComponentsPerPixel();

  fourDVolume->SetRegions(volSize);
  fourDVolume->SetOrigin(volOrigin);
  fourDVolume->SetSpacing(volSpacing);
  fourDVolume->SetDirection(volDirection);
  fourDVolume->Allocate();

  const Volume4DType::IndexType::IndexValueType vecLength = vector3DVolume->GetNumberOfComponentsPerPixel();

  VectorVolumeType::IndexType vecIndex;
  Volume4DType::IndexType       volIndex;
  // convert from vector image to 4D volume image
  for( volIndex[3] = 0; volIndex[3] < vecLength; ++volIndex[3] )
  {
    for( volIndex[2] = 0; volIndex[2] < static_cast<Volume4DType::IndexType::IndexValueType>( inputSize[2] ); ++volIndex[2] )
    {
      vecIndex[2] = volIndex[2];
      for( volIndex[1] = 0; volIndex[1] < static_cast<Volume4DType::IndexType::IndexValueType>( inputSize[1] ); ++volIndex[1] )
      {
        vecIndex[1] = volIndex[1];
        for( volIndex[0] = 0; volIndex[0] < static_cast<Volume4DType::IndexType::IndexValueType>( inputSize[0] ); ++volIndex[0] )
        {
          vecIndex[0] = volIndex[0];
          fourDVolume->SetPixel(volIndex, vector3DVolume->GetPixel(vecIndex)[volIndex[3]]);
        }
      }
    }
  }

  fourDVolume->SetMetaDataDictionary(vector3DVolume->GetMetaDataDictionary());
  return fourDVolume;
}

void
NRRDDWIConverter::LoadFromDisk()
{
  //HACK: TODO:
  const bool allowLossyConversion=false;

  const std::string nrrdNRRDFile = m_InputFileNames[0];

  VectorVolumeType::Pointer vector3DVolume;
  if( ReadVolume<VectorVolumeType>( vector3DVolume, nrrdNRRDFile, allowLossyConversion ) != EXIT_SUCCESS )
  {
    itkGenericExceptionMacro(<< "ERROR Reading NRRD File : " << nrrdNRRDFile << std::endl;);
  }

  //Conert vector 3D volume to 4DVolume
  Volume4DType::Pointer                 fourDVolume = CreateVolume(vector3DVolume);
  this->m_SlicesPerVolume = fourDVolume->GetLargestPossibleRegion().GetSize()[2];
  this->m_NVolume = fourDVolume->GetLargestPossibleRegion().GetSize()[3];
  this->m_NSlice = this->m_SlicesPerVolume * this->m_NVolume;
  this->m_Volume = FourDToThreeDImage(fourDVolume);
}

void
NRRDDWIConverter::ExtractDWIData()
{
  RecoverBVectors<Volume3DUnwrappedType>(this->m_Volume.GetPointer(), this->m_DiffusionVectors);
  RecoverBValues<Volume3DUnwrappedType>(this->m_Volume.GetPointer(), this->m_DiffusionVectors, this->m_BValues);
}
