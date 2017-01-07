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
}


ScalarImage4DType::Pointer
NRRDDWIConverter::CreateVolume(VectorImage3DType::Pointer & vector3DVolume)
{
  VectorImage3DType::SizeType      inputSize = vector3DVolume->GetLargestPossibleRegion().GetSize();
  VectorImage3DType::SpacingType   inputSpacing = vector3DVolume->GetSpacing();
  VectorImage3DType::PointType     inputOrigin = vector3DVolume->GetOrigin();
  VectorImage3DType::DirectionType inputDirection = vector3DVolume->GetDirection();

  ScalarImage4DType::Pointer fourDVolume = ScalarImage4DType::New();
  ScalarImage4DType::SizeType      volSize;
  ScalarImage4DType::SpacingType   volSpacing;
  ScalarImage4DType::PointType     volOrigin;
  ScalarImage4DType::DirectionType volDirection;

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

  const ScalarImage4DType::IndexType::IndexValueType vecLength = vector3DVolume->GetNumberOfComponentsPerPixel();

  VectorImage3DType::IndexType vecIndex;
  ScalarImage4DType::IndexType       volIndex;
  // convert from vector image to 4D volume image
  for( volIndex[3] = 0; volIndex[3] < vecLength; ++volIndex[3] )
  {
    for( volIndex[2] = 0; volIndex[2] < static_cast<ScalarImage4DType::IndexType::IndexValueType>( inputSize[2] ); ++volIndex[2] )
    {
      vecIndex[2] = volIndex[2];
      for( volIndex[1] = 0; volIndex[1] < static_cast<ScalarImage4DType::IndexType::IndexValueType>( inputSize[1] ); ++volIndex[1] )
      {
        vecIndex[1] = volIndex[1];
        for( volIndex[0] = 0; volIndex[0] < static_cast<ScalarImage4DType::IndexType::IndexValueType>( inputSize[0] ); ++volIndex[0] )
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
  const std::string nrrdNRRDFile = m_InputFileNames[0];

  /*VectorImage3DType::Pointer vector3DVolume;
  if( ReadVectorVolume<VectorImage3DType>( vector3DVolume, nrrdNRRDFile, this->m_allowLossyConversion ) != EXIT_SUCCESS )
  {
    itkGenericExceptionMacro(<< "ERROR Reading NRRD File : " << nrrdNRRDFile << std::endl;);
  }

  //Conert vector 3D volume to 4DVolume
  ScalarImage4DType::Pointer                 fourDVolume = CreateVolume(vector3DVolume);
  this->m_Vector3DVolume = Convert4DVolumeTo3DVectorVolume(fourDVolume);*/

  // modified by Hui Xie Jan 6th, 2016
  if( ReadVectorVolume<VectorImage3DType>( m_Vector3DVolume, nrrdNRRDFile, this->m_allowLossyConversion ) != EXIT_SUCCESS )
  {
      itkGenericExceptionMacro(<< "ERROR Reading NRRD File : " << nrrdNRRDFile << std::endl;);
  }


}

void
NRRDDWIConverter::ExtractDWIData()
{
  RecoverMeasurementFrame<VectorImage3DType>(this->m_Vector3DVolume.GetPointer(), this->m_MeasurementFrame);
  RecoverBVectors<VectorImage3DType>(this->m_Vector3DVolume.GetPointer(), this->m_DiffusionVectors);
  RecoverBValues<VectorImage3DType>(this->m_Vector3DVolume.GetPointer(), this->m_DiffusionVectors, this->m_BValues);
}

DWIConverter::CommonDicomFieldMapType NRRDDWIConverter::GetCommonDicomFieldsMap() const
{
  return CommonDicomFieldMapType();
}