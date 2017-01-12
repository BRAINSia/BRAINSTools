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
#include "DWIConvertUtils.h"

typedef short                               PixelValueType;

ScalarImage4DType::Pointer CreateVolume(VectorImage3DType::Pointer & inputVol)
{
  VectorImage3DType::SizeType inputSize =
    inputVol->GetLargestPossibleRegion().GetSize();
  VectorImage3DType::SpacingType   inputSpacing = inputVol->GetSpacing();
  VectorImage3DType::PointType     inputOrigin = inputVol->GetOrigin();
  VectorImage3DType::DirectionType inputDirection = inputVol->GetDirection();

  ScalarImage4DType::Pointer niftiVolume =
    ScalarImage4DType::New();
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
  volSize[3] = inputVol->GetNumberOfComponentsPerPixel();

  niftiVolume->SetRegions(volSize);
  niftiVolume->SetOrigin(volOrigin);
  niftiVolume->SetSpacing(volSpacing);
  niftiVolume->SetDirection(volDirection);
  niftiVolume->Allocate();
  return niftiVolume;
}


//
// strip the NIfTI suffix from a filename.
static std::string
StripNIfTIName(const std::string &niftiName)
{
  const std::string niigz(".nii.gz");
  if(niftiName.size() > niigz.size() &&
     niftiName.substr(niftiName.size() - 7) == niigz)
    {
    return niftiName.substr(0,niftiName.size() - 7);
    }
  const std::string nii(".nii");
  if(niftiName.size() > nii.size() &&
     niftiName.substr(niftiName.size() - 4) == nii)
    {
    return niftiName.substr(0,niftiName.size() - 4);
    }
  // if it isn't one of the standard suffixes, all you can do is
  // return the whole name.
  return niftiName;
}

int NrrdToFSL(const std::string & inputVolume,
              const std::string & outputVolume,
              const std::string & outputBValues,
              const std::string & outputBVectors,
              bool allowLossyConversion)
{
  if( CheckArg<std::string>("Input Volume", inputVolume, "") == EXIT_FAILURE ||
      CheckArg<std::string>("Output Volume", outputVolume, "") == EXIT_FAILURE)
    {
    return EXIT_FAILURE;
    }

  std::string _outputBValues, _outputBVectors;
  if(CheckArg<std::string>("B Values", outputBValues, "") == EXIT_FAILURE)
    {
    _outputBValues = StripNIfTIName(outputVolume) + ".bval";
    }
  else
    {
    _outputBValues = outputBValues;
    }
  if(CheckArg<std::string>("B Vectors", outputBVectors, "") == EXIT_FAILURE)
    {
    _outputBVectors = StripNIfTIName(outputVolume) + ".bvec";
    }
  else
    {
    _outputBVectors = outputBVectors;
    }

  VectorImage3DType::Pointer inputVol;
  if( ReadVolume<VectorImage3DType>( inputVol, inputVolume, allowLossyConversion ) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  ScalarImage4DType::Pointer                         niftiVolume = CreateVolume(inputVol);
  const VectorImage3DType::SizeType            inputSize( inputVol->GetLargestPossibleRegion().GetSize() );
  const ScalarImage4DType::IndexType::IndexValueType vecLength = inputVol->GetNumberOfComponentsPerPixel();

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
          niftiVolume->SetPixel(volIndex, inputVol->GetPixel(vecIndex)[volIndex[3]]);
          }
        }
      }
    }
  if( WriteVolume<ScalarImage4DType>(niftiVolume, outputVolume) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  DWIMetaDataDictionaryValidator::GradientTableType bVectors;
  if( RecoverBVectors<VectorImage3DType>(inputVol, bVectors) != EXIT_SUCCESS )
    {
    std::cerr << "No gradient vectors found in "
              << inputVolume << std::endl;
    return EXIT_FAILURE;
    }

  if( WriteBVectors(bVectors, _outputBVectors) != EXIT_SUCCESS )
    {
    std::cerr << "Failed to write " << _outputBVectors << std::endl;
    return EXIT_FAILURE;
    }

  std::vector<double> bValues;
  RecoverBValues<VectorImage3DType>(inputVol, bVectors, bValues);

  if( WriteBValues(bValues, _outputBValues) != EXIT_SUCCESS )
    {
    std::cerr << "Failed to write " << _outputBValues << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
