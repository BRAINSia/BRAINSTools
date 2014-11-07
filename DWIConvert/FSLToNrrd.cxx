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
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "DWIConvertUtils.h"
#include "itksys/SystemTools.hxx"
#include "vnl/vnl_math.h"
typedef short                               PixelValueType;
typedef itk::Image<PixelValueType, 4>       VolumeType;
typedef itk::VectorImage<PixelValueType, 3> VectorVolumeType;

#include "itkMetaDataObject.h"

typedef short                               PixelValueType;
typedef itk::Image<PixelValueType, 4>       VolumeType;
typedef itk::VectorImage<PixelValueType, 3> VectorVolumeType;

int
FSLToNrrd(const std::string & inputVolume,
          const std::string & outputVolume,
          const std::string & fslNIFTIFile,
          const std::string & inputBValues,
          const std::string & inputBVectors)
{
  if( (CheckArg<std::string>("Input Volume", inputVolume, "") == EXIT_FAILURE &&
       CheckArg<std::string>("Input Volume", fslNIFTIFile, "") == EXIT_FAILURE) ||
      CheckArg<std::string>("Output Volume", outputVolume, "") == EXIT_FAILURE)
    {
    return EXIT_FAILURE;
    }

  VolumeType::Pointer inputVol;

  // string to use as template if no bval or bvec filename is given.
  std::string inputVolumeNameTemplate = inputVolume;
  if(fslNIFTIFile.size() > 0)
    {
    if( ReadVolume<VolumeType>(inputVol, fslNIFTIFile) != EXIT_SUCCESS )
      {
      return EXIT_FAILURE;
      }
    inputVolumeNameTemplate = fslNIFTIFile;
    }
  else if( inputVolume.size() == 0 || ReadVolume<VolumeType>(inputVol, inputVolume) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }

  std::string _inputBValues = inputBValues;
  if( CheckArg<std::string>("B Values", inputBValues, "") == EXIT_FAILURE )
    {
    _inputBValues = itksys::SystemTools::GetFilenameWithoutExtension(inputVolumeNameTemplate) +
      ".bval";
    }
  std::string _inputBVectors = inputBVectors;
  if( CheckArg<std::string>("B Vectors", inputBVectors, "") == EXIT_FAILURE )
    {
    _inputBVectors = itksys::SystemTools::GetFilenameWithoutExtension(inputVolumeNameTemplate) +
      ".bvec";
    }

  std::vector<double>               BVals;
  std::vector<std::vector<double> > BVecs;
  unsigned int                      bValCount = 0;
  unsigned int                      bVecCount = 0;
  double                            maxBValue(0.0);
  if( ReadBVals(BVals, bValCount, _inputBValues, maxBValue) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  if( ReadBVecs(BVecs, bVecCount, _inputBVectors) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  if( bValCount != bVecCount )
    {
    std::cerr << "Mismatch between count of B Vectors ("
              << bVecCount << ") and B Values ("
              << bValCount << ")" << std::endl;
    return EXIT_FAILURE;
    }

  // As suggested by Martin Styner:
  // The implicit normalization is implementing by dividing
  // all the gradients (or B-matrices) by the maximal
  // gradient (or B-matrix) magnitude. The magnitude of a
  // gradient direction vector is determined by the usual
  // L^2 norm, and the magnitude of a B-matrix is via the
  // [Frobenius Norm]. It is after this magnitude rescaling
  // that the nominal b-value (given via
  // "DWMRI_b-value:=b") applies... If the nominal b-value
  // is 1000 (via DWMRI_b-value:=1000), then to represent a
  // DWI with b=500, use a gradient vector (or B-matrix) who's norm is sqrt(1/2) the norm for the b=1000 DWI. "
  //
  // So, since all the b-vecs from FSL have norm 1, all the
  // gradients with maximal b-value (which is your Nrrd
  // DWMRI_b-value) should have norm 1 (i.e. for those you
  // can simply copy over the b-vec info). For all other
  // gradients (i.e. those with b-values below the maximal
  // b-value), you need to scale the coordinates of the
  // b-vector by sqrt(this-b-value/max-b-value).
  std::vector<double>::const_iterator bValIt = BVals.begin(),
    bValsEnd = BVals.end();
  std::vector<std::vector<double> >::iterator bVecIt = BVecs.begin(),
    bVecsEnd = BVecs.end();

  for(; bVecIt != bVecsEnd && bValIt != bValsEnd; ++bVecIt, ++bValIt)
    {
    if((*bValIt) == maxBValue)
      {
      continue;
      }
    double scale = vcl_sqrt((*bValIt) / maxBValue);
    std::vector<double> &cur = *bVecIt;
    for(unsigned int i = 0; i < 3; ++i)
      {
      cur[i] *= scale;
      }
    }

  VolumeType::SizeType inputSize =
    inputVol->GetLargestPossibleRegion().GetSize();

  const unsigned int volumeCount = inputSize[3];
  if( volumeCount != bValCount )
    {
    std::cerr << "Mismatch between BVector count ("
              << bVecCount << ") and image volume count ("
              << volumeCount << ")" << std::endl;
    return EXIT_SUCCESS;
    }

  // convert from image series to vector voxels
  VolumeType::SpacingType   inputSpacing = inputVol->GetSpacing();
  VolumeType::PointType     inputOrigin = inputVol->GetOrigin();
  VolumeType::DirectionType inputDirection = inputVol->GetDirection();

  std::ofstream header;
  // std::string headerFileName = outputDir + "/" + outputFileName;

  header.open(outputVolume.c_str(), std::ios::out | std::ios::binary);
  header.precision(17);
  header << "NRRD0005" << std::endl;
  header << "# This file was created by DWIConvert version 1.0" << std::endl
         << "# https://github.com/BRAINSia/BRAINSTools" << std::endl
         << "# part of the BRAINSTools package." << std::endl
         << "type: short" << std::endl
         << "dimension: 4" << std::endl;

  // need to check
  header << "space: left-posterior-superior" << std::endl;
  // in nrrd, size array is the number of pixels in 1st, 2nd, 3rd, ... dimensions
  header << "sizes: " << inputSize[0] << " "
         << inputSize[1] << " "
         << inputSize[2] << " "
         << inputSize[3] << std::endl;
  header << "thicknesses:  NaN  NaN " << inputSpacing[2] << " NaN" << std::endl;
  double spaceDirections[3][3];
  for( unsigned i = 0; i < 3; ++i )
    {
    for( unsigned j = 0; j < 3; ++j )
      {
      spaceDirections[i][j] = inputDirection[i][j];
      if( i == j )
        {
        spaceDirections[i][j] *= inputSpacing[i];
        }
      }
    }
  // need to check
  header << "space directions: "
         << "(" << spaceDirections[0][0] << ","
         << spaceDirections[1][0] << ","
         << spaceDirections[2][0] << ") "
         << "(" << spaceDirections[0][1] << ","
         << spaceDirections[1][1] << ","
         << spaceDirections[2][1] << ") "
         << "(" << spaceDirections[0][2] << ","
         << spaceDirections[1][2] << ","
         << spaceDirections[2][2] << ") none"
         << std::endl;
  header << "centerings: cell cell cell ???" << std::endl;
  header << "kinds: space space space list" << std::endl;

  header << "endian: little" << std::endl;
  header << "encoding: raw" << std::endl;
  header << "space units: \"mm\" \"mm\" \"mm\"" << std::endl;
  header << "space origin: "
         << "(" << inputOrigin[0]
         << "," << inputOrigin[1]
         << "," << inputOrigin[2] << ") " << std::endl;
  header << "measurement frame: "
         << "(" << 1 << "," << 0 << "," << 0 << ") "
         << "(" << 0 << "," << 1 << "," << 0 << ") "
         << "(" << 0 << "," << 0 << "," << 1 << ")"
         << std::endl;

  header << "modality:=DWMRI" << std::endl;
  // this is the norminal BValue, i.e. the largest one.
  header << "DWMRI_b-value:=" << maxBValue << std::endl;
  for( unsigned int i = 0; i < bVecCount; ++i )
    {
    header << "DWMRI_gradient_" << std::setw(4) << std::setfill('0')
           << i << ":="
           << BVecs[i][0] << "   "
           << BVecs[i][1] << "   "
           << BVecs[i][2]
           << std::endl;
    }

  // write data in the same file is .nrrd was chosen
  header << std::endl;;
  unsigned long nVoxels = inputVol->GetLargestPossibleRegion().GetNumberOfPixels();
  header.write( reinterpret_cast<char *>(inputVol->GetBufferPointer() ),
                nVoxels * sizeof(short) );
  header.close();
  return EXIT_SUCCESS;
}
