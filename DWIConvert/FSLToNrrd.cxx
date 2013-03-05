#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "DWIConvertUtils.h"

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
          const std::string & inputBValues,
          const std::string & inputBVectors)
{
  if( CheckArg<std::string>("Input Volume", inputVolume, "") == EXIT_FAILURE ||
      CheckArg<std::string>("Output Volume", outputVolume, "") == EXIT_FAILURE ||
      CheckArg<std::string>("B Values", inputBValues, "") == EXIT_FAILURE ||
      CheckArg<std::string>("B Vectors", inputBVectors, "") )
    {
    return EXIT_FAILURE;
    }

  VolumeType::Pointer inputVol;
  if( ReadVolume<VolumeType>(inputVol, inputVolume) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  std::vector<double>               BVals;
  std::vector<std::vector<double> > BVecs;
  unsigned int                      bValCount = 0;
  unsigned int                      bVecCount = 0;
  double                            maxBValue(0.0);
  if( ReadBVals(BVals, bValCount, inputBValues, maxBValue) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  if( ReadBVecs(BVecs, bVecCount, inputBVectors) != EXIT_SUCCESS )
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
  header << "type: short" << std::endl;
  header << "dimension: 4" << std::endl;

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
