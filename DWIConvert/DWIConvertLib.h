//
// Created by Hui Xie on 12/26/16.
//

#ifndef BRAINSTOOLS_DWICONVERTLIB_H
#define BRAINSTOOLS_DWICONVERTLIB_H

#include "DWIConverter.h"
#include "DWIConverterFactory.h"


struct DWIConvertParameters {
  std::string inputVolume;
  std::string inputDicomDirectory;
  std::string inputBValues;
  std::string inputBVectors;
  std::string gradientVectorFile; //deprecated
  double smallGradientThreshold;

  std::string conversionMode; // only one of ["DicomToNrrd", "DicomToFSL", "NrrdToFSL", "FSLToNrrd"]
  bool fMRIOutput;
  bool transpose;
  bool allowLossyConversion;
  bool useIdentityMeasurementFrame;
  bool useBMatrixGradientDirections;

  std::string outputVolume;
  std::string outputDirectory;
  std::string outputBValues;
  std::string outputBVectors;
};

//rewrite DWIConvert1, there is a test failure needing to debug
//Please use DWIConvert2 interface
int DWIConvert1(const DWIConvertParameters& params);

DWIConverter * CreateDicomConverter(
        const std::string inputDicomDirectory,
        const bool useBMatrixGradientDirections,
        const bool transpose,
        const double smallGradientThreshold,
        const bool allowLossyConversion);

//according Hans's code, encapsulate the DWIConvert2
int DWIConvert2(const DWIConvertParameters& params);


#endif //BRAINSTOOLS_DWICONVERTLIB_H