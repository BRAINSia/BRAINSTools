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
  std::string gradientVectorFile;
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

int DWIConvert(const DWIConvertParameters& params);

#endif //BRAINSTOOLS_DWICONVERTLIB_H