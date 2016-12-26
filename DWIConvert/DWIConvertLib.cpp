//
// Created by Hui Xie on 12/26/16.
//

#include "DWIConvertLib.h"

int DWIConvert(const DWIConvertParameters &params) {
  const std::string version = "Jan 2017";
  std::vector<std::string> convertModeVector;
  convertModeVector.reserve(4);
  convertModeVector.push_back("DicomToNrrd");
  convertModeVector.push_back("DicomToFSL");
  convertModeVector.push_back("NrrdToFSL");
  convertModeVector.push_back("FSLToNrrd");
  bool useIdentityMeasurementFrame = params.useIdentityMeasurementFrame;

  //check parameters
  if (params.fMRIOutput) {
    std::cerr << "Deprecated feature no longer supported: --fMRIOutput" << std::endl;
    return EXIT_FAILURE;
  }
  if (params.gradientVectorFile != "") {
    std::cerr << "Deprecated feature no longer supported: --gradientVectorFile" << std::endl;
    return EXIT_FAILURE;
  }

  if (params.outputVolume == "") {
    std::cerr << "Missing output volume name" << std::endl;
    return EXIT_FAILURE;
  }

  //check conversion mode
  bool correctConvertMode = false;
  for (int i = 0; i < 4; ++i) {
    if (0 == convertModeVector[i].compare(params.conversionMode)) {
      correctConvertMode = true;
      break;
    }
  }
  if (!correctConvertMode) {
    std::cerr << "DWI Conversion mode is invalid. " << std::endl;
    return EXIT_FAILURE;
  }

  //Create DWIConverter according different conversion mode
  DWIConverter *converter = NULL;
  DWIConverter::FileNamesContainer filesList;
  filesList.clear();
  filesList.push_back(params.inputVolume);

  try {
    if (0 == params.conversionMode.compare("FSLToNrrd")) {
      converter = new FSLDWIConverter(filesList, params.inputBValues, params.inputBVectors, params.transpose);
    } else if (0 == params.conversionMode.compare("NrrdToFSL")) {
      converter = new NRRDDWIConverter(filesList, params.transpose);
      useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
    } else  // "DicomToNrrd" or "DicomToFSL"
    {
      if (params.inputDicomDirectory == "") {
        std::cerr << "Missing Dicom input directory path" << std::endl;
        return EXIT_FAILURE;
      }
      if (params.conversionMode == "DicomToFSL") {
        useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
      }

      DWIConverterFactory converterFactory(params.inputDicomDirectory,
                                           params.useBMatrixGradientDirections,
                                           params.transpose,
                                           params.smallGradientThreshold);
      converter = converterFactory.New();
    }

    //extract DWI data
    converter->SetAllowLossyConversion(params.allowLossyConversion);
    converter->LoadFromDisk();
    converter->ExtractDWIData();
  }
  catch (std::exception &e) {
    std::cerr << "Exception extracting DWI data" << e.what() << std::endl;
    delete converter;
    return EXIT_FAILURE;
  }

  //Write output
  if (useIdentityMeasurementFrame) {
    converter->ConvertBVectorsToIdentityMeasurementFrame();
  }
  std::string outputVolumeHeaderName(params.outputVolume);
  // concatenate with outputDirectory
  if (params.outputVolume.find("/") == std::string::npos
      && params.outputVolume.find("\\") == std::string::npos
      && outputVolumeHeaderName.size() != 0) {
    outputVolumeHeaderName = params.outputDirectory;
    outputVolumeHeaderName += "/";
    outputVolumeHeaderName += params.outputVolume;
  }

  if (params.conversionMode == "DicomToFSL" || params.conversionMode == "NrrdToFSL") {

    Volume4DType::Pointer img4D = converter->OrientForFSLConventions();
    // write the image */
    converter->WriteFSLFormattedFileSet(outputVolumeHeaderName, params.outputBValues, params.outputBVectors, img4D);
  } else {
    //NRRD requires scaledDiffusionVectors, so find max BValue
    converter->ConvertToSingleBValueScaledDiffusionVectors();
    const std::string commentSection = converter->MakeFileComment(version, params.useBMatrixGradientDirections,
                                                                  useIdentityMeasurementFrame,
                                                                  params.smallGradientThreshold, params.conversionMode);

    converter->ManualWriteNRRDFile(outputVolumeHeaderName, commentSection);
    std::cout << "Wrote file: " << outputVolumeHeaderName << std::endl;
  }
  return EXIT_SUCCESS;
}