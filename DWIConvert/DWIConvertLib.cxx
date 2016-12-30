//
// Created by Hui Xie on 12/26/16.
//

#include "DWIConvertLib.h"

#include "dcmtk/oflog/helpers/loglog.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmjpeg/djdecode.h"
#include "dcmtk/dcmjpls/djdecode.h"
#include "dcmtk/dcmdata/dcrledrg.h"

DWIConverter * DWIConvert::CreateDicomConverter(
        const std::string inputDicomDirectory,
        const bool useBMatrixGradientDirections,
        const bool transpose,
        const double smallGradientThreshold,
        const bool allowLossyConversion)
{
// check for required parameters
  if( inputDicomDirectory == "" )
  {
    std::cerr << "Missing DICOM input directory path" << std::endl;
    return ITK_NULLPTR;
  }

// use the factor to instantiate a converter object based on the vender.
  DWIConverterFactory converterFactory(inputDicomDirectory,
                                       useBMatrixGradientDirections,
                                       transpose,
                                       smallGradientThreshold);
  DWIConverter * converter;
  try
  {
    converter = converterFactory.New();
  }
  catch( itk::ExceptionObject &excp)
  {
    std::cerr << "Exception creating converter " << excp << std::endl;
    return ITK_NULLPTR;
  }

// read Dicom directory
  try
  {
    converter->SetAllowLossyConversion(allowLossyConversion);
    converter->LoadFromDisk();
  }
  catch( itk::ExceptionObject &excp)
  {
    std::cerr << "Exception creating converter " << excp << std::endl;
    delete converter;
    return ITK_NULLPTR;
  }
  // extract the DWI data
  try
  {
    converter->ExtractDWIData();
  }
  catch( itk::ExceptionObject &excp)
  {
    std::cerr << "Exception extracting gradient vectors " << excp << std::endl;
    delete converter;
    return ITK_NULLPTR;
  }
// this is a punt, it will still write out the volume image
// even if we don't know how to extract gradients.
  if(converterFactory.GetVendor() == "GENERIC")
  {
    std::cerr << "Can't extract DWI data from files created by vendor "
              << converterFactory.GetVendor() << std::endl;
    delete converter;
    exit(EXIT_FAILURE);
  }
  return converter;
}

DWIConvert::DWIConvert()
{
    m_inputVolume = "";
    m_inputDicomDirectory = "";
    m_inputBValues = "";  //default: ""
    m_inputBVectors = ""; //default: ""
    m_gradientVectorFile = ""; //deprecated
    m_smallGradientThreshold = 0.2; //default = 0.2

    //only one of ["DicomToNrrd", "DicomToFSL", "NrrdToFSL", "FSLToNrrd",
    //                                          "NrrdToNrrd", "FSLToFSL"]
    m_conversionMode = "";

    m_fMRIOutput = false; //default: false
    m_transpose = false; //default:false
    m_allowLossyConversion = false; //defualt: false
    m_useIdentityMeasurementFrame = false; //default: false
    m_useBMatrixGradientDirections = false; //default: false

    m_outputVolume = "";
    m_outputDirectory = ".";  //default: "."
    m_outputBValues = ""; //default: ""
    m_outputBVectors = "";//default: ""

}

int DWIConvert::DWIConvert1() {
  const std::string version = "4.8.0";
  std::vector<std::string> convertModeVector;
  convertModeVector.reserve(4);
  convertModeVector.push_back("DicomToNrrd");
  convertModeVector.push_back("DicomToFSL");
  convertModeVector.push_back("NrrdToFSL");
  convertModeVector.push_back("FSLToNrrd");
  bool useIdentityMeasurementFrame = m_useIdentityMeasurementFrame;

  //check parameters
  if (m_fMRIOutput) {
    std::cerr << "Deprecated feature no longer supported: --fMRIOutput" << std::endl;
    return EXIT_FAILURE;
  }
  if (m_gradientVectorFile != "") {
    std::cerr << "Deprecated feature no longer supported: --gradientVectorFile" << std::endl;
    return EXIT_FAILURE;
  }

  if (m_outputVolume == "") {
    std::cerr << "Missing output volume name" << std::endl;
    return EXIT_FAILURE;
  }

  //check conversion mode
  bool correctConvertMode = false;
  for (int i = 0; i < 4; ++i) {
    if (0 == convertModeVector[i].compare(getConversionMode())) {
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
  filesList.push_back(m_inputVolume);

  try {
    if (0 == getConversionMode().compare("FSLToNrrd")) {
      converter = new FSLDWIConverter(filesList, m_inputBValues, m_inputBVectors, m_transpose);
    } else if (0 == getConversionMode().compare("NrrdToFSL")) {
      converter = new NRRDDWIConverter(filesList, m_transpose);
      useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
    } else  // "DicomToNrrd" or "DicomToFSL"
    {
      if (m_inputDicomDirectory == "") {
        std::cerr << "Missing Dicom input directory path" << std::endl;
        return EXIT_FAILURE;
      }
      if (getConversionMode() == "DicomToFSL") {
        useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
      }

      DWIConverterFactory converterFactory(m_inputDicomDirectory,
                                           m_useBMatrixGradientDirections,
                                           m_transpose,
                                           m_smallGradientThreshold);
      converter = converterFactory.New();

      // this is a punt, it will still write out the volume image
      // even if we don't know how to extract gradients.
      if(converterFactory.GetVendor() == "GENERIC")
      {
        std::cerr << "Can't extract DWI data from files created by vendor Generic"<< std::endl;
        delete converter;
        return EXIT_FAILURE;
      }
    }
  }
  catch (std::exception &e) {
    std::cerr << "Exception in creating converter: " << e.what() << std::endl;
    if (NULL != converter) delete converter;
    return EXIT_FAILURE;
  }

  try {
      //extract DWI data
      converter->SetAllowLossyConversion(m_allowLossyConversion);
  }
  catch (std::exception &e) {
      std::cerr << "Exception in SetAllowLossyConversion: " << e.what() << std::endl;
      if (NULL != converter) delete converter;
      return EXIT_FAILURE;
  }

  try {
    converter->LoadFromDisk();
  }
  catch (std::exception &e) {
      std::cerr << "Exception in LoadFromDisk: " << e.what() << std::endl;
      if (NULL != converter) delete converter;
      return EXIT_FAILURE;
  }

  try {
      converter->ExtractDWIData();
  }
  catch (std::exception &e) {
      std::cerr << "Exception in extracting DWI data: " << e.what() << std::endl;
      if (NULL != converter) delete converter;
      return EXIT_FAILURE;
  }

  //Write output
  if (useIdentityMeasurementFrame) {
    converter->ConvertBVectorsToIdentityMeasurementFrame();
  }
  std::string outputVolumeHeaderName(m_outputVolume);
  // concatenate with outputDirectory
  if (m_outputVolume.find("/") == std::string::npos
      && m_outputVolume.find("\\") == std::string::npos
      && outputVolumeHeaderName.size() != 0) {
    outputVolumeHeaderName = m_outputDirectory;
    outputVolumeHeaderName += "/";
    outputVolumeHeaderName += m_outputVolume;
  }

  if (getConversionMode() == "DicomToFSL" || getConversionMode() == "NrrdToFSL") {

    Volume4DType::Pointer img4D = converter->OrientForFSLConventions();
    // write the image */
    converter->WriteFSLFormattedFileSet(outputVolumeHeaderName, m_outputBValues, m_outputBVectors, img4D);
  } else {
    //NRRD requires scaledDiffusionVectors, so find max BValue
    converter->ConvertToSingleBValueScaledDiffusionVectors();
    const std::string commentSection = converter->MakeFileComment(version, m_useBMatrixGradientDirections,
                                                                  useIdentityMeasurementFrame,
                                                                  m_smallGradientThreshold,getConversionMode());

    converter->ManualWriteNRRDFile(outputVolumeHeaderName, commentSection);
    std::cout << "Wrote file: " << outputVolumeHeaderName << std::endl;
  }
  if (NULL != converter) delete converter;
  return EXIT_SUCCESS;
}


int DWIConvert::DWIConvert2()
{
  const std::string version = "4.8.0";
  dcmtk::log4cplus::helpers::LogLog::getLogLog()->setQuietMode(true);

  // register DCMTK codecs, otherwise they will not be available when
  // `itkDCMTKSeriesFileNames` is used to build a list of filenames,
  // so reading series with JPEG transfer syntax will fail.
  DJDecoderRegistration::registerCodecs();
  DcmRLEDecoderRegistration::registerCodecs();


  bool useIdentityMeasurementFrame = m_useIdentityMeasurementFrame;
  if(m_fMRIOutput)
  {
    std::cerr << "Deprecated feature no longer supported: --fMRIOutput" << std::endl;
    return EXIT_FAILURE;
  }
  if( m_gradientVectorFile != "" ) {
    std::cerr << "Deprecated feature no longer supported: --gradientVectorFile" << std::endl;
    return EXIT_FAILURE;
  }

  if( m_outputVolume == "" )
  {
    std::cerr << "Missing output volume name" << std::endl;
    return EXIT_FAILURE;
  }

  DWIConverter *converter;
  // build a NRRD file out of FSL output, which is two text files
  // for gradients and b values plus a NIfTI file for the gradient volumes.
  if( getConversionMode() == "FSLToNrrd" )
  {
    DWIConverter::FileNamesContainer filesList;
    filesList.clear();
    filesList.push_back(m_inputVolume);

    std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
    FSLDWIConverter * FSLconverter = new FSLDWIConverter(filesList, m_inputBValues, m_inputBVectors,m_transpose);
    try
    {
      FSLconverter->SetAllowLossyConversion(m_allowLossyConversion);
      FSLconverter->LoadFromDisk();
      FSLconverter->ExtractDWIData();
    }
    catch( itk::ExceptionObject &excp)
    {
      itkGenericExceptionMacro( << "Exception creating converter " << excp << std::endl);
      delete FSLconverter;
    }
    converter = FSLconverter;
  }
    // make FSL file set from a NRRD file.
  else if( getConversionMode() == "NrrdToFSL" )
  {
    DWIConverter::FileNamesContainer filesList;
    filesList.clear();
    filesList.push_back(m_inputVolume);

    std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
    NRRDDWIConverter * NRRDconverter = new NRRDDWIConverter(filesList, m_transpose);
    useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
    try
    {
      NRRDconverter->SetAllowLossyConversion(m_allowLossyConversion);
      NRRDconverter->LoadFromDisk();
      NRRDconverter->ExtractDWIData();
    }
    catch( itk::ExceptionObject &excp)
    {
      itkGenericExceptionMacro( << "Exception creating converter " << excp << std::endl);
      delete NRRDconverter;
    }
    converter = NRRDconverter;
  }
  else if( getConversionMode() == "DicomToNrrd" || getConversionMode() == "DicomToFSL")
  {
    converter = CreateDicomConverter(m_inputDicomDirectory,m_useBMatrixGradientDirections, m_transpose,
                                     m_smallGradientThreshold,m_allowLossyConversion);
    if (getConversionMode() == "DicomToFSL")
    {
      useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
    }
  }
  else
  {
    std::cerr << "Invalid conversion mode" << std::endl;
    exit(-1);
  }

//#if 0 //This should use the bvec and bval file formats
//  // NEED TO ADD --forceGradientOverwrite, and then read bvec and bval files
//  // A test needs to be written for this case
//  // ^^^^^^^^^^^^^^^^^^^^^^^ Done Reading Above this line
//  //Overwrite gradient directions
//  if( gradientVectorFile != "" )
//  {
//    converter->readOverwriteGradientVectorFile(gradientVectorFile);
//  }
//#endif
  //^^^^^^^^^^^^^^^^^^^^^^^^^Done modifying above this line vvvvvvvvvvvvvvvvvvvvv Write outputs
  if (useIdentityMeasurementFrame)
  {
    converter->ConvertBVectorsToIdentityMeasurementFrame();
  }



  std::string outputVolumeHeaderName(m_outputVolume);
  { // concatenate with outputDirectory
    if( m_outputVolume.find("/") == std::string::npos &&
            m_outputVolume.find("\\") == std::string::npos )
    {
      if( outputVolumeHeaderName.size() != 0 )
      {
        outputVolumeHeaderName = m_outputDirectory;
        outputVolumeHeaderName += "/";
        outputVolumeHeaderName += m_outputVolume;
      }
    }
  }

  if( getConversionMode() == "DicomToFSL"  || getConversionMode() == "NrrdToFSL" )
  {

    Volume4DType::Pointer img4D = converter->OrientForFSLConventions();
    // write the image */
    converter->WriteFSLFormattedFileSet(outputVolumeHeaderName, m_outputBValues, m_outputBVectors, img4D);
  }
  else
  {
    //NRRD requires scaledDiffusionVectors, so find max BValue
    converter->ConvertToSingleBValueScaledDiffusionVectors();
    //////////////////////////////////////////////
    // write header file
    // This part follows a DWI NRRD file in NRRD format 5.
    // There should be a better way using itkNRRDImageIO.
    const std::string commentSection = converter->MakeFileComment(version,m_useBMatrixGradientDirections,
                                                                  useIdentityMeasurementFrame, m_smallGradientThreshold, getConversionMode());

    converter->ManualWriteNRRDFile(outputVolumeHeaderName, commentSection);
    std::cout << "Wrote file: " << outputVolumeHeaderName << std::endl;
  }
  delete converter;
  return EXIT_SUCCESS;


}

//return: only one of ["DicomToNrrd", "DicomToFSL", "NrrdToFSL", "FSLToNrrd",
//                                                  "NrrdToNrrd", "FSLToFSL"]
// if return "", invalidate input parameters.
std::string DWIConvert::setConversionMode()
{
  std::string conversionMode= "";
  std::string outputExt = findFilenameExt(m_outputVolume);
  if ("" == m_inputDicomDirectory){
     std::string inputExt = findFilenameExt(m_inputVolume);
     if (".nii" == inputExt && ".nrrd" == outputExt) conversionMode = "FSLToNrrd";
     if (".nrrd" == inputExt && ".nii" == outputExt) conversionMode = "NrrdToFSL";
     if (".nrrd" == inputExt && ".nrrd" == outputExt) conversionMode = "NrrdToNrrd";
     if (".nii" == inputExt && ".nii" == outputExt) conversionMode = "FSLToFSL";
  }
  else if ("" == m_inputVolume)
  {
      if (".nrrd" == outputExt) conversionMode = "DicomToNrrd";
      if (".nii" == outputExt) conversionMode = "DicomToFSL";
  }
  else{
      std::cerr <<"Error: One of inputVoume and inputDicomDirectory must be empty string"<<std::endl;
  }
  if ("" == conversionMode){
      std::cerr <<"Error: invalidate input or output file tyep"<<std::endl;
  }
  m_conversionMode =  conversionMode;
  return conversionMode;
}

void DWIConvert::setConversionMode(const std::string conversionMode){
    m_conversionMode =  conversionMode;
}

std::string DWIConvert::getConversionMode()
{
  if ("" == m_conversionMode) setConversionMode();
  return m_conversionMode;
}

//{ ".nii", ".nii.gz", ".nhdr", ".nrrd"}
std::string DWIConvert::findFilenameExt(const std::string filename){
    std::string::size_type pos = filename.rfind(".");
    std::string subStr = filename.substr(pos);
    if (".gz" == subStr){
        subStr = filename.substr(pos-4,4);
    }
    return subStr;
}

const std::string &DWIConvert::getinputVolume() const {
  return m_inputVolume;
}

void DWIConvert::setM_inputVolume(const std::string &m_inputVolume) {
  DWIConvert::m_inputVolume = m_inputVolume;
}

const std::string &DWIConvert::getinputDicomDirectory() const {
  return m_inputDicomDirectory;
}

void DWIConvert::setM_inputDicomDirectory(const std::string &m_inputDicomDirectory) {
  DWIConvert::m_inputDicomDirectory = m_inputDicomDirectory;
}

const std::string &DWIConvert::getinputBValues() const {
  return m_inputBValues;
}

void DWIConvert::setM_inputBValues(const std::string &m_inputBValues) {
  DWIConvert::m_inputBValues = m_inputBValues;
}

const std::string &DWIConvert::getInputBVectors() const {
  return m_inputBVectors;
}

void DWIConvert::setM_inputBVectors(const std::string &m_inputBVectors) {
  DWIConvert::m_inputBVectors = m_inputBVectors;
}

const std::string &DWIConvert::getgradientVectorFile() const {
  return m_gradientVectorFile;
}

void DWIConvert::setM_gradientVectorFile(const std::string &m_gradientVectorFile) {
  DWIConvert::m_gradientVectorFile = m_gradientVectorFile;
}

double DWIConvert::getSmallGradientThreshold() const {
  return m_smallGradientThreshold;
}

void DWIConvert::setSmallGradientThreshold(double m_smallGradientThreshold) {
  DWIConvert::m_smallGradientThreshold = m_smallGradientThreshold;
}

bool DWIConvert::isM_fMRIOutput() const {
  return m_fMRIOutput;
}

void DWIConvert::setM_fMRIOutput(bool m_fMRIOutput) {
  DWIConvert::m_fMRIOutput = m_fMRIOutput;
}

bool DWIConvert::isM_transpose() const {
  return m_transpose;
}

void DWIConvert::setM_transpose(bool m_transpose) {
  DWIConvert::m_transpose = m_transpose;
}

bool DWIConvert::isM_allowLossyConversion() const {
  return m_allowLossyConversion;
}

void DWIConvert::setM_allowLossyConversion(bool m_allowLossyConversion) {
  DWIConvert::m_allowLossyConversion = m_allowLossyConversion;
}

bool DWIConvert::isM_useIdentityMeasurementFrame() const {
  return m_useIdentityMeasurementFrame;
}

void DWIConvert::setM_useIdentityMeasurementFrame(bool m_useIdentityMeasurementFrame) {
  DWIConvert::m_useIdentityMeasurementFrame = m_useIdentityMeasurementFrame;
}

bool DWIConvert::isM_useBMatrixGradientDirections() const {
  return m_useBMatrixGradientDirections;
}

void DWIConvert::setM_useBMatrixGradientDirections(bool m_useBMatrixGradientDirections) {
  DWIConvert::m_useBMatrixGradientDirections = m_useBMatrixGradientDirections;
}

const std::string &DWIConvert::getoutputVolume() const {
  return m_outputVolume;
}

void DWIConvert::setM_outputVolume(const std::string &m_outputVolume) {
  DWIConvert::m_outputVolume = m_outputVolume;
}

const std::string &DWIConvert::getoutputDirectory() const {
  return m_outputDirectory;
}

void DWIConvert::setM_outputDirectory(const std::string &m_outputDirectory) {
  DWIConvert::m_outputDirectory = m_outputDirectory;
}

const std::string &DWIConvert::getoutputBValues() const {
  return m_outputBValues;
}

void DWIConvert::setM_outputBValues(const std::string &m_outputBValues) {
  DWIConvert::m_outputBValues = m_outputBValues;
}

const std::string &DWIConvert::getoutputBVectors() const {
  return m_outputBVectors;
}

void DWIConvert::setM_outputBVectors(const std::string &m_outputBVectors) {
  DWIConvert::m_outputBVectors = m_outputBVectors;
}