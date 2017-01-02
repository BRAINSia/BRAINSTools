//
// Created by Hui Xie on 12/26/16.
//

#include "DWIConvertLib.h"

#include "dcmtk/oflog/helpers/loglog.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmjpeg/djdecode.h"
#include "dcmtk/dcmjpls/djdecode.h"
#include "dcmtk/dcmdata/dcrledrg.h"



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

    m_converter = NULL;

}

DWIConvert::DWIConvert(const std::string& outputVolume, const std::string& inputVolume, const std::string& inputDicomDirectory)
{
   DWIConvert();
   setOutputVolume(outputVolume);
   setInputVolume(inputVolume);
   setInputDicomDirectory(inputDicomDirectory);
   setConversionMode();
}

DWIConvert::~DWIConvert(){
  delete m_converter;
}

/*int DWIConvert::DWIConvert1() {
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
    // write the image *//*
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
 */


int DWIConvert::read()
{
  if ("" == getConversionMode()) setConversionMode();

  dcmtk::log4cplus::helpers::LogLog::getLogLog()->setQuietMode(true);
  // register DCMTK codecs, otherwise they will not be available when
  // `itkDCMTKSeriesFileNames` is used to build a list of filenames,
  // so reading series with JPEG transfer syntax will fail.
  DJDecoderRegistration::registerCodecs();
  DcmRLEDecoderRegistration::registerCodecs();

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

  // build a NRRD file out of FSL output, which is two text files
  // for gradients and b values plus a NIfTI file for the gradient volumes.
  if( getConversionMode() == "FSLToNrrd" || "FSLToFSL" == getConversionMode())
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
    m_converter = FSLconverter;
  }
    // make FSL file set from a NRRD file.
  else if( getConversionMode() == "NrrdToFSL" || "NrrdToNrrd" ==  getConversionMode() )
  {
    DWIConverter::FileNamesContainer filesList;
    filesList.clear();
    filesList.push_back(m_inputVolume);

    std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
    NRRDDWIConverter * NRRDconverter = new NRRDDWIConverter(filesList, m_transpose);
    m_useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
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
    m_converter = NRRDconverter;
  }
  else if( getConversionMode() == "DicomToNrrd" || getConversionMode() == "DicomToFSL")
  {
    m_converter = CreateDicomConverter(m_inputDicomDirectory,m_useBMatrixGradientDirections, m_transpose,
                                     m_smallGradientThreshold,m_allowLossyConversion);
    if (getConversionMode() == "DicomToFSL")
    {
      m_useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
    }
  }
  else
  {
    std::cerr << "Invalid conversion mode" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}

int DWIConvert::write()
{
  const std::string version = "4.8.0";
  if (m_useIdentityMeasurementFrame)
  {
    m_converter->ConvertBVectorsToIdentityMeasurementFrame();
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

  if( getConversionMode() == "DicomToFSL"  || getConversionMode() == "NrrdToFSL"  || "FSLToFSL" == getConversionMode())
  {

    Volume4DType::Pointer img4D = m_converter->OrientForFSLConventions();
    // write the image */
    m_converter->WriteFSLFormattedFileSet(outputVolumeHeaderName, m_outputBValues, m_outputBVectors, img4D);
  }
  else if (getConversionMode() == "DicomToNrrd"  || getConversionMode() == "NrrdToNrrd"  || "FSLToNrrd" == getConversionMode())
  {
      //NRRD requires scaledDiffusionVectors, so find max BValue
      m_converter->ConvertToSingleBValueScaledDiffusionVectors();
      //////////////////////////////////////////////
      // write header file
      // This part follows a DWI NRRD file in NRRD format 5.
      // There should be a better way using itkNRRDImageIO.
      const std::string commentSection = m_converter->MakeFileComment(version, m_useBMatrixGradientDirections,
                                                                      m_useIdentityMeasurementFrame,
                                                                      m_smallGradientThreshold, getConversionMode());
      m_converter->ManualWriteNRRDFile(outputVolumeHeaderName, commentSection);
  }
  else
  {
     std::cerr<<"Invalidate conversion mode" << std::endl;
     return EXIT_FAILURE;
  }
  std::cout << "Wrote file: " << outputVolumeHeaderName << std::endl;
  return EXIT_SUCCESS;
}

int DWIConvert::readWrite()
{
  int result = read();
  if (EXIT_SUCCESS == result ){
    return write();
  }
  else return result;
}


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

//one of ["DicomToNrrd", "DicomToFSL", "NrrdToFSL", "FSLToNrrd", "NrrdToNrrd", "FSLToFSL"]
//if return "", invalidate input parameters.
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
    //["DicomToNrrd", "DicomToFSL", "NrrdToFSL", "FSLToNrrd","NrrdToNrrd", "FSLToFSL"]
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

const std::string &DWIConvert::getInputVolume() const {
  return m_inputVolume;
}

void DWIConvert::setInputVolume(const std::string &inputVolume) {
  m_inputVolume = inputVolume;
}

const std::string &DWIConvert::getInputDicomDirectory() const {
  return m_inputDicomDirectory;
}

void DWIConvert::setInputDicomDirectory(const std::string &inputDicomDirectory) {
  m_inputDicomDirectory = inputDicomDirectory;
}

const std::string &DWIConvert::getInputBValues() const {
  return m_inputBValues;
}

void DWIConvert::setInputBValues(const std::string &inputBValues) {
  m_inputBValues = inputBValues;
}

const std::string &DWIConvert::getInputBVectors() const {
  return m_inputBVectors;
}

void DWIConvert::setInputBVectors(const std::string &inputBVectors) {
  m_inputBVectors = inputBVectors;
}

const std::string &DWIConvert::getGradientVectorFile() const {
  return m_gradientVectorFile;
}

void DWIConvert::setGradientVectorFile(const std::string &gradientVectorFile) {
  m_gradientVectorFile = gradientVectorFile;
}

double DWIConvert::getSmallGradientThreshold() const {
  return m_smallGradientThreshold;
}

void DWIConvert::setSmallGradientThreshold(double smallGradientThreshold) {
  m_smallGradientThreshold = smallGradientThreshold;
}

bool DWIConvert::isfMRIOutput() const {
  return m_fMRIOutput;
}

void DWIConvert::setfMRIOutput(bool fMRIOutput) {
  m_fMRIOutput = fMRIOutput;
}

bool DWIConvert::isTranspose() const {
  return m_transpose;
}

void DWIConvert::setTranspose(bool transpose) {
  m_transpose = transpose;
}

bool DWIConvert::isAllowLossyConversion() const {
  return m_allowLossyConversion;
}

void DWIConvert::setAllowLossyConversion(bool allowLossyConversion) {
  m_allowLossyConversion = allowLossyConversion;
}

bool DWIConvert::isUseIdentityMeasurementFrame() const {
  return m_useIdentityMeasurementFrame;
}

void DWIConvert::setUseIdentityMeasurementFrame(bool useIdentityMeasurementFrame) {
  m_useIdentityMeasurementFrame = useIdentityMeasurementFrame;
}

bool DWIConvert::isUseBMatrixGradientDirections() const {
  return m_useBMatrixGradientDirections;
}

void DWIConvert::setUseBMatrixGradientDirections(bool useBMatrixGradientDirections) {
  m_useBMatrixGradientDirections = useBMatrixGradientDirections;
}

const std::string &DWIConvert::getOutputVolume() const {
  return m_outputVolume;
}

void DWIConvert::setOutputVolume(const std::string &outputVolume) {
  m_outputVolume = outputVolume;
}

const std::string &DWIConvert::getOutputDirectory() const {
  return m_outputDirectory;
}

void DWIConvert::setOutputDirectory(const std::string &outputDirectory) {
  m_outputDirectory = outputDirectory;
}

const std::string &DWIConvert::getOutputBValues() const {
  return m_outputBValues;
}

void DWIConvert::setOutputBValues(const std::string &outputBValues) {
  m_outputBValues = outputBValues;
}

const std::string &DWIConvert::getOutputBVectors() const {
  return m_outputBVectors;
}

void DWIConvert::setOutputBVectors(const std::string &outputBVectors) {
  m_outputBVectors = outputBVectors;
}