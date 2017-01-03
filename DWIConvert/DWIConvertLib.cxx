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
    m_inputFileType = "";
    m_outputFileType = "";

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
   setInputFileType(inputVolume,inputDicomDirectory);
   setOutputFileType(outputVolume);
}

DWIConvert::DWIConvert(const std::string& inputVolume, const std::string& inputDicomDirectory)
{
  DWIConvert();
  setInputFileType(inputVolume,inputDicomDirectory);
}

DWIConvert::~DWIConvert(){
  delete m_converter;
}

int DWIConvert::read()
{
  if ("" == getInputFileType()){
    std::cerr << "illegal input file type, exit" << std::endl;
    return EXIT_FAILURE;
  }

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
  if( "FSL" == getInputFileType())
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
  else if( "Nrrd" ==  getInputFileType() )
  {
    DWIConverter::FileNamesContainer filesList;
    filesList.clear();
    filesList.push_back(m_inputVolume);

    std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
    NRRDDWIConverter * NRRDconverter = new NRRDDWIConverter(filesList, m_transpose);
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
  else if( "Dicom" ==  getInputFileType())
  {
    m_converter = CreateDicomConverter(m_inputDicomDirectory,m_useBMatrixGradientDirections, m_transpose,
                                     m_smallGradientThreshold,m_allowLossyConversion);
  }
  else
  {
    std::cerr << "Invalid conversion mode" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}

int DWIConvert::write(const std::string& outputVolume)
{
  setOutputFileType(outputVolume);
  const std::string version = "4.8.0";
  if ("FSL" == getOutputFileType())
  {
    m_useIdentityMeasurementFrame = true; //Only true is valid for writing FSL
  }

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

  if( "FSL" == getOutputFileType())
  {

    Volume4DType::Pointer img4D = m_converter->OrientForFSLConventions();
    // write the image */
    m_converter->WriteFSLFormattedFileSet(outputVolumeHeaderName, m_outputBValues, m_outputBVectors, img4D);
  }
  else if ("Nrrd" == getOutputFileType())
  {
      //NRRD requires scaledDiffusionVectors, so find max BValue
      m_converter->ConvertToSingleBValueScaledDiffusionVectors();
      //////////////////////////////////////////////
      // write header file
      // This part follows a DWI NRRD file in NRRD format 5.
      // There should be a better way using itkNRRDImageIO.
      const std::string commentSection = m_converter->MakeFileComment(version, m_useBMatrixGradientDirections,
                                                                      m_useIdentityMeasurementFrame,
                                                                      m_smallGradientThreshold, getInputFileType());
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

int  DWIConvert::write(){
  if ("" != m_outputVolume){
      return write(m_outputVolume);
  }
  else{
    std::cout<<"illegal output Volume name. exit"<<std::endl;
    return EXIT_FAILURE;
  }
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

void DWIConvert::setInputFileType(const std::string& inputVolume, const std::string& inputDicomDirectory){
  m_inputVolume = inputVolume;
  m_inputDicomDirectory = inputDicomDirectory;
  if ("" == m_inputDicomDirectory){
    std::string inputExt = findFilenameExt(m_inputVolume);
    if (".nii" == inputExt) m_inputFileType = "FSL";
    if (".nrrd" == inputExt || ".nhdr" == inputExt) m_inputFileType = "Nrrd";
  }
  else if ("" == m_inputVolume)
  {
    m_inputFileType = "Dicom";
  }
  else{
    std::cerr <<"Error: One of inputVoume and inputDicomDirectory must be empty string"<<std::endl;
  }
}

void DWIConvert::setOutputFileType(const std::string& outputVolume){
  m_outputVolume = outputVolume;
  std::string outputExt = findFilenameExt(m_outputVolume);
  if (".nii" == outputExt) m_outputFileType = "FSL";
  else if (".nrrd" == outputExt || ".nhdr" == outputExt) m_outputFileType = "Nrrd";
  else{
    std::cerr <<"Error: the output file type is not supported currently"<<std::endl;
  }
}

std::string DWIConvert::getInputFileType()
{
  return m_inputFileType;
}

std::string DWIConvert::getOutputFileType()
{
  return m_outputFileType;
}

/*
void DWIConvert::setConversionMode(const std::string conversionMode){
    //["DicomToNrrd", "DicomToFSL", "NrrdToFSL", "FSLToNrrd","NrrdToNrrd", "FSLToFSL"]
    m_conversionMode =  conversionMode;
}

std::string DWIConvert::getConversionMode()
{
  if ("" == m_conversionMode) setConversionMode();
  return m_conversionMode;
}

 */

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

DWIConverter *DWIConvert::getConverter() const {
    return m_converter;
}
