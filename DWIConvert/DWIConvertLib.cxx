//
// Created by Hui Xie on 12/26/16.
//

#include "DWIConvertLib.h"

#include "dcmtk/oflog/helpers/loglog.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmjpeg/djdecode.h"
#include "dcmtk/dcmjpls/djdecode.h"
#include "dcmtk/dcmdata/dcrledrg.h"
#include "itksys/SystemTools.hxx"

const std::string emptyString(""); //A named empty string

DWIConvert::DWIConvert()
{
    m_inputVolume = emptyString;
    m_inputDicomDirectory = emptyString;
    m_inputBValues = emptyString;  //default: emptyString
    m_inputBVectors = emptyString; //default: emptyString
    m_gradientVectorFile = emptyString; //deprecated
    m_smallGradientThreshold = 0.2; //default = 0.2


    m_inputFileType = emptyString;
    m_outputFileType = emptyString;

    m_fMRIOutput = false; //default: false
    m_transpose = false; //default:false
    m_allowLossyConversion = false; //defualt: false
    m_useIdentityMeasurementFrame = false; //default: false
    m_useBMatrixGradientDirections = false; //default: false

    m_outputVolume = emptyString;
    m_outputDirectory = ".";  //default: "."
    m_outputBValues = emptyString; //default: emptyString
    m_outputBVectors = emptyString;//default: emptyString

    m_converter = nullptr;

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
  if (emptyString == getInputFileType()){
    std::cerr << "illegal input file type, exit" << std::endl;
    return EXIT_FAILURE;
  }

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
  if( m_gradientVectorFile != emptyString ) {
    std::cerr << "Deprecated feature no longer supported: --gradientVectorFile" << std::endl;
    return EXIT_FAILURE;
  }

  if( m_outputVolume == emptyString )
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
  return (nullptr == m_converter ? EXIT_FAILURE : EXIT_SUCCESS);

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
    m_converter->ConvertToMutipleBValuesUnitScaledBVectors();
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
  if (emptyString != m_outputVolume){
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
  if( inputDicomDirectory == emptyString )
  {
    std::cerr << "Missing DICOM input directory path" << std::endl;
    return nullptr;
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
    return nullptr;
  }
  if (nullptr == converter)
  {
    std::cerr << "Unable to create converter!" << std::endl;
    return nullptr;
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
    return nullptr;
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
    return nullptr;
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


//currently supported file types: { ".nii", ".nii.gz", ".nhdr", ".nrrd"}
void DWIConvert::setInputFileType(const std::string& inputVolume, const std::string& inputDicomDirectory){
  m_inputVolume = inputVolume;
  m_inputDicomDirectory = inputDicomDirectory;

  // prefer the inputVolume field if available
  if ( (!m_inputVolume.empty()) ||
       m_inputDicomDirectory.empty() )
  {
    const std::string inputExt = itksys::SystemTools::GetFilenameExtension(m_inputVolume);
    if ( std::string::npos != inputExt.rfind(".nii"))
    {
      m_inputFileType = "FSL";
    }
    else if (std::string::npos != inputExt.rfind(".nrrd") || std::string::npos != inputExt.rfind(".nhdr"))
    {
      m_inputFileType = "Nrrd";
    }
    else
    {
      std::cerr <<"Error: file type of inputVoume is not supported currently"<<std::endl;
    }
  }
  else if (emptyString != m_inputDicomDirectory)
  {
    m_inputFileType = "Dicom";
  }
  else{
    std::cerr <<"Error: One of inputVoume and inputDicomDirectory must be empty string"<<std::endl;
  }
}

std::string DWIConvert::detectOuputVolumeType(const std::string& outputVolume){
  std::string outputFileType;
  const std::string outputExt = itksys::SystemTools::GetFilenameExtension(outputVolume);
  if (std::string::npos != outputExt.rfind(".nii")) outputFileType = "FSL";
  else if (std::string::npos != outputExt.rfind(".nrrd")|| std::string::npos != outputExt.rfind(".nhdr")) outputFileType = "Nrrd";
  else{
    std::cerr <<"Error: the output file type is not supported currently"<<std::endl;
  }
  return outputFileType;

}


void DWIConvert::setOutputFileType(const std::string& outputVolume){
  m_outputVolume = outputVolume;
  m_outputFileType = detectOuputVolumeType(outputVolume);
}

std::string DWIConvert::getInputFileType()
{
  return m_inputFileType;
}

std::string DWIConvert::getOutputFileType()
{
  return m_outputFileType;
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

int DWIConvert::convertInputVolumeVectorToNrrdOrNifti(const std::string targetType,
                                          const std::vector<std::string> inputVolumeVector,
                                          std::vector<std::string>& targetVolumeVector){
  targetVolumeVector.clear();
  int nSize = inputVolumeVector.size();
  for (int i = 0; i< nSize; ++i){
    std::string targetVolume;
    if (0 == convertInputVolumeToNrrdOrNifti(targetType, inputVolumeVector.at(i),targetVolume)){
      targetVolumeVector.push_back(targetVolume);
    } else{
      return 1;
     }
  }
  return 0;
}

int DWIConvert::convertInputVolumeToNrrdOrNifti(const std::string targetType,
                                                const std::string inputVolume, std::string &targetVolume) {

  // AVAILABLE_MAC_OS_X_VERSION_10_0_AND_LATER
  // On Windows,judge file or directory are different with Linux. It does not support Windows currently.
  struct stat inputInfor;
  if (-1 == stat(inputVolume.c_str(), &inputInfor)) {
    std::cout << "Error: inputVolume is illegal file description. " << std::endl;
    return 1;
  } else {
    if (inputInfor.st_mode & S_IFDIR) {
      setInputFileType("", inputVolume);
    }
    else if (inputInfor.st_mode & S_IFREG) {
      setInputFileType(inputVolume, "");
    } else {
      std::cout << "Error: the inputVolume is neither file nor directory." << std::endl;
      return -1;
    }
  }

  if (targetType == getInputFileType()) {
    targetVolume = inputVolume;
  }
  else {
    if ("FSL" == targetType){
      targetVolume = itksys::SystemTools::GetFilenameWithoutExtension(inputVolume)+"_convert.nii";
    }
    else if("Nrrd" == targetType){
      targetVolume = itksys::SystemTools::GetFilenameWithoutExtension(inputVolume)+"_convert.nrrd";
    }
    else{
      std::cout <<"Error: the targetType of DWIConvert is not supported. "<<std::endl;
      return -1;
    }

    setOutputFileType(targetVolume);
    int result = read();
    if (EXIT_SUCCESS == result) {
      write(targetVolume);
    }
    else {
      std::cout<<"Error: read inputVolume failed."<<std::endl;
      return 1;
    }
  }
  return 0;
}