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

const std::string emptyString(""); // A named empty string

DWIConvert::DWIConvert()
  : m_inputFileType{}
  , m_outputFileType{}
  , m_inputVolume{}
  , m_inputDicomDirectory{}
  , m_inputBValues{}
  , m_inputBVectors{}
  , m_gradientVectorFile{} // deprecated
  ,

  m_outputVolume{}
  , m_outputDirectory{ "." }
  , m_outputBValues{}
  , m_outputBVectors{}


{}

DWIConvert::DWIConvert(std::string inputVolume, std::string outputVolume)
  : DWIConvert()
{
  SetInputFileName(inputVolume);
  SetOutputFileName(outputVolume);
}

DWIConvert::~DWIConvert()
{
  delete m_converter;
}

int
DWIConvert::read()
{
  if (emptyString == getInputFileType())
  {
    std::cerr << "illegal input file type, exit" << std::endl;
    return EXIT_FAILURE;
  }

  // register DCMTK codecs, otherwise they will not be available when
  // `itkDCMTKSeriesFileNames` is used to build a list of filenames,
  // so reading series with JPEG transfer syntax will fail.
  DJDecoderRegistration::registerCodecs();
  DcmRLEDecoderRegistration::registerCodecs();

  if (m_fMRIOutput)
  {
    std::cerr << "Deprecated feature no longer supported: --fMRIOutput" << std::endl;
    return EXIT_FAILURE;
  }
  if (m_gradientVectorFile != emptyString)
  {
    std::cerr << "Deprecated feature no longer supported: --gradientVectorFile" << std::endl;
    return EXIT_FAILURE;
  }

  if (m_outputVolume == emptyString)
  {
    std::cerr << "Missing output volume name" << std::endl;
    return EXIT_FAILURE;
  }

  // build a NRRD file out of FSL output, which is two text files
  // for gradients and b values plus a NIfTI file for the gradient volumes.
  if ("FSL" == getInputFileType())
  {
    DWIConverter::FileNamesContainer filesList;
    filesList.clear();
    filesList.push_back(m_inputVolume);

    std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
    auto * FSLconverter = new FSLDWIConverter(filesList, m_inputBValues, m_inputBVectors);
    try
    {
      FSLconverter->SetAllowLossyConversion(m_allowLossyConversion);
      FSLconverter->LoadFromDisk();
      FSLconverter->SetThicknessFromSpacing(); // When converting from FSL, use thinkness as spacing[2]
      FSLconverter->ExtractDWIData();
    }
    catch (itk::ExceptionObject & excp)
    {
      itkGenericExceptionMacro(<< "Exception creating converter " << excp << std::endl);
      delete FSLconverter;
    }
    m_converter = FSLconverter;
  }
  // make FSL file set from a NRRD file.
  else if ("Nrrd" == getInputFileType())
  {
    DWIConverter::FileNamesContainer filesList;
    filesList.clear();
    filesList.push_back(m_inputVolume);

    std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
    auto * NRRDconverter = new NRRDDWIConverter(filesList);
    try
    {
      NRRDconverter->SetAllowLossyConversion(m_allowLossyConversion);
      NRRDconverter->LoadFromDisk();
      NRRDconverter->ExtractDWIData();
    }
    catch (itk::ExceptionObject & excp)
    {
      itkGenericExceptionMacro(<< "Exception creating converter " << excp << std::endl);
      delete NRRDconverter;
    }
    m_converter = NRRDconverter;
  }
  else if ("Dicom" == getInputFileType())
  {
    m_converter = CreateDicomConverter(
      m_inputDicomDirectory, m_useBMatrixGradientDirections, m_smallGradientThreshold, m_allowLossyConversion);
  }
  else
  {
    std::cerr << "Invalid conversion mode" << std::endl;
    return EXIT_FAILURE;
  }
  return (nullptr == m_converter ? EXIT_FAILURE : EXIT_SUCCESS);
}

int
DWIConvert::write(const std::string & outputVolume)
{
  SetOutputFileName(outputVolume);
  const std::string version = "5.0.0";
  if ("FSL" == getOutputFileType())
  {
    m_useIdentityMeasurementFrame = true; // Only true is valid for writing FSL
  }

  if (m_useIdentityMeasurementFrame)
  {
    m_converter->ConvertBVectorsToIdentityMeasurementFrame();
  }

  std::string outputVolumeHeaderName(m_outputVolume);
  { // concatenate with outputDirectory
    if (m_outputVolume.find('/') == std::string::npos && m_outputVolume.find('\\') == std::string::npos)
    {
      if (!outputVolumeHeaderName.empty())
      {
        outputVolumeHeaderName = m_outputDirectory;
        outputVolumeHeaderName += "/";
        outputVolumeHeaderName += m_outputVolume;
      }
    }
  }

  if ("FSL" == getOutputFileType())
  {
    m_converter->ConvertToMutipleBValuesUnitScaledBVectors();
    Volume4DType::Pointer img4D = m_converter->OrientForFSLConventions();
    // write the image */
    m_converter->WriteFSLFormattedFileSet(outputVolumeHeaderName, m_outputBValues, m_outputBVectors, img4D);
  }
  else if ("Nrrd" == getOutputFileType())
  {
    // NRRD requires scaledDiffusionVectors, so find max BValue
    m_converter->ConvertToSingleBValueScaledDiffusionVectors();
    //////////////////////////////////////////////
    // write header file
    // This part follows a DWI NRRD file in NRRD format 5.
    // There should be a better way using itkNRRDImageIO.
    const std::string commentSection = m_converter->MakeFileComment(version,
                                                                    m_useBMatrixGradientDirections,
                                                                    m_useIdentityMeasurementFrame,
                                                                    m_smallGradientThreshold,
                                                                    getInputFileType());
    m_converter->ManualWriteNRRDFile(outputVolumeHeaderName, commentSection);
  }
  else
  {
    std::cerr << "Invalidate conversion mode" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Wrote file: " << outputVolumeHeaderName << std::endl;
  return EXIT_SUCCESS;
}

int
DWIConvert::write()
{
  if (emptyString != m_outputVolume)
  {
    return write(m_outputVolume);
  }
  else
  {
    std::cout << "illegal output Volume name. exit" << std::endl;
    return EXIT_FAILURE;
  }
}

DWIConverter *
DWIConvert::CreateDicomConverter(const std::string inputDicomDirectory,
                                 const bool        useBMatrixGradientDirections,
                                 const double      smallGradientThreshold,
                                 const bool        allowLossyConversion)
{
  // check for required parameters
  if (inputDicomDirectory == emptyString)
  {
    std::cerr << "Missing DICOM input directory path" << std::endl;
    return nullptr;
  }

  // use the factor to instantiate a converter object based on the vender.
  DWIConverterFactory converterFactory(inputDicomDirectory, useBMatrixGradientDirections, smallGradientThreshold);
  DWIConverter *      converter;
  try
  {
    converter = converterFactory.New();
  }
  catch (itk::ExceptionObject & excp)
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
  catch (itk::ExceptionObject & excp)
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
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Exception extracting gradient vectors " << excp << std::endl;
    delete converter;
    return nullptr;
  }
  // this is a punt, it will still write out the volume image
  // even if we don't know how to extract gradients.
  if (converterFactory.GetVendor() == "GENERIC")
  {
    std::cerr << "Can't extract DWI data from files created by vendor " << converterFactory.GetVendor() << std::endl;
    delete converter;
    exit(EXIT_FAILURE);
  }
  return converter;
}


void
DWIConvert::SetInputFileName(std::string inputFilePath)
{

  const auto isDirectory = itksys::SystemTools::FileIsDirectory(inputFilePath);
  if (isDirectory)
  {
    m_inputDicomDirectory = inputFilePath;
  }
  else
  {
    const auto isRegularFile = itksys::SystemTools::FileExists(inputFilePath, true);
    if (isRegularFile)
    {
      m_inputVolume = inputFilePath;
    }
    else
    {
      std::cerr << "Error: the inputVolume is neither file nor directory." << std::endl;
    }
  }

  // prefer the inputVolume field if available
  if ((!m_inputVolume.empty()) || m_inputDicomDirectory.empty())
  {
    const std::string inputExt = itksys::SystemTools::GetFilenameExtension(m_inputVolume);
    if (std::string::npos != inputExt.rfind(".nii"))
    {
      m_inputFileType = "FSL";
    }
    else if (std::string::npos != inputExt.rfind(".nrrd") || std::string::npos != inputExt.rfind(".nhdr"))
    {
      m_inputFileType = "Nrrd";
    }
    else
    {
      std::cerr << "Error: file type of inputVoume is not supported currently" << std::endl;
    }
  }
  else if (emptyString != m_inputDicomDirectory)
  {
    m_inputFileType = "Dicom";
  }
  else
  {
    std::cerr << "Error: One of inputVoume and inputDicomDirectory must be empty string" << std::endl;
  }
}


void
DWIConvert::SetOutputFileName(std::string outputFilePath)
{
  m_outputVolume = outputFilePath;
  m_outputFileType = detectOuputVolumeType(outputFilePath);
}

std::string
DWIConvert::getInputFileType()
{
  return m_inputFileType;
}

std::string
DWIConvert::getOutputFileType()
{
  return m_outputFileType;
}


const std::string &
DWIConvert::getInputVolume() const
{
  return m_inputVolume;
}

void
DWIConvert::setInputVolume(const std::string & inputVolume)
{
  m_inputVolume = inputVolume;
}

const std::string &
DWIConvert::getInputDicomDirectory() const
{
  return m_inputDicomDirectory;
}

void
DWIConvert::setInputDicomDirectory(const std::string & inputDicomDirectory)
{
  m_inputDicomDirectory = inputDicomDirectory;
}

const std::string &
DWIConvert::getInputBValues() const
{
  return m_inputBValues;
}

void
DWIConvert::setInputBValues(const std::string & inputBValues)
{
  m_inputBValues = inputBValues;
}

const std::string &
DWIConvert::getInputBVectors() const
{
  return m_inputBVectors;
}

void
DWIConvert::setInputBVectors(const std::string & inputBVectors)
{
  m_inputBVectors = inputBVectors;
}

const std::string &
DWIConvert::getGradientVectorFile() const
{
  return m_gradientVectorFile;
}

void
DWIConvert::setGradientVectorFile(const std::string & gradientVectorFile)
{
  m_gradientVectorFile = gradientVectorFile;
}

double
DWIConvert::getSmallGradientThreshold() const
{
  return m_smallGradientThreshold;
}

void
DWIConvert::setSmallGradientThreshold(double smallGradientThreshold)
{
  m_smallGradientThreshold = smallGradientThreshold;
}

bool
DWIConvert::isfMRIOutput() const
{
  return m_fMRIOutput;
}

void
DWIConvert::setfMRIOutput(bool fMRIOutput)
{
  m_fMRIOutput = fMRIOutput;
}

bool
DWIConvert::isAllowLossyConversion() const
{
  return m_allowLossyConversion;
}

void
DWIConvert::setAllowLossyConversion(bool allowLossyConversion)
{
  m_allowLossyConversion = allowLossyConversion;
}

bool
DWIConvert::isUseIdentityMeasurementFrame() const
{
  return m_useIdentityMeasurementFrame;
}

void
DWIConvert::setUseIdentityMeasurementFrame(bool useIdentityMeasurementFrame)
{
  m_useIdentityMeasurementFrame = useIdentityMeasurementFrame;
}

bool
DWIConvert::isUseBMatrixGradientDirections() const
{
  return m_useBMatrixGradientDirections;
}

void
DWIConvert::setUseBMatrixGradientDirections(bool useBMatrixGradientDirections)
{
  m_useBMatrixGradientDirections = useBMatrixGradientDirections;
}

const std::string &
DWIConvert::getOutputVolume() const
{
  return m_outputVolume;
}

void
DWIConvert::setOutputVolume(const std::string & outputVolume)
{
  m_outputVolume = outputVolume;
}

const std::string &
DWIConvert::getOutputDirectory() const
{
  return m_outputDirectory;
}

void
DWIConvert::setOutputDirectory(const std::string & outputDirectory)
{
  m_outputDirectory = outputDirectory;
}

const std::string &
DWIConvert::getOutputBValues() const
{
  return m_outputBValues;
}

void
DWIConvert::setOutputBValues(const std::string & outputBValues)
{
  m_outputBValues = outputBValues;
}

const std::string &
DWIConvert::getOutputBVectors() const
{
  return m_outputBVectors;
}

void
DWIConvert::setOutputBVectors(const std::string & outputBVectors)
{
  m_outputBVectors = outputBVectors;
}

DWIConverter *
DWIConvert::getConverter() const
{
  return m_converter;
}


// public utility interface

std::string
detectOuputVolumeType(const std::string & outputVolume)
{
  std::string       outputFileType = "";
  const std::string outputExt = itksys::SystemTools::GetFilenameExtension(outputVolume);
  if (std::string::npos != outputExt.rfind(".nii"))
  {
    outputFileType = "FSL";
  }
  else if (std::string::npos != outputExt.rfind(".nrrd") || std::string::npos != outputExt.rfind(".nhdr"))
  {
    outputFileType = "Nrrd";
  }
  else
  {
    std::cerr << "Error: the output file type is not supported currently" << std::endl;
  }
  return outputFileType;
}

bool
convertInputVolumeVectorToNrrdOrNifti(const std::string &              targetType,
                                      const std::vector<std::string> & inputVolumeVector,
                                      std::vector<std::string> &       targetVolumeVector)
{
  targetVolumeVector.clear();
  int nSize = inputVolumeVector.size();
  for (int i = 0; i < nSize; ++i)
  {
    std::string targetVolume;
    if (convertInputVolumeToNrrdOrNifti(targetType, inputVolumeVector.at(i), targetVolume))
    {
      targetVolumeVector.push_back(targetVolume);
    }
    else
    {
      return false;
    }
  }
  return true;
}

bool
convertInputVolumeToNrrdOrNifti(const std::string & targetType,
                                const std::string & inputVolume,
                                std::string &       targetVolume)
{

  DWIConvert dwiConvert;
  dwiConvert.SetInputFileName(inputVolume);

  if (targetType == dwiConvert.getInputFileType())
  {
    targetVolume = inputVolume;
  }
  else
  {
    if ("FSL" == targetType)
    {
      targetVolume = itksys::SystemTools::GetFilenameWithoutExtension(inputVolume) + "_convert.nii";
    }
    else if ("Nrrd" == targetType)
    {
      targetVolume = itksys::SystemTools::GetFilenameWithoutExtension(inputVolume) + "_convert.nrrd";
    }
    else
    {
      std::cout << "Error: the targetType of DWIConvert is not supported. " << std::endl;
      return false;
    }

    dwiConvert.SetOutputFileName(targetVolume);
    int result = dwiConvert.read();
    if (EXIT_SUCCESS == result)
    {
      dwiConvert.write(targetVolume);
    }
    else
    {
      std::cout << "Error: read inputVolume failed." << std::endl;
      return false;
    }
  }
  return true;
}
