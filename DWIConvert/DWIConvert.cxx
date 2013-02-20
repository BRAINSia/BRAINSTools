/*=========================================================================
Computing (NAMIC), funded by the National Institutes of Health
through the NIH Roadmap for Medical Research, Grant U54 EB005149.

See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

 ***
 This program converts Diffusion weighted MR images in Dicom format into
 NRRD format.

Assumptions:

1) Uses left-posterior-superior (Dicom default) as default space for philips and siemens.
This is the default space for NRRD header.
2) For GE data, Dicom data are arranged in volume interleaving order.
3) For Siemens data, images are arranged in mosaic form.
4) For oblique collected Philips data, the measurement frame for the
gradient directions is the same as the ImageOrientationPatient

Reference materials:
DICOM Data Dictionary: http://medical.nema.org/Dicom/2011/11_06pu.pdf
=========================================================================*/
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include "DWIConvertCLP.h"

#include "itkMacro.h"
#include "itkIntTypes.h"
#include "itkDCMTKSeriesFileNames.h"
#undef HAVE_SSTREAM
#include "itkDCMTKImageIO.h"
#include "itkRawImageIO.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itksys/Directory.hxx"
#include "itksys/SystemTools.hxx"
#include "itksys/Base64.h"
#undef HAVE_SSTREAM
#include "itkDCMTKFileReader.h"
#include "djdecode.h"
#include "StringContains.h"
#include "DWIConvertUtils.h"
#include "DoubleToString.h"

unsigned int ConvertFromCharPtr(const char *s)
{
  unsigned int rval = 0;

  // assume little-endian
  for( unsigned i = 0; i < sizeof(unsigned int); ++i )
    {
    rval += ( (unsigned int)s[i]) << (i * 8);
    }
// this makes no sense, according to what I've read, but apparently,
// the uint32 numbers in the CSA header are little-endian even in
// files with BigEndian transfer syntax.#if 0
#if 0

  switch( xferSyntax )
    {
    case EXS_LittleEndianImplicit:
    case EXS_LittleEndianExplicit:
      {
      itk::ByteSwapper<unsigned int>::SwapFromSystemToLittleEndian(&rval);
      }
      break;
    case EXS_BigEndianImplicit:
    case EXS_BigEndianExplicit:
      {
      itk::ByteSwapper<unsigned int>::SwapFromSystemToBigEndian(&rval);
      }
      break;
    default:
      {
      }
      break;
    }
#endif
  return rval;
}

/** pull data out of Siemens scans.
 *
 *  Siemens sticks most of the DTI information into a single
 *  OB-format entry.  This is actually rigidly structured, but
 *  this function depends on the needed data living at fixed offset
 *  from the beginning of the name of each tag, and ignores the
 *  internal structure documented in the Siemens Dicom Compliance
 *  document.
 */
unsigned int
ExtractSiemensDiffusionInformation(const std::string & tagString,
                                   const std::string & nameString,
                                   std::vector<double>& valueArray)
{
  ::size_t atPosition = tagString.find( nameString );

  if( atPosition == std::string::npos )
    {
    return 0;
    }
  while( true )  // skip nameString inside a quotation
    {
    std::string nextChar = tagString.substr( atPosition + nameString.size(), 1 );

    if( nextChar.c_str()[0] == 0 )
      {
      break;
      }
    else
      {
      atPosition = tagString.find( nameString, atPosition + 2 );
      }
    }

  if( atPosition == std::string::npos )
    {
    return 0;
    }
  std::string  infoAsString = tagString.substr( atPosition, tagString.size() - atPosition + 1 );
  const char * infoAsCharPtr = infoAsString.c_str();

  unsigned int vm = ConvertFromCharPtr(infoAsCharPtr + 64);
    {
    std::string vr = infoAsString.substr( 68, 2 );
    int         syngodt = ConvertFromCharPtr(infoAsCharPtr + 72);
    int         nItems = ConvertFromCharPtr(infoAsCharPtr + 76);
    int         localDummy = ConvertFromCharPtr(infoAsCharPtr + 80);

    // std::cout << "\tName String: " << nameString << std::endl;
    // std::cout << "\tVR: " << vr << std::endl;
    // std::cout << "\tVM: " << vm << std::endl;
    // std::cout << "Local String: " << infoAsString.substr(0,80) << std::endl;

    /* This hack is required for some Siemens VB15 Data */
    if( ( nameString == "DiffusionGradientDirection" ) && (vr != "FD") )
      {
      bool loop = true;
      while( loop )
        {
        atPosition = tagString.find( nameString, atPosition + 26 );
        if( atPosition == std::string::npos )
          {
          // std::cout << "\tFailed to find DiffusionGradientDirection Tag - returning" << vm << std::endl;
          return 0;
          }
        infoAsString = tagString.substr( atPosition, tagString.size() - atPosition + 1 );
        infoAsCharPtr = infoAsString.c_str();
        // std::cout << "\tOffset to new position" << std::endl;
        // std::cout << "\tNew Local String: " << infoAsString.substr(0,80) << std::endl;
        vm = ConvertFromCharPtr(infoAsCharPtr + 64);
        vr = infoAsString.substr( 68, 2 );
        if( vr == "FD" )
          {
          loop = false;
          }
        syngodt = ConvertFromCharPtr(infoAsCharPtr + 72);
        nItems = ConvertFromCharPtr(infoAsCharPtr + 76);
        localDummy = ConvertFromCharPtr(infoAsCharPtr + 80);
        // std::cout << "\tVR: " << vr << std::endl;
        // std::cout << "\tVM: " << vm << std::endl;
        }
      }
    else
      {
      // std::cout << "\tUsing initial position" << std::endl;
      }
    // std::cout << "\tArray Length: " << vm << std::endl;
    }

  unsigned int offset = 84;
  for( unsigned int k = 0; k < vm; ++k )
    {
    const int itemLength = ConvertFromCharPtr(infoAsCharPtr + offset + 4);
    const int strideSize = static_cast<int>(ceil(static_cast<double>(itemLength) / 4) * 4);
    if( infoAsString.length() < offset + 16 + itemLength )
      {
      // data not available or incomplete
      return 0;
      }
    const std::string valueString = infoAsString.substr( offset + 16, itemLength );
    valueArray.push_back( atof(valueString.c_str() ) );
    offset += 16 + strideSize;
    }
  return vm;
}

/**
 *  Add private tags to the Dicom Dictionary
 */
void
AddFlagsToDictionary()
{
  // these have to be dynamically allocated because otherwise there's
  // a malloc error after main exits.

  // relevant GE tags
  DcmDictEntry *GEDictBValue = new DcmDictEntry(0x0043, 0x1039, DcmVR(EVR_IS),
                                                "B Value of diffusion weighting", 1, 1, 0, true,
                                                "dicomtonrrd");
  DcmDictEntry *GEDictXGradient = new DcmDictEntry(0x0019, 0x10bb, DcmVR(EVR_DS),
                                                   "X component of gradient direction", 1, 1, 0, true,
                                                   "dicomtonrrd");
  DcmDictEntry *GEDictYGradient = new DcmDictEntry(0x0019, 0x10bc, DcmVR(EVR_DS),
                                                   "Y component of gradient direction", 1, 1, 0, true,
                                                   "dicomtonrrd");
  DcmDictEntry *GEDictZGradient = new DcmDictEntry(0x0019, 0x10bd, DcmVR(EVR_DS),
                                                   "Z component of gradient direction", 1, 1, 0, true,
                                                   "dicomtonrrd");

  // relevant Siemens private tags
  DcmDictEntry *SiemensMosiacParameters = new DcmDictEntry(0x0051, 0x100b, DcmVR(EVR_IS),
                                                           "Mosiac Matrix Size", 1, 1, 0, true,
                                                           "dicomtonrrd");
  DcmDictEntry *SiemensDictNMosiac = new DcmDictEntry(0x0019, 0x100a, DcmVR(EVR_US),
                                                      "Number of Images In Mosaic", 1, 1, 0, true,
                                                      "dicomtonrrd");
  DcmDictEntry *SiemensDictBValue = new DcmDictEntry(0x0019, 0x100c, DcmVR(EVR_IS),
                                                     "B Value of diffusion weighting", 1, 1, 0, true,
                                                     "dicomtonrrd");
  DcmDictEntry *SiemensDictDiffusionDirection = new DcmDictEntry(0x0019, 0x100e, DcmVR(EVR_FD),
                                                                 "Diffusion Gradient Direction", 3, 3, 0, true,
                                                                 "dicomtonrrd");
  DcmDictEntry *SiemensDictDiffusionMatrix = new DcmDictEntry(0x0019, 0x1027, DcmVR(EVR_FD),
                                                              "Diffusion Matrix", 6, 6, 0, true,
                                                              "dicomtonrrd");
  DcmDictEntry *SiemensDictShadowInfo = new DcmDictEntry(0x0029, 0x1010, DcmVR(EVR_OB),
                                                         "Siemens DWI Info", 1, 1, 0, true,
                                                         "dicomtonrrd");

  // relevant Philips private tags
  DcmDictEntry *PhilipsDictBValue  = new DcmDictEntry(0x2001, 0x1003, DcmVR(EVR_FL),
                                                      "B Value of diffusion weighting", 1, 1, 0, true,
                                                      "dicomtonrrd");
  DcmDictEntry *PhilipsDictDiffusionDirection   = new DcmDictEntry(0x2001, 0x1004, DcmVR(EVR_CS),
                                                                   "Diffusion Gradient Direction", 1, 1, 0, true,
                                                                   "dicomtonrrd");
  DcmDictEntry *PhilipsDictDiffusionDirectionRL = new DcmDictEntry(0x2005, 0x10b0, DcmVR(EVR_FL),
                                                                   "Diffusion Direction R/L", 4, 4, 0, true,
                                                                   "dicomtonrrd");
  DcmDictEntry *PhilipsDictDiffusionDirectionAP = new DcmDictEntry(0x2005, 0x10b1, DcmVR(EVR_FL),
                                                                   "Diffusion Direction A/P", 4, 4, 0, true,
                                                                   "dicomtonrrd");
  DcmDictEntry *PhilipsDictDiffusionDirectionFH = new DcmDictEntry(0x2005, 0x10b2, DcmVR(EVR_FL),
                                                                   "Diffusion Direction F/H", 4, 4, 0, true,
                                                                   "dicomtonrrd");

  itk::DCMTKFileReader::AddDictEntry(GEDictBValue);
  itk::DCMTKFileReader::AddDictEntry(GEDictXGradient);
  itk::DCMTKFileReader::AddDictEntry(GEDictYGradient);
  itk::DCMTKFileReader::AddDictEntry(GEDictZGradient);

  // relevant Siemens private tags
  itk::DCMTKFileReader::AddDictEntry(SiemensMosiacParameters);
  itk::DCMTKFileReader::AddDictEntry(SiemensDictNMosiac);
  itk::DCMTKFileReader::AddDictEntry(SiemensDictBValue);
  itk::DCMTKFileReader::AddDictEntry(SiemensDictDiffusionDirection);
  itk::DCMTKFileReader::AddDictEntry(SiemensDictDiffusionMatrix);
  itk::DCMTKFileReader::AddDictEntry(SiemensDictShadowInfo);

  // relevant Philips private tags
  itk::DCMTKFileReader::AddDictEntry(PhilipsDictBValue);
  itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirection);
  itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionRL);
  itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionAP);
  itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionFH);
}

/** Free the headers for all dicom files.
 *  also calls the DJDecoder cleanup.
 */
void FreeHeaders(std::vector<itk::DCMTKFileReader *> & allHeaders)
{
  for( std::vector<itk::DCMTKFileReader *>::iterator it = allHeaders.begin();
       it != allHeaders.end(); ++it )
    {
    delete (*it);
    }
}

typedef short                         PixelValueType;
typedef itk::Image<PixelValueType, 3> VolumeType;

int
Write4DVolume( VolumeType::Pointer & img, int nVolumes, const std::string & fname )
{
  typedef itk::Image<PixelValueType, 4> Volume4DType;

  VolumeType::SizeType      size3D(img->GetLargestPossibleRegion().GetSize() );
  VolumeType::DirectionType direction3D(img->GetDirection() );
  VolumeType::SpacingType   spacing3D(img->GetSpacing() );
  VolumeType::PointType     origin3D(img->GetOrigin() );

  Volume4DType::SizeType size4D;
  size4D[0] = size3D[0];
  size4D[1] = size3D[1];
  size4D[2] = size3D[2] / nVolumes;
  size4D[3] = nVolumes;

  if( (size4D[2] * nVolumes) != size3D[2] )
    {
    std::cerr << "#of slices in volume not evenly divisible by"
              << " the number of volumes: slices = " << size3D[2]
              << " volumes = " << nVolumes << " left-over slices = "
              << size3D[2] % nVolumes << std::endl;
    }
  Volume4DType::DirectionType direction4D;
  Volume4DType::SpacingType   spacing4D;
  Volume4DType::PointType     origin4D;
  for( unsigned i = 0; i < 3; ++i )
    {
    for( unsigned j = 0; j < 3; ++j )
      {
      direction4D[i][j] = direction3D[i][j];
      }
    direction4D[3][i] = 0.0;
    direction4D[i][3] = 0.0;
    spacing4D[i] = spacing3D[i];
    origin4D[i] = origin3D[i];
    }
  direction4D[3][3] = 1.0;
  spacing4D[3] = 1.0;
  origin4D[3] = 0.0;

  Volume4DType::Pointer img4D = Volume4DType::New();
  img4D->SetRegions(size4D);
  img4D->SetDirection(direction4D);
  img4D->SetSpacing(spacing4D);
  img4D->SetOrigin(origin4D);

  img4D->Allocate();
  size_t bytecount = img4D->GetLargestPossibleRegion().GetNumberOfPixels();
  bytecount *= sizeof(PixelValueType);
  memcpy(img4D->GetBufferPointer(), img->GetBufferPointer(), bytecount);
#if 0
    {
    itk::ImageFileWriter<VolumeType>::Pointer writer = itk::ImageFileWriter<VolumeType>::New();
    writer->SetFileName( "dwi3dconvert.nii.gz");
    writer->SetInput( img );
    writer->Update();
    }
#endif

  itk::ImageFileWriter<Volume4DType>::Pointer imgWriter =
    itk::ImageFileWriter<Volume4DType>::New();

  imgWriter->SetInput( img4D );
  imgWriter->SetFileName( fname.c_str() );
  try
    {
    imgWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing "
              << fname << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

void
DeInterleaveVolume(VolumeType::Pointer & volume,
                   size_t SlicesPerVolume,
                   size_t NSlices)
{
  size_t NVolumes = NSlices / SlicesPerVolume;

  VolumeType::RegionType R = volume->GetLargestPossibleRegion();

  R.SetSize(2, 1);
  std::vector<VolumeType::PixelType> v(NSlices);
  std::vector<VolumeType::PixelType> w(NSlices);

  itk::ImageRegionIteratorWithIndex<VolumeType> I( volume, R );
  // permute the slices by extracting the 1D array of voxels for
  // a particular {x,y} position, then re-ordering the voxels such
  // that all the voxels for a particular volume are adjacent
  for( I.GoToBegin(); !I.IsAtEnd(); ++I )
    {
    VolumeType::IndexType idx = I.GetIndex();
    // extract all values in one "column"
    for( unsigned int k = 0; k < NSlices; ++k )
      {
      idx[2] = k;
      v[k] = volume->GetPixel( idx );
      }
    // permute
    for( unsigned int k = 0; k < NVolumes; ++k )
      {
      for( unsigned int m = 0; m < SlicesPerVolume; ++m )
        {
        w[(k * SlicesPerVolume) + m] = v[(m * NVolumes) + k];
        }
      }
    // put things back in order
    for( unsigned int k = 0; k < NSlices; ++k )
      {
      idx[2] = k;
      volume->SetPixel( idx, w[k] );
      }
    }
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  // just need one instance to do double to string conversions
  DoubleToString DoubleConvert;

  if( conversionMode == "FSLToNrrd" )
    {
    extern int FSLToNrrd(const std::string & inputVolume,
                         const std::string & outputVolume,
                         const std::string & inputBValues,
                         const std::string & inputBVectors);

    return FSLToNrrd(inputVolume, outputVolume,
                     inputBValues, inputBVectors);
    }
  if( conversionMode == "NrrdToFSL" )
    {
    extern int NrrdToFSL(const std::string & inputVolume,
                         const std::string & outputVolume,
                         const std::string & outputBValues,
                         const std::string & outputBVectors);
    return NrrdToFSL(inputVolume, outputVolume,
                     outputBValues, outputBVectors);
    }

  typedef itk::ImageSeriesReader<VolumeType> ReaderType;
  typedef itk::ImageFileReader<VolumeType>   SingleFileReaderType;
  typedef itk::DCMTKSeriesFileNames          InputNamesGeneratorType;

  AddFlagsToDictionary();
  bool nrrdFormat(true);
  // check for required parameters
  if( inputDicomDirectory == "" )
    {
    std::cerr << "Missing DICOM input directory path" << std::endl;
    return EXIT_FAILURE;
    }

  if( outputVolume == "" )
    {
    std::cerr << "Missing DICOM output volume name" << std::endl;
    return EXIT_FAILURE;
    }

  std::string outputVolumeHeaderName(outputVolume);
  if( outputVolume.find("/") == std::string::npos &&
      outputVolume.find("\\") == std::string::npos )
    {
    if( outputVolumeHeaderName.size() != 0 )
      {
      outputVolumeHeaderName = outputDirectory;
      outputVolumeHeaderName += "/";
      outputVolumeHeaderName += outputVolume;
      }
    }

  // decide whether the output is a single file or
  // header/raw pair
  std::string outputVolumeDataName;
  std::string outputFSLBValFilename;
  std::string outputFSLBVecFilename;
  if( conversionMode != "DicomToFSL" )
    {
    const size_t extensionPos = outputVolumeHeaderName.find(".nhdr");
    if( extensionPos != std::string::npos )
      {
      outputVolumeDataName = outputVolumeHeaderName.substr(0, extensionPos);
      outputVolumeDataName += ".raw";
      nrrdFormat = false;
      }
    }
  else
    {
    // FSL output of gradients & BValues
    std::string fslPrefix;
    size_t      extensionPos;
    extensionPos = outputVolumeHeaderName.find(".nii.gz");
    if( extensionPos == std::string::npos )
      {
      extensionPos = outputVolumeHeaderName.find(".nii");
      if( extensionPos == std::string::npos )
        {
        std::cerr << "FSL Format output chosen, "
                  << "but output Volume not a recognized "
                  << "NIfTI filename " << outputVolumeHeaderName
                  << std::endl;
        exit(1);
        }
      }
    if( outputBValues == "" )
      {
      outputFSLBValFilename = outputVolumeHeaderName.substr(0, extensionPos);
      outputFSLBValFilename += ".bval";
      }
    else
      {
      outputFSLBValFilename = outputBValues;
      }
    if( outputBVectors == "" )
      {
      outputFSLBVecFilename = outputVolumeHeaderName.substr(0, extensionPos);
      outputFSLBVecFilename += ".bvec";
      }
    else
      {
      outputFSLBVecFilename = outputBVectors;
      }
    }

  ReaderType::FileNamesContainer inputFileNames;
  // get the names of all slices in the directory
  InputNamesGeneratorType::Pointer inputNames = InputNamesGeneratorType::New();
  if( itksys::SystemTools::FileIsDirectory(inputDicomDirectory.c_str() ) )
    {
    inputNames->SetUseSeriesDetails( true);
    inputNames->SetLoadSequences( true );
    inputNames->SetLoadPrivateTags( true );
    inputNames->SetInputDirectory(inputDicomDirectory);
    inputFileNames = inputNames->GetInputFileNames();
    }
  else if( itksys::SystemTools::FileExists(inputDicomDirectory.c_str() ) )
  // or, if it isn't a directory, maybe it is
  // a multi-frame DICOM file.
    {
    inputFileNames.push_back(inputDicomDirectory);
    }

  const size_t nFiles = inputFileNames.size();

  if( nFiles < 1 )
    {
    std::cerr << "Error: no DICOMfiles found in inputDirectory: " << inputDicomDirectory
              << std::endl;
    return EXIT_FAILURE;
    }

  //////////////////////////////////////////////////
  // load all files in the dicom series.
  //////////////////////////////////////////////////

  std::vector<itk::DCMTKFileReader *> allHeaders(inputFileNames.size() );
  int                                 headerCount = 0;
  for( unsigned i = 0; i < allHeaders.size(); ++i )
    {
    itk::DCMTKFileReader *curReader = new itk::DCMTKFileReader;
    curReader->SetFileName(inputFileNames[i]);
    try
      {
      curReader->LoadFile();
      }
    catch( ... )
      {
      std::cerr << "Error reading slice" << inputFileNames[i] << std::endl;
      FreeHeaders(allHeaders);
      return EXIT_FAILURE;
      }
    // check for pixel data.
    if( !curReader->HasPixelData() )
      {
      delete curReader;
      }
    else
      {
      allHeaders[headerCount] = curReader;
      headerCount++;
      }
    }

  if( headerCount == 0 )
    {
    std::cerr << "No pixel data in series" << std::endl;
    return EXIT_FAILURE;
    }
  // reorder the filename list
  inputFileNames.resize( 0 );
  for( unsigned i = 0; i < allHeaders.size(); ++i )
    {
    std::string fname = allHeaders[i]->GetFileName();
    inputFileNames.push_back(fname);
    }

  // find vendor and modality
  std::string vendor;
  try
    {
    allHeaders[0]->GetElementLO(0x0008, 0x0070, vendor);
    strupper(vendor);
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Can't get vendor name from DICOM file" << excp << std::endl;
    FreeHeaders(allHeaders);
    return EXIT_FAILURE;
    }

  std::string modality;
  try
    {
    allHeaders[0]->GetElementCS(0x0008, 0x0060, modality);
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Can't find modality in DICOM file" << excp << std::endl;
    FreeHeaders(allHeaders);
    return EXIT_FAILURE;
    }

  // IF it's a PET or SPECT file, just write it out as a float image.
  if( StringContains(modality, "PT") || StringContains(modality, "ST") )
    {
    typedef itk::Image<float, 3> USVolumeType;
    itk::ImageSeriesReader<USVolumeType>::Pointer seriesReader =
      itk::ImageSeriesReader<USVolumeType>::New();
    seriesReader->SetFileNames( inputFileNames );

    itk::ImageFileWriter<USVolumeType>::Pointer nrrdImageWriter =
      itk::ImageFileWriter<USVolumeType>::New();

    nrrdImageWriter->SetFileName( outputVolumeHeaderName );
    nrrdImageWriter->SetInput( seriesReader->GetOutput() );
    try
      {
      nrrdImageWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      FreeHeaders(allHeaders);
      return EXIT_FAILURE;
      }
    FreeHeaders(allHeaders);
    return EXIT_SUCCESS;
    }

  // any failure in extracting data is an error
  // so encapsulate rest of program in exception block
  try
    {
    bool SliceMosaic(false);
    // decide if it is a mosaic
    if( StringContains(vendor, "SIEMENS") )
      {
      std::string ImageType;
      allHeaders[0]->GetElementCS(0x0008, 0x0008, ImageType);
      if( ImageType.find("MOSAIC") != std::string::npos )
        {
        SliceMosaic = true;
        }
      }
    else if( !StringContains(vendor, "GE") && !StringContains(vendor, "PHILIPS") )
      {
      std::cerr << "Unrecognized scanner vendor |"
                << vendor << "|" << std::endl;
      }

    //////////////////////////////////////////////////
    // 1) Read the input series as an array of slices
    unsigned int               nSlice;
    VolumeType::Pointer        readerOutput;
    itk::DCMTKImageIO::Pointer dcmtkIO = itk::DCMTKImageIO::New();
    bool                       multiSliceVolume;
    if( inputFileNames.size() > 1 )
      {
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetImageIO( dcmtkIO );
      reader->SetFileNames( inputFileNames );
      nSlice
        = inputFileNames.size();
      try
        {
        reader->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cerr << "Exception thrown while reading the series" << std::endl;
        std::cerr << excp << std::endl;
        FreeHeaders(allHeaders);
        return EXIT_FAILURE;
        }
      readerOutput = reader->GetOutput();
      multiSliceVolume = false;
      }
    else
      {
      SingleFileReaderType::Pointer reader =
        SingleFileReaderType::New();
      reader->SetImageIO( dcmtkIO );
      reader->SetFileName( inputFileNames[0] );
      nSlice = inputFileNames.size();
      try
        {
        reader->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cerr << "Exception thrown while reading the series" << std::endl;
        std::cerr << excp << std::endl;
        FreeHeaders(allHeaders);
        return EXIT_FAILURE;
        }
      readerOutput = reader->GetOutput();
      multiSliceVolume = true;
      }

    // get image dims and resolution
    unsigned short nRows, nCols;

    allHeaders[0]->GetElementUS(0x0028, 0x0010, nRows);
    allHeaders[0]->GetElementUS(0x0028, 0x0011, nCols);

    double sliceSpacing;
    double xRes, yRes;
      {
      double spacing[3];
      allHeaders[0]->GetSpacing(spacing);
      yRes = spacing[1];
      xRes = spacing[0];
      sliceSpacing = spacing[2];
      }

    itk::Vector<double, 3> ImageOrigin;
      {
      double origin[3];
      allHeaders[0]->GetOrigin(origin);
      ImageOrigin[0] = origin[0];
      ImageOrigin[1] = origin[1];
      ImageOrigin[2] = origin[2];
      }

    unsigned int               numberOfSlicesPerVolume = 0;
    std::map<std::string, int> sliceLocations;
    if( !multiSliceVolume )
      {
      // Make a hash of the sliceLocations in order to get the correct
      // count.  This is more reliable since SliceLocation may not be available.
      std::vector<int>         sliceLocationIndicator;
      std::vector<std::string> sliceLocationStrings;

      sliceLocationIndicator.resize( nSlice );
      for( unsigned int k = 0; k < nSlice; ++k )
        {
        std::string originString;

        allHeaders[k]->GetElementDS(0x0020, 0x0032, originString );
        sliceLocationStrings.push_back( originString );
        sliceLocations[originString]++;
        // std::cerr << inputFileNames[k] << " " << originString << std::endl;
        }
      for( unsigned int k = 0; k < nSlice; ++k )
        {
        std::map<std::string, int>::iterator it = sliceLocations.find( sliceLocationStrings[k] );
        sliceLocationIndicator[k] = distance( sliceLocations.begin(), it );
        }

      numberOfSlicesPerVolume = sliceLocations.size();
      std::cout << "=================== numberOfSlicesPerVolume:" << numberOfSlicesPerVolume << std::endl;

      // if the numberOfSlicesPerVolume == 1, de-interleaving won't do
      // anything so there's no point in doing it.
      if( nSlice >= 2 && numberOfSlicesPerVolume > 1 )
        {
        if( sliceLocationIndicator[0] != sliceLocationIndicator[1] )
          {
          std::cout << "Dicom images are ordered in a volume interleaving way." << std::endl;
          }
        else
          {
          std::cout << "Dicom images are ordered in a slice interleaving way." << std::endl;
          // reorder slices into a volume interleaving manner
          DeInterleaveVolume(readerOutput, numberOfSlicesPerVolume, nSlice);
#if 0
          itk::ImageFileWriter<VolumeType>::Pointer writer = itk::ImageFileWriter<VolumeType>::New();
          writer->SetFileName( "deinterleave.nrrd");
          writer->SetInput( readerOutput );
          writer->Update();
#endif
          }
        }
      }
    itk::Matrix<double, 3, 3> MeasurementFrame;
    MeasurementFrame.SetIdentity();

    // check ImageOrientationPatient and figure out slice direction in
    // L-P-I (right-handed) system.
    // In Dicom, the coordinate frame is L-P by default. Look at
    // http://medical.nema.org/dicom/2007/07_03pu.pdf ,  page 301
    itk::Matrix<double, 3, 3> LPSDirCos;

      {
      double dirCosArray[6];
      // 0020,0037 -- Image Orientation (Patient)
      allHeaders[0]->GetDirCosArray(dirCosArray);
      double *dirCosArrayP = dirCosArray;
      for( unsigned i = 0; i < 2; ++i )
        {
        for( unsigned j = 0; j < 3; ++j, ++dirCosArrayP )
          {
          LPSDirCos[j][i] = *dirCosArrayP;
          }
        }
      }

    // Cross product, this gives I-axis direction
    LPSDirCos[0][2] = (LPSDirCos[1][0] * LPSDirCos[2][1] - LPSDirCos[2][0] * LPSDirCos[1][1]);
    LPSDirCos[1][2] = (LPSDirCos[2][0] * LPSDirCos[0][1] - LPSDirCos[0][0] * LPSDirCos[2][1]);
    LPSDirCos[2][2] = (LPSDirCos[0][0] * LPSDirCos[1][1] - LPSDirCos[1][0] * LPSDirCos[0][1]);

    std::cout << "ImageOrientationPatient (0020:0037): ";
    std::cout << "LPS Orientation Matrix" << std::endl;
    std::cout << LPSDirCos << std::endl;

    itk::Matrix<double, 3, 3> SpacingMatrix;
    SpacingMatrix.Fill(0.0);
    SpacingMatrix[0][0] = xRes;
    SpacingMatrix[1][1] = yRes;
    SpacingMatrix[2][2] = sliceSpacing;
    std::cout << "SpacingMatrix" << std::endl;
    std::cout << SpacingMatrix << std::endl;

    itk::Matrix<double, 3, 3> OrientationMatrix;
    OrientationMatrix.SetIdentity();

    itk::Matrix<double, 3, 3> NRRDSpaceDirection;
    std::string               nrrdSpaceDefinition = "left-posterior-superior";;
    NRRDSpaceDirection = LPSDirCos * OrientationMatrix * SpacingMatrix;

    std::cout << "NRRDSpaceDirection" << std::endl;
    std::cout << NRRDSpaceDirection << std::endl;

    unsigned int mMosaic = 0;   // number of raws in each mosaic block;
    unsigned int nMosaic = 0;   // number of columns in each mosaic block
    unsigned int nSliceInVolume = 0;
    unsigned int nVolume = 0;
    bool         SliceOrderIS(true);

    // figure out slice order and mosaic arrangement.
    if( StringContains(vendor, "GE") || (StringContains(vendor, "SIEMENS") && !SliceMosaic) )
      {
      if( StringContains(vendor, "GE") )
        {
        MeasurementFrame = LPSDirCos;
        }
      else // SIEMENS data assumes a measurement frame that is the identity matrix.
        {
        MeasurementFrame.SetIdentity();
        }
      // has the measurement frame represented as an identity matrix.
      double image0Origin[3];
      allHeaders[0]->GetElementDS(0x0020, 0x0032, 3, image0Origin);
      std::cout << "Slice 0: " << image0Origin[0] << " " << image0Origin[1] << " " << image0Origin[2] << std::endl;

      // assume volume interleaving, i.e. the second dicom file stores
      // the second slice in the same volume as the first dicom file
      double image1Origin[3];
      allHeaders[1]->GetElementDS(0x0020, 0x0032, 3, image1Origin);
      std::cout << "Slice 0: " << image1Origin[0] << " " << image1Origin[1] << " " << image1Origin[2] << std::endl;

      image1Origin[0] -= image0Origin[0];
      image1Origin[1] -= image0Origin[1];
      image1Origin[2] -= image0Origin[2];
      double x1 = image1Origin[0] * (NRRDSpaceDirection[0][2])
        + image1Origin[1] * (NRRDSpaceDirection[1][2])
        + image1Origin[2] * (NRRDSpaceDirection[2][2]);
      if( x1 < 0 )
        {
        SliceOrderIS = false;
        }
      }
    else if( StringContains(vendor, "SIEMENS") && SliceMosaic )
      {
      MeasurementFrame.SetIdentity(); // The DICOM version of SIEMENS that uses private tags
      // has the measurement frame represented as an identity matrix.
      std::cout << "Siemens SliceMosaic......" << std::endl;

      SliceOrderIS = false;

      // for siemens mosaic image, figure out mosaic slice order from 0029|1010
      // copy information stored in 0029,1010 into a string for parsing
      std::string tag;
      allHeaders[0]->GetElementOB(0x0029, 0x1010, tag);
      // parse SliceNormalVector from 0029,1010 tag
      std::vector<double> valueArray(0);
      int                 nItems = ExtractSiemensDiffusionInformation(tag, "SliceNormalVector", valueArray);
      if( nItems != 3 )  // did not find enough information
        {
        std::cout << "Warning: Cannot find complete information on SliceNormalVector in 0029|1010" << std::endl;
        std::cout << "         Slice order may be wrong." << std::endl;
        }
      else if( valueArray[2] > 0 )
        {
        SliceOrderIS = true;
        }

      // parse NumberOfImagesInMosaic from 0029,1010 tag
      valueArray.resize(0);
      nItems = ExtractSiemensDiffusionInformation(tag, "NumberOfImagesInMosaic", valueArray);
      if( nItems == 0 )  // did not find enough information
        {
        std::cout << "Warning: Cannot find complete information on NumberOfImagesInMosaic in 0029|1010" << std::endl;
        std::cout << "         Resulting image may contain empty slices." << std::endl;
        }
      else
        {
        nSliceInVolume = static_cast<int>(valueArray[0]);
        mMosaic = static_cast<int>(ceil(sqrt(valueArray[0]) ) );
        nMosaic = mMosaic;
        }
      std::cout << "Mosaic in " << mMosaic << " X " << nMosaic
                << " blocks (total number of blocks = " << valueArray[0] << ")." << std::endl;
      }
    else if( !multiSliceVolume &&  StringContains(vendor, "PHILIPS") && nSlice > 1 )
    // so this is not a philips multi-frame single dicom file
      {
      MeasurementFrame = LPSDirCos; // Philips oblique scans list the gradients with respect to the
                                    // ImagePatientOrientation.
      SliceOrderIS = true;

      nSliceInVolume = numberOfSlicesPerVolume;
      nVolume = nSlice / nSliceInVolume;

      double image0Origin[3];
      allHeaders[0]->GetElementDS(0x0020, 0x0032, 3, image0Origin);
      std::cout << "Slice 0: " << image0Origin[0] << " " << image0Origin[1] << " " << image0Origin[2] << std::endl;

      // assume volume interleaving, i.e. the second dicom file stores
      // the second slice in the same volume as the first dicom file
      double image1Origin[3];
      allHeaders[nVolume]->GetElementDS(0x0020, 0x0032, 3, image1Origin);
      std::cout << "Slice " << nVolume << ": " << image1Origin[0] << " "
                << image1Origin[1] << " " << image1Origin[2] << std::endl;

      image1Origin[0] -= image0Origin[0];
      image1Origin[1] -= image0Origin[1];
      image1Origin[2] -= image0Origin[2];
      double x0 = image1Origin[0] * (NRRDSpaceDirection[0][2])
        + image1Origin[1] * (NRRDSpaceDirection[1][2])
        + image1Origin[2] * (NRRDSpaceDirection[2][2]);

      // VAM - This needs more investigation -
      // Should we default to false and change based on slice order
      if( x0 < 0 )
        {
        SliceOrderIS = false;
        }
      }
    else if( StringContains(vendor, "PHILIPS") && nSlice == 1 )
      {
      // special handling for philips multi-frame dicom later.
      }
    else
      {
      std::cout << " Warning: vendor type not valid" << std::endl;
      // treate the dicom series as an ordinary image and write a straight nrrd file.
      WriteVolume<VolumeType>( readerOutput, outputVolumeHeaderName );
      FreeHeaders(allHeaders);
      return EXIT_SUCCESS;
      }

    if( SliceOrderIS )
      {
      std::cout << "Slice order is IS" << std::endl;
      }
    else
      {
      std::cout << "Slice order is SI" << std::endl;
      NRRDSpaceDirection[0][2] = -NRRDSpaceDirection[0][2];
      NRRDSpaceDirection[1][2] = -NRRDSpaceDirection[1][2];
      NRRDSpaceDirection[2][2] = -NRRDSpaceDirection[2][2];
      }

    std::cout << "Row: " << (NRRDSpaceDirection[0][0])  << ", "
              << (NRRDSpaceDirection[1][0]) << ", "
              << (NRRDSpaceDirection[2][0]) << std::endl;
    std::cout << "Col: " << (NRRDSpaceDirection[0][1])
              << ", " << (NRRDSpaceDirection[1][1])
              << ", " << (NRRDSpaceDirection[2][1]) << std::endl;
    std::cout << "Sli: " << (NRRDSpaceDirection[0][2])
              << ", " << (NRRDSpaceDirection[1][2]) << ", "
              << (NRRDSpaceDirection[2][2]) << std::endl;

    const float orthoSliceSpacing = fabs( (NRRDSpaceDirection[2][2]) );

    int              nIgnoreVolume = 0; // Used for Philips Trace like images
    std::vector<int> useVolume;

    std::vector<float> bValues(0);
    float              maxBvalue = 0;
    int                nBaseline = 0;

    // UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem is only of debug purposes.
    std::vector<vnl_vector_fixed<double, 3> > DiffusionVectors;
    std::vector<vnl_vector_fixed<double, 3> > UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem;
    std::vector<int>                          ignorePhilipsSliceMultiFrame;

    ////////////////////////////////////////////////////////////
    // vendor dependent tags.
    // read in gradient vectors and determin nBaseline and nMeasurement

    if( StringContains(vendor, "GE") )
      {
#if 0
      // don't even try to convert DTI 6 Direction files
      std::string seriesDescription;
      if( allHeaders[0]->GetElementLO(0x0008, 0x103e, seriesDescription, false) == EXIT_SUCCESS &&
          StringContains(seriesDescription, "6 Directions") )
        {
        std::cerr << "Can't recover B-value & diffusion directions from DTI - 6 Directions scans" << std::endl;
        FreeHeaders(allHeaders);
        return EXIT_FAILURE;
        }
#endif
      std::string ModelName;
      // OK, so there is an accomdation made on the basis of one site
      // having garbage BVal/GVectors.  It has to do with variations
      // of behavior of the Signa HDxt scanner.
      // In all cases, the data is thus:
      // BVal = [0043,1039]
      // GVec[0] = [0019.10bb] GVec[1] = [0019,10bc] GVec[2] = [0019,10bd]
      // there are 3 possible encodings of this data for GE scanners:
      // 1. As IS/DS -- integer an decimal strings -- the normal
      // behavior
      // 2. As OB, but the byte data is binary and may need byte
      // swapping.
      // 3. As OB, but it's actuall as case 1 -- numbers represented
      // as strings.
      // I'm accounting for these cases by looking specifically for
      // the Signa HDxt scanner, and if it doesn't find IS/DS data,
      // look for char strings in the OB data.
      // Honestly this is not an optimal way to handle this
      // situation. In an ideal world we'd have accurate knowledge of
      // what each Scanner/Software Version is doing in these tags,
      // and handle them accordingly. But we don't live in that world.
      bool isSignaHDxt(false);
      if( allHeaders[0]->GetElementLO(0x0008, 0x001090, ModelName, false) == EXIT_SUCCESS &&
          ModelName == "Signa HDxt" )
        {
        isSignaHDxt = true;
        }
      nSliceInVolume = numberOfSlicesPerVolume;
      nVolume = nSlice / nSliceInVolume;

      // assume volume interleaving
      std::cout << "Number of Slices: " << nSlice << std::endl;
      std::cout << "Number of Volume: " << nVolume << std::endl;
      std::cout << "Number of Slices in each volume: " << nSliceInVolume << std::endl;
      for( unsigned int k = 0; k < nSlice; k += nSliceInVolume )
        {
        // parsing bvalue and gradient directions
        vnl_vector_fixed<double, 3> vect3d;
        vect3d.fill( 0 );
        // for some weird reason this item in the GE dicom
        // header is stored as an IS (Integer String) element.
        ::itk::int32_t intb;
        if( !isSignaHDxt )
          {
          allHeaders[k]->GetElementISorOB(0x0043, 0x1039, intb);
          }
        else
          {
          if( allHeaders[k]->GetElementIS(0x0043, 0x1039, intb, false) != EXIT_SUCCESS )
            {
            std::string val;
            allHeaders[k]->GetElementOB(0x0043, 0x1039, val);
            size_t slashpos = val.find('\\');
            val = val.substr(0, slashpos);
            std::stringstream s(val);
            s >> intb;
            }
          }
        float b = static_cast<float>(intb);
        for( unsigned elementNum = 0x10bb; elementNum <= 0x10bd; ++elementNum )
          {
          int vecI(elementNum - 0x10bb);
          if( !isSignaHDxt )
            {
            allHeaders[k]->GetElementDSorOB(0x0019, elementNum, vect3d[vecI]);
            }
          else
            {
            if( allHeaders[k]->GetElementDS(0x0019, elementNum, 1, &vect3d[vecI], false) != EXIT_SUCCESS )
              {
              std::string val;
              allHeaders[k]->GetElementOB(0x0019, elementNum, val);
              std::stringstream s(val);
              s >> vect3d[vecI];
              }
            }
          }

        vect3d[0] = -vect3d[0];
        vect3d[1] = -vect3d[1];

        bValues.push_back( b );
        if( b == 0 )
          {
          vect3d.fill( 0 );
          UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
          DiffusionVectors.push_back(vect3d);
          }
        else
          {
          UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
          // vect3d.normalize();
          DiffusionVectors.push_back(vect3d);
          }

        std::cout << "B-value: " << b
                  << "; diffusion direction: "
                  << DoubleConvert(vect3d[0])
                  << ", "
                  << DoubleConvert(vect3d[1])
                  << ", "
                  << DoubleConvert(vect3d[2]) << std::endl;
        }
      }
    else if( StringContains(vendor, "PHILIPS") )
      {
      if( nSlice > 1 )
        {
        // assume volume interleaving
        std::cout << "Number of Slices: " << nSlice << std::endl;
        std::cout << "Number of Volumes: " << nVolume << std::endl;
        std::cout << "Number of Slices in each volume: " << nSliceInVolume << std::endl;

        std::string tmpString = "";
        // NOTE:  Philips interleaves the directions, so the all gradient directions can be
        // determined in the first "nVolume" slices which represents the first slice from each
        // of the gradient volumes.
        for( unsigned int k = 0; k < nVolume; ++k )
          {
          std::string DiffusionDirectionality;
          bool        useSupplement49Definitions(false);
          if( allHeaders[k]->GetElementCSorOB(0x0018, 0x9075, DiffusionDirectionality, false) == EXIT_SUCCESS )
            {
            useSupplement49Definitions = true;
            }

          bool   B0FieldFound = false;
          double b = 0.0;
          if( useSupplement49Definitions == true )
            {
            B0FieldFound = allHeaders[k]->GetElementFD(0x0018, 0x9087, b, false) == EXIT_SUCCESS;
            }
          else
            {
            float floatB;
            if( allHeaders[k]->GetElementFLorOB(0x2001, 0x1003, floatB, false) == EXIT_SUCCESS )
              {
              B0FieldFound = true;
              }
            if( B0FieldFound )
              {
              b = static_cast<double>(floatB);
              }
            std::string tag;
            allHeaders[k]->GetElementCSorOB(0x2001, 0x1004, tag, false );
            if( StringContains(tag, "I") && b != 0 )
              {
              DiffusionDirectionality = "ISOTROPIC";
              }
            }

          vnl_vector_fixed<double, 3> vect3d;
          vect3d.fill( 0 );
          if( StringContains(DiffusionDirectionality, "ISOTROPIC") )
            { // Deal with images that are to be ignored
            ++nIgnoreVolume;
            useVolume.push_back(0);
            continue;
            }
          else if( ( !B0FieldFound || b == 0 ) || StringContains(DiffusionDirectionality, "NONE") )
            { // Deal with b0 images
            bValues.push_back(b);
            UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
            DiffusionVectors.push_back(vect3d);
            useVolume.push_back(1);
            continue;
            }
          else if( StringContains(DiffusionDirectionality, "DIRECTIONAL") || ( DiffusionDirectionality == "" ) )
            { // Deal with gradient direction images
            bValues.push_back(b);
            useVolume.push_back(1);
            if( useSupplement49Definitions == true )
              {
              double doubleArray[3];
              // Use alternate method to get value out of a sequence header (Some Phillips Data).
              if( allHeaders[k]->GetElementFD(0x0018, 0x9089, 3, doubleArray, false) != EXIT_SUCCESS )
                {
                // std::cout << "Looking for  0018|9089 in sequence 0018,9076" << std::endl;
                // gdcm::SeqEntry *
                // DiffusionSeqEntry=allHeaders[k]->GetSeqEntry(0x0018,0x9076);
                itk::DCMTKSequence DiffusionSeqEntry;
                allHeaders[k]->GetElementSQ(0x0018, 0x9076, DiffusionSeqEntry);
                // const unsigned int
                // n=DiffusionSeqEntry->GetNumberOfSQItems();
                unsigned int n = DiffusionSeqEntry.card();
                if( n == 0 )
                  {
                  std::cout << "ERROR:  Sequence entry 0018|9076 has no items." << std::endl;
                  FreeHeaders(allHeaders);
                  return EXIT_FAILURE;
                  }
                DiffusionSeqEntry.GetElementFD(0x0018, 0x9089, 3, doubleArray);
                }
              vect3d[0] = doubleArray[0];
              vect3d[1] = doubleArray[1];
              vect3d[2] = doubleArray[2];
              std::cout << "===== gradient orientations:" << k << " "
                        << inputFileNames[k] << " (0018,9089) " << " " << vect3d << std::endl;
              }
            else
              {
              float tmp[3];
              /*const bool b0exist =*/
              allHeaders[k]->GetElementFLorOB( 0x2005, 0x10b0, tmp[0] );
              allHeaders[k]->GetElementFLorOB( 0x2005, 0x10b1, tmp[1] );
              allHeaders[k]->GetElementFLorOB( 0x2005, 0x10b2, tmp[2] );
              vect3d[0] = static_cast<double>(tmp[0]);
              vect3d[1] = static_cast<double>(tmp[1]);
              vect3d[2] = static_cast<double>(tmp[2]);
              }

            UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
            // vect3d.normalize();
            DiffusionVectors.push_back(vect3d);
            }
          else // Have no idea why we'd be here so error out
            {
            std::cout << "ERROR: DiffusionDirectionality was "
                      << DiffusionDirectionality << "  Don't know what to do with that..." << std::endl;
            FreeHeaders(allHeaders);
            return EXIT_FAILURE;
            }

          std::cout << "B-value: " << b
                    << "; diffusion direction: "
                    << DoubleConvert(vect3d[0]) << ", "
                    << DoubleConvert(vect3d[1]) << ", "
                    << DoubleConvert(vect3d[2])
                    << std::endl;
          }
        }
      else
        {
        // multi-frame file, everything is inside
        std::map<std::vector<double>, double> gradientDirectionAndBValue;
        ignorePhilipsSliceMultiFrame.clear();

        sliceLocations.clear();
        bValues.clear();
        DiffusionVectors.clear();
        useVolume.clear();

        itk::DCMTKSequence perFrameFunctionalGroup;
        itk::DCMTKSequence innerSeq;
        double             dwbValue;

        allHeaders[0]->GetElementSQ(0x5200, 0x9230, perFrameFunctionalGroup);
        int nItems = perFrameFunctionalGroup.card();

        // have to determine if volume slices are interleaved
        std::string origins[2];
        for( unsigned long i = 0;
             i < static_cast<unsigned long>(perFrameFunctionalGroup.card() ); ++i )
          {
          itk::DCMTKItem curItem;
          perFrameFunctionalGroup.GetElementItem(i, curItem);

          // index slice locations with string origin
          itk::DCMTKSequence originSeq;
          curItem.GetElementSQ(0x0020, 0x9113, originSeq);
          std::string originString;
          originSeq.GetElementDS(0x0020, 0x0032, originString);
          ++sliceLocations[originString];
          // save origin of first 2 slices to compare and see if the
          // volume is interleaved.
          if( i < 2 )
            {
            origins[i] = originString;
            }

          itk::DCMTKSequence mrDiffusionSeq;
          curItem.GetElementSQ(0x0018, 0x9117, mrDiffusionSeq);

          std::string dirValue;
          mrDiffusionSeq.GetElementCSorOB(0x0018, 0x9075, dirValue);

          if( StringContains(dirValue, "ISO") )
            {
            useVolume.push_back(0);
            ignorePhilipsSliceMultiFrame.push_back( i );
            }
          else if( StringContains(dirValue, "NONE") )
            {
            useVolume.push_back(1);
            std::vector<double> v(3);
            v[0] = 0; v[1] = 0; v[2] = 0;
            unsigned int nOld = gradientDirectionAndBValue.size();
            gradientDirectionAndBValue[v] = 0;
            unsigned int nNew = gradientDirectionAndBValue.size();

            if( nOld != nNew )
              {
              vnl_vector_fixed<double, 3> vect3d;
              vect3d.fill( 0 );
              DiffusionVectors.push_back( vect3d );
              UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
              bValues.push_back( 0 );
              }
            }
          else
            {
            useVolume.push_back(1);
            if( mrDiffusionSeq.GetElementDSorOB(0x0018, 0x9087, dwbValue, false) != EXIT_SUCCESS )
              {
              mrDiffusionSeq.GetElementFD(0x0018, 0x9087, dwbValue);
              }
            itk::DCMTKSequence volSeq;
            mrDiffusionSeq.GetElementSQ(0x0018, 0x9076, volSeq);
            double dwgVal[3];
            if( volSeq.GetElementDSorOB<double>(0x0018, 0x9089, 3, dwgVal, false) != EXIT_SUCCESS )
              {
              volSeq.GetElementFD(0x0018, 0x9089, 3, dwgVal);
              }
            std::vector<double> v(3);
            v[0] = dwgVal[0];
            v[1] = dwgVal[1];
            v[2] = dwgVal[2];
            unsigned int nOld = gradientDirectionAndBValue.size();
            gradientDirectionAndBValue[v] = dwbValue;
            unsigned int nNew = gradientDirectionAndBValue.size();

            if( nOld != nNew )
              {
              vnl_vector_fixed<double, 3> vect3d;
              vect3d[0] = v[0]; vect3d[1] = v[1]; vect3d[2] = v[2];
              UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
              // vect3d.normalize();
              DiffusionVectors.push_back( vect3d );

              bValues.push_back( dwbValue);
              }
            }
          }

        numberOfSlicesPerVolume = sliceLocations.size();

        // de-interleave slices if the origins of the first 2 slices
        // are the same.
        if( origins[0] == origins[1] )
          {
          // interleaved image
          DeInterleaveVolume(readerOutput, numberOfSlicesPerVolume, perFrameFunctionalGroup.card() );
          }

        std::cout << "LPS Matrix: " << std::endl << LPSDirCos << std::endl;
        std::cout << "Volume Origin: " << std::endl << ImageOrigin[0] << ","
                  << ImageOrigin[1] << ","  << ImageOrigin[2] << "," << std::endl;
        std::cout << "Number of slices per volume: " << numberOfSlicesPerVolume << std::endl;
        std::cout << "Slice matrix size: " << nRows << " X " << nCols << std::endl;
        std::cout << "Image resolution: " << xRes << ", " << yRes << ", " << sliceSpacing << std::endl;

        NRRDSpaceDirection = LPSDirCos * OrientationMatrix * SpacingMatrix;

        MeasurementFrame = LPSDirCos;

        nSliceInVolume = sliceLocations.size();
        nVolume = nItems / nSliceInVolume;
        nIgnoreVolume = ignorePhilipsSliceMultiFrame.size() / nSliceInVolume;
        for( unsigned int k2 = 0; k2 < bValues.size(); ++k2 )
          {
          std::cout << k2 << ": direction: "
                    << DoubleConvert(DiffusionVectors[k2][0]) << ", "
                    << DoubleConvert(DiffusionVectors[k2][1]) << ", "
                    << DoubleConvert(DiffusionVectors[k2][2])
                    << ", b-value: " << bValues[k2] << std::endl;
          }
        }
      }
    else if( StringContains(vendor, "SIEMENS") )
      {
      int nStride = 1;

      if( !SliceMosaic )
        {
        std::cout << orthoSliceSpacing << std::endl;
        nSliceInVolume = numberOfSlicesPerVolume;
        nVolume = nSlice / nSliceInVolume;
        std::cout << "Number of Slices: " << nSlice << std::endl;
        std::cout << "Number of Volume: " << nVolume << std::endl;
        std::cout << "Number of Slices in each volume: " << nSliceInVolume << std::endl;
        nStride = nSliceInVolume;
        }
      else
        {
        std::cout << "Data in Siemens Mosaic Format" << std::endl;
        nVolume = nSlice;
        std::cout << "Number of Volume: " << nVolume << std::endl;
        std::cout << "Number of Slices in each volume: " << nSliceInVolume << std::endl;
        nStride = 1;
        }
      // JTM - Determine bvalues from all gradients
      for( unsigned int k = 0; k < nSlice; k += nStride )
        {
        // in Siemens, this entry is a 'CSA Header' which is blob
        // of mixed text & binary data.  Pretty annoying but there you
        // have it.
        std::string diffusionInfoString;;
        allHeaders[k]->GetElementOB( 0x0029, 0x1010, diffusionInfoString );

        // parse B_value from 0029,1010 tag
        std::vector<double> valueArray(0);

        int nItems = ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", valueArray);
        if( nItems != 1 )
          {
          vnl_vector_fixed<double, 3> vect3d;
          // B_Value is missing -- the punt position is to count this
          // volume as having a B_value & Gradient Direction of zero
          std::cout << "Warning: Cannot find complete information on B_value in 0029|1010" << std::endl;
          bValues.push_back( 0.0 );
          vect3d.fill( 0.0 );
          UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
          DiffusionVectors.push_back(vect3d);
          continue;
          }

        // we got a 'valid' B-value
        // If we're trusting the gradient directions in the header,
        // then all we need to do here is save the bValue.
        if( !useBMatrixGradientDirections )
          {
          valueArray.resize(0);
          ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", valueArray);

          bValues.push_back( valueArray[0] );
          }
        else
          {
          // JTM - Patch from UNC: fill the nhdr header with the gradient directions and
          // bvalues computed out of the BMatrix
          valueArray.resize(0);
          nItems = ExtractSiemensDiffusionInformation(diffusionInfoString, "B_matrix", valueArray);
          vnl_matrix_fixed<double, 3, 3> bMatrix;

          if( nItems == 6 )
            {
            std::cout << "=============================================" << std::endl;
            std::cout << "BMatrix calculations..." << std::endl;
            // UNC comments: We get the value of the b-value tag in the header.
            // We won't use it as is, but just to locate the B0 images.
            // This check must be added, otherwise the bmatrix of the B0 is not
            // read properly (it's not an actual field in the DICOM header of the B0).
            std::vector<double> bval_tmp(0);
            bool                b0_image = false;

            // UNC comments: Get the bvalue
            nItems = ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", bval_tmp);
            if( bval_tmp[0] == 0 )
              {
              b0_image = true;
              }

            // UNC comments: The principal eigenvector of the bmatrix is to be extracted as
            // it's the gradient direction and trace of the matrix is the b-value

            // UNC comments: Fill out the 3x3 bmatrix with the 6 components read from the
            // DICOM header.
            bMatrix[0][0] = valueArray[0];
            bMatrix[0][1] = valueArray[1];
            bMatrix[0][2] = valueArray[2];
            bMatrix[1][1] = valueArray[3];
            bMatrix[1][2] = valueArray[4];
            bMatrix[2][2] = valueArray[5];
            bMatrix[1][0] = bMatrix[0][1];
            bMatrix[2][0] = bMatrix[0][2];
            bMatrix[2][1] = bMatrix[1][2];

            // UNC comments: Computing the decomposition
            vnl_svd<double> svd(bMatrix);

            // UNC comments: Extracting the principal eigenvector i.e. the gradient direction
            vnl_vector_fixed<double, 3> vect3d;
            vect3d[0] = svd.U(0, 0);
            vect3d[1] = svd.U(1, 0);
            vect3d[2] = svd.U(2, 0);

            std::cout << "BMatrix: " << std::endl;
            std::cout << bMatrix[0][0] << std::endl;
            std::cout << bMatrix[0][1] << "\t" << bMatrix[1][1] << std::endl;
            std::cout << bMatrix[0][2] << "\t" << bMatrix[1][2] << "\t" << bMatrix[2][2] << std::endl;

            // UNC comments: The b-value si the trace of the bmatrix
            const double bvalue = bMatrix[0][0] + bMatrix[1][1] + bMatrix[2][2];
            std::cout << bvalue << std::endl;
            // UNC comments: Even if the bmatrix is null, the svd decomposition set the 1st eigenvector
            // to (1,0,0). So we force the gradient direction to 0 if the bvalue is null
            if( (b0_image == true) || (bvalue == 0) )
              {
              std::cout << "B0 image detected: gradient direction and bvalue forced to 0" << std::endl;
              vect3d[0] = 0;
              vect3d[1] = 0;
              vect3d[2] = 0;
              std::cout << "Gradient coordinates: " << vect3d[0] << " " << vect3d[1] << " " << vect3d[2] << std::endl;
              bValues.push_back(0);
              }
            else
              {
              std::cout << "Gradient coordinates: " << vect3d[0] << " " << vect3d[1] << " " << vect3d[2] << std::endl;
              bValues.push_back(bvalue);
              }
            DiffusionVectors.push_back(vect3d);
            }
          else
            {
            valueArray.resize(0);
            ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", valueArray);
            vnl_vector_fixed<double, 3> vect3d;
            bValues.push_back( valueArray[0] );
            vect3d[0] = 0;
            vect3d[1] = 0;
            vect3d[2] = 0;
            DiffusionVectors.push_back(vect3d);
            }
          }
        }

      if( useBMatrixGradientDirections == false )
        {
        for( unsigned int k = 0; k < nSlice; k += nStride )
          {
          std::cout << "=======================================" << std::endl << std::endl;
          std::string diffusionInfoString;
          allHeaders[k]->GetElementOB(0x0029, 0x1010, diffusionInfoString );

          std::vector<double>         valueArray;
          vnl_vector_fixed<double, 3> vect3d;

          // parse DiffusionGradientDirection from 0029,1010 tag
          valueArray.resize(0);
          int nItems =
            ExtractSiemensDiffusionInformation(diffusionInfoString, "DiffusionGradientDirection", valueArray);
          if( nItems != 3 )  // did not find enough information
            {
            std::cout << "Warning: Cannot find complete information on DiffusionGradientDirection in 0029|1010"
                      << std::endl;
            vect3d.fill( 0 );
            UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
            DiffusionVectors.push_back(vect3d);
            }
          else
            {
            std::cout << "Number of Directions : " << nItems << std::endl;
            std::cout << "   Directions 0: " << valueArray[0] << std::endl;
            std::cout << "   Directions 1: " << valueArray[1] << std::endl;
            std::cout << "   Directions 2: " << valueArray[2] << std::endl;
            double DiffusionVector_magnitude;
            vect3d[0] = valueArray[0];
            vect3d[1] = valueArray[1];
            vect3d[2] = valueArray[2];

            DiffusionVector_magnitude = sqrt(
                (vect3d[0] * vect3d[0]) + (vect3d[1] * vect3d[1]) + (vect3d[2] * vect3d[2]) );

            std::cout << "DiffusionVector_magnitude " << DiffusionVector_magnitude << std::endl;
            if( DiffusionVector_magnitude <= smallGradientThreshold )
              {
              std::cout << "ERROR: Gradient vector with unreasonably small magnitude exists." << std::endl;
              std::cout << "Gradient #" << k << " with magnitude " << DiffusionVector_magnitude << std::endl;
              std::cout << "Please set useBMatrixGradientDirections to calculate gradient directions "
                        << "from the scanner B Matrix to alleviate this problem." << std::endl;
              FreeHeaders(allHeaders);
              return EXIT_FAILURE;
              }

            UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.push_back(vect3d);
            // vect3d.normalize();
            DiffusionVectors.push_back(vect3d);
            int p = bValues.size();
            std::cout << "Image#: " << k
                      << " BV: " << bValues[p - 1] << " GD: "
                      << DoubleConvert(DiffusionVectors[k / nStride][0]) << ","
                      << DoubleConvert(DiffusionVectors[k / nStride][1]) << ","
                      << DoubleConvert(DiffusionVectors[k / nStride][2])
                      << std::endl;
            }
          }
        }
      }
    else
      {
      std::cout << "ERROR: Unknown scanner vendor " << vendor << std::endl;
      std::cout << "       this dti file format is properly handled." << std::endl;
      FreeHeaders(allHeaders);
      return EXIT_FAILURE;
      }

    // dmImage is a constructed volume.
    VolumeType::Pointer dmImage;

    // imgWriter is used to write out image in case it is not a dicom DWI image
    itk::ImageFileWriter<VolumeType>::Pointer imgWriter = itk::ImageFileWriter<VolumeType>::New();

    ///////////////////////////////////////////////
    // Update the number of volumes based on the
    // number to ignore from the header information
    const unsigned int nUsableVolumes = nVolume - nIgnoreVolume;
    std::cout << "Number of usable volumes: " << nUsableVolumes << std::endl;

    if( StringContains(vendor, "GE") ||
        (StringContains(vendor, "SIEMENS") && !SliceMosaic) )
      {
      dmImage = readerOutput;
      }
    else if( StringContains(vendor, "SIEMENS") && SliceMosaic )
      {
      // de-mosaic
      nRows /= mMosaic;
      nCols /= nMosaic;

      // center the volume since the image position patient given in the
      // dicom header was useless
      ImageOrigin[0] = -(nRows * (NRRDSpaceDirection[0][0])
                         + nCols * (NRRDSpaceDirection[0][1])
                         + nSliceInVolume * (NRRDSpaceDirection[0][2]) ) / 2.0;
      ImageOrigin[1] = -(nRows * (NRRDSpaceDirection[1][0])
                         + nCols * (NRRDSpaceDirection[1][1])
                         + nSliceInVolume * (NRRDSpaceDirection[1][2]) ) / 2.0;
      ImageOrigin[2] = -(nRows * (NRRDSpaceDirection[2][0])
                         + nCols * (NRRDSpaceDirection[2][1])
                         + nSliceInVolume * (NRRDSpaceDirection[2][2]) ) / 2.0;

      VolumeType::Pointer img = readerOutput;

      VolumeType::RegionType region = img->GetLargestPossibleRegion();
      VolumeType::SizeType   size = region.GetSize();

      VolumeType::SizeType dmSize = size;
      unsigned int         original_slice_number = dmSize[2] * nSliceInVolume;
      dmSize[0] /= mMosaic;
      dmSize[1] /= nMosaic;
      dmSize[2] = nUsableVolumes * nSliceInVolume;

      region.SetSize( dmSize );
      dmImage = VolumeType::New();
      dmImage->CopyInformation( img );
      dmImage->SetRegions( region );
      dmImage->Allocate();

      VolumeType::RegionType dmRegion = dmImage->GetLargestPossibleRegion();
      dmRegion.SetSize(2, 1);
      region.SetSize(0, dmSize[0]);
      region.SetSize(1, dmSize[1]);
      region.SetSize(2, 1);
      //    int rawMosaic = 0;
      //    int colMosaic = 0;
//      bool bad_slice = false;
//      unsigned int bad_slice_counter = 0;
      for( unsigned int k = 0; k < original_slice_number; ++k )
        {
#if 0
        // bad_gradient_indices never used, so this is pointless
        for( unsigned int j = 0; j < bad_gradient_indices.size(); ++j )
          {
          unsigned int start_bad_slice_number = bad_gradient_indices[j] * nSliceInVolume;
          unsigned int end_bad_slice_number = start_bad_slice_number + (nSliceInVolume - 1);

          if( k >= start_bad_slice_number && k <= end_bad_slice_number )
            {
            bad_slice = true;
            ++bad_slice_counter;
            break;
            }
          else
            {
            bad_slice = false;
            }
          }
        if( bad_slice == false )
#endif
          {
          unsigned int new_k = k /* - bad_slice_counter */;

          dmRegion.SetIndex(2, new_k);
          itk::ImageRegionIteratorWithIndex<VolumeType> dmIt( dmImage, dmRegion );

          // figure out the mosaic region for this slice
          int sliceIndex = k;

          // int nBlockPerSlice = mMosaic*nMosaic;
          int slcMosaic = sliceIndex / (nSliceInVolume);
          sliceIndex -= slcMosaic * nSliceInVolume;
          int colMosaic = sliceIndex / mMosaic;
          int rawMosaic = sliceIndex - mMosaic * colMosaic;
          region.SetIndex( 0, rawMosaic * dmSize[0] );
          region.SetIndex( 1, colMosaic * dmSize[1] );
          region.SetIndex( 2, slcMosaic );

          itk::ImageRegionConstIteratorWithIndex<VolumeType> imIt( img, region );
          for( dmIt.GoToBegin(), imIt.GoToBegin(); !dmIt.IsAtEnd(); ++dmIt, ++imIt )
            {
            dmIt.Set( imIt.Get() );
            }
          }
        }
      }
    else if( StringContains(vendor, "PHILIPS") )
      {
#if 0
      // this code, if you use it, scrambles the image data.
      VolumeType::Pointer img = readerOutput;

      VolumeType::RegionType region = img->GetLargestPossibleRegion();
      VolumeType::SizeType   size = region.GetSize();

      VolumeType::SizeType dmSize = size;
      dmSize[2] = nSliceInVolume * (nUsableVolumes);

      region.SetSize( dmSize );
      dmImage = VolumeType::New();
      dmImage->CopyInformation( img );
      dmImage->SetRegions( region );
      dmImage->Allocate();

      VolumeType::RegionType dmRegion = dmImage->GetLargestPossibleRegion();
      dmRegion.SetSize(2, 1);
      region.SetSize(0, dmSize[0]);
      region.SetSize(1, dmSize[1]);
      region.SetSize(2, 1);

      unsigned int count = 0;
      for( unsigned int i = 0; i < nVolume; ++i )
        {
        if( useVolume[i] == 1 )
          {
          for( unsigned int k = 0; k < nSliceInVolume; ++k )
            {
            dmRegion.SetIndex(0, 0);
            dmRegion.SetIndex(1, 0);
            dmRegion.SetIndex(2, count * (nSliceInVolume) + k);
            itk::ImageRegionIteratorWithIndex<VolumeType> dmIt( dmImage, dmRegion );

            // figure out the region for this slice
            const int sliceIndex = k * nVolume + i;
            region.SetIndex( 0, 0 );
            region.SetIndex( 1, 0 );
            region.SetIndex( 2, sliceIndex );

            itk::ImageRegionConstIteratorWithIndex<VolumeType> imIt( img, region );
            for( dmIt.GoToBegin(), imIt.GoToBegin(); !dmIt.IsAtEnd(); ++dmIt, ++imIt )
              {
              dmIt.Set( imIt.Get() );
              }
            }
          ++count;
          }
        }
      // Verify sizes
      if( count != bValues.size() )
        {
        std::cout << "ERROR:  bValues are the wrong size." <<  count << " != " << bValues.size() << std::endl;
        FreeHeaders(allHeaders);
        return EXIT_FAILURE;
        }
      if( count != DiffusionVectors.size() )
        {
        std::cout << "ERROR:  DiffusionVectors are the wrong size." <<  count << " != " << DiffusionVectors.size()
                  << std::endl;
        FreeHeaders(allHeaders);
        return EXIT_FAILURE;
        }
      if( count != UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.size() )
        {
        std::cout << "ERROR:  UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem are the wrong size."
                  <<  count << " != " << UnmodifiedDiffusionVectorsInDicomLPSCoordinateSystem.size() << std::endl;
        FreeHeaders(allHeaders);
        return EXIT_FAILURE;
        }
#endif
      // INSANE VERY BAD NO GOOD HACK! All the code above is terrible
      // and generates a garbage image. I only noticed it because I
      // though I was actually supposed to write it out to the NRRD
      // file, when in fact the image with skipped volumes is built
      // but never used in the original program.
      dmImage = readerOutput;
      }
    else
      {
      std::cout << "Warning:  invalid vendor found." << std::endl;
      WriteVolume<VolumeType>( readerOutput, outputVolumeHeaderName );
      FreeHeaders(allHeaders);
      return EXIT_SUCCESS;
      }

    // FSLOutput requires a NIfT file
    if( conversionMode != "DicomToFSL" )
      {
      if( !nrrdFormat )
        {
        itk::ImageFileWriter<VolumeType>::Pointer   rawWriter = itk::ImageFileWriter<VolumeType>::New();
        itk::RawImageIO<PixelValueType, 3>::Pointer rawIO = itk::RawImageIO<PixelValueType, 3>::New();
        rawWriter->SetImageIO( rawIO );
        rawIO->SetByteOrderToLittleEndian();
        rawWriter->SetFileName( outputVolumeDataName.c_str() );
        rawWriter->SetInput( dmImage );
        try
          {
          rawWriter->Update();
          }
        catch( itk::ExceptionObject & excp )
          {
          std::cerr << "Exception thrown while reading the series" << std::endl;
          std::cerr << excp << std::endl;
          FreeHeaders(allHeaders);
          return EXIT_FAILURE;
          }
        }
      else if( nUsableVolumes == 1 )
        {
        int rval = WriteVolume<VolumeType>(dmImage, outputVolumeHeaderName);
        // A single usable volume indicates the input is not a DWI file
        // and therefore DWIConvert is simply that -- it
        // converts a DICOM volume to whatever format you specify by way
        // of the output filename.
        FreeHeaders(allHeaders);
        return rval;
        }
      }
    else
      {
      // copy the computed reference frame to the image so that NIFTI
      // writes the correct stuff out. NIFTI writer
      // expects RAS?
      itk::Matrix<double, 3, 3> NIfTIDirCos = LPSDirCos;
      for( unsigned i = 0; i < 3; ++i )
        {
        NIfTIDirCos[i][2] *= -1.0;
        }
      dmImage->SetDirection(NIfTIDirCos);
      VolumeType::SpacingType spacing;
      spacing[0] = xRes; spacing[1] = yRes; spacing[2] = sliceSpacing;
      dmImage->SetSpacing(spacing);
      VolumeType::PointType origin;
      origin[0] = ImageOrigin[0]; origin[1] = ImageOrigin[1]; origin[2] = ImageOrigin[2];
      dmImage->SetOrigin(origin);
      if( Write4DVolume(dmImage, nUsableVolumes, outputVolumeHeaderName) != EXIT_SUCCESS )
        {
        FreeHeaders(allHeaders);
        return EXIT_FAILURE;
        }
      }
    const vnl_matrix_fixed<double, 3, 3> InverseMeasurementFrame = MeasurementFrame.GetInverse();
    //  float bValue = 0;
    for( unsigned int k = 0; k < nUsableVolumes; ++k )
      {
      if( bValues[k] > maxBvalue )
        {
        maxBvalue = bValues[k];
        }
      }

    // construct vector of gradients
    std::vector<std::vector<double> > gradientVectors;
    if( gradientVectorFile != "" )
      {
      // override gradients embedded in file with an external file.
      // FORMAT:
      // <num_gradients>
      // x y z
      // x y z
      // etc
      std::ifstream gradientFile(gradientVectorFile.c_str(), std::ifstream::in);
      unsigned int  numGradients;
      gradientFile >> numGradients;
      if( numGradients != nUsableVolumes )
        {
        std::cerr << "number of Gradients doesn't match number of volumes" << std::endl;
        return EXIT_FAILURE;
        }
      for( unsigned int imageCount = 0; !gradientFile.eof(); ++imageCount )
        {
        std::vector<double> vec(3);
        for( unsigned i = 0; !gradientFile.eof() &&  i < 3; ++i )
          {
          gradientFile >> vec[i];
          }
        gradientVectors.push_back(vec);
        }
      }
    else
      {
      unsigned int original_volume_number = nUsableVolumes /* + bad_gradient_indices.size() */;
      for( unsigned int k = 0; k < original_volume_number; ++k )
        {
        // bool print_gradient = true;
        // for (unsigned int j = 0; j < bad_gradient_indices.size(); ++j)
        //   {
        //   if (k == bad_gradient_indices[j])
        //     {
        //     print_gradient = false;
        //     break;
        //     }
        float scaleFactor = 0;
        //   }
        if( maxBvalue > 0 )
          {
          scaleFactor = sqrt( bValues[k] / maxBvalue );
          }
        std::cout << "For Multiple BValues: " << k << " -- " << bValues[k] << " / " << maxBvalue << " = "
                  << scaleFactor << std::endl;
        std::vector<double> vec(3);
        // if (print_gradient == true)
          {
          if( useIdentityMeaseurementFrame )
            {
            vnl_vector_fixed<double,
                             3> RotatedDiffusionVectors = InverseMeasurementFrame * (DiffusionVectors[k - nBaseline]);
            for( unsigned ind = 0; ind < 3; ++ind )
              {
              vec[ind] = RotatedDiffusionVectors[ind] * scaleFactor;
              }
            }
          else
            {
            if( useBMatrixGradientDirections )
              {
              for( unsigned ind = 0; ind < 3; ++ind )
                {
                vec[ind] = DiffusionVectors[k][ind];
                }
              }
            else
              {
              for( unsigned ind = 0; ind < 3; ++ind )
                {
                vec[ind] = DiffusionVectors[k - nBaseline][ind] * scaleFactor;
                }
              }
            }
          gradientVectors.push_back(vec);
          }
        // else
        //   {
        //   std::cout << "Gradient " << k << " was removed and will not be printed in the NRRD header file." <<
        // std::endl;
        //   }
        }
      }

    //////////////////////////////////////////////
    // write header file
    // This part follows a DWI NRRD file in NRRD format 5.
    // There should be a better way using itkNRRDImageIO.
    if( conversionMode != "DicomToFSL" )
      {
      std::ofstream header;
      // std::string headerFileName = outputDir + "/" + outputFileName;

      header.open(outputVolumeHeaderName.c_str(), std::ios::out | std::ios::binary);
      header << "NRRD0005" << std::endl
             << std::setprecision(17) << std::scientific;

      if( !nrrdFormat )
        {
        header << "content: exists(" << itksys::SystemTools::GetFilenameName(outputVolumeDataName) << ",0)"
               << std::endl;
        }
      header << "type: short" << std::endl;
      header << "dimension: 4" << std::endl;

      // need to check
      header << "space: " << nrrdSpaceDefinition << "" << std::endl;
      // in nrrd, size array is the number of pixels in 1st, 2nd, 3rd, ... dimensions
      header << "sizes: " << nCols << " " << nRows << " " << nSliceInVolume << " " << nUsableVolumes << std::endl;
      header << "thicknesses:  NaN  NaN " << DoubleConvert(sliceSpacing) << " NaN" << std::endl;
      // need to check
      header << "space directions: "
             << "("
             << DoubleConvert(NRRDSpaceDirection[0][0]) << ","
             << DoubleConvert(NRRDSpaceDirection[1][0]) << ","
             << DoubleConvert(NRRDSpaceDirection[2][0])
             << ") "
             << "("
             << DoubleConvert(NRRDSpaceDirection[0][1]) << ","
             << DoubleConvert(NRRDSpaceDirection[1][1]) << ","
             << DoubleConvert(NRRDSpaceDirection[2][1]) << ") "
             << "("
             << DoubleConvert(NRRDSpaceDirection[0][2]) << ","
             << DoubleConvert(NRRDSpaceDirection[1][2]) << ","
             << DoubleConvert(NRRDSpaceDirection[2][2])
             << ") none" << std::endl;
      header << "centerings: cell cell cell ???" << std::endl;
      header << "kinds: space space space list" << std::endl;

      header << "endian: little" << std::endl;
      header << "encoding: raw" << std::endl;
      header << "space units: \"mm\" \"mm\" \"mm\"" << std::endl;
      header << "space origin: "
             << "(" << DoubleConvert(ImageOrigin[0])
             << "," << DoubleConvert(ImageOrigin[1])
             << "," << DoubleConvert(ImageOrigin[2]) << ") " << std::endl;
      if( !nrrdFormat )
        {
        header << "data file: " << itksys::SystemTools::GetFilenameName(outputVolumeDataName) << std::endl;
        }

      // For scanners, the measurement frame for the gradient directions is the same as the
      // Excerpt from http://teem.sourceforge.net/nrrd/format.html definition of "measurement frame:"
      // There is also the possibility that a measurement frame
      // should be recorded for an image even though it is storing
      // only scalar values (e.g., a sequence of diffusion-weighted MR
      // images has a measurement frame for the coefficients of
      // the diffusion-sensitizing gradient directions, and
      // the measurement frame field is the logical store
      // this information).
      // It was noticed on oblique Philips DTI scans that the prescribed protocol directions were
      // rotated by the ImageOrientationPatient amount and recorded in the DICOM header.
      // In order to compare two different scans to determine if the same protocol was prosribed,
      // it is necessary to multiply each of the recorded diffusion gradient directions by
      // the inverse of the LPSDirCos.
      if( useIdentityMeaseurementFrame )
        {
        header << "measurement frame: "
               << "(" << 1 << "," << 0 << "," << 0 << ") "
               << "(" << 0 << "," << 1 << "," << 0 << ") "
               << "(" << 0 << "," << 0 << "," << 1 << ")"
               << std::endl;
        }
      else
        {
        header << "measurement frame: "
               << "(" << DoubleConvert(MeasurementFrame[0][0]) << ","
               << DoubleConvert(MeasurementFrame[1][0]) << ","
               << DoubleConvert(MeasurementFrame[2][0]) << ") "
               << "(" << DoubleConvert(MeasurementFrame[0][1]) << ","
               << DoubleConvert(MeasurementFrame[1][1]) << ","
               << DoubleConvert(MeasurementFrame[2][1]) << ") "
               << "(" << DoubleConvert(MeasurementFrame[0][2]) << ","
               << DoubleConvert(MeasurementFrame[1][2]) << ","
               << DoubleConvert(MeasurementFrame[2][2]) << ")"
               << std::endl;
        }

      header << "modality:=DWMRI" << std::endl;
      // this is the norminal BValue, i.e. the largest one.
      header << "DWMRI_b-value:=" << DoubleConvert(maxBvalue) << std::endl;

      //  the following three lines are for older NRRD format, where
      //  baseline images are always in the begining.
      //  header << "DWMRI_gradient_0000:=0  0  0" << std::endl;
      //  header << "DWMRI_NEX_0000:=" << nBaseline << std::endl;
      //  need to check
      if( gradientVectorFile != "" )
        {
        for( unsigned int imageCount = 0; imageCount < nUsableVolumes; ++imageCount )
          {
          header << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << imageCount << ":="
                 << DoubleConvert(gradientVectors[imageCount][0]) << "   "
                 << DoubleConvert(gradientVectors[imageCount][1]) << "   "
                 << DoubleConvert(gradientVectors[imageCount][2])
                 << std::endl;
          }
        }
      else
        {
        unsigned int original_volume_number = nUsableVolumes /* + bad_gradient_indices.size() */;

        // above, where the gradientVectors is filled out, the same
        // logic is used to skip bad gradients.  Only good gradients are
        // put into the gradientVectors.  The gradient rejection logic
        // is duplicated below, so the gradient numbers will be correct
        // in the nrrd header
        unsigned int gradientVecIndex = 0;
        for( unsigned int k = 0; k < original_volume_number; ++k )
          {
          // bool print_gradient = true;

          // for (unsigned int j = 0; j < bad_gradient_indices.size(); ++j)
          //   {
          //   if (k == bad_gradient_indices[j])
          //     {
          //     print_gradient = false;
          //     break;
          //     }
          //   }
          // if (print_gradient == true)
            {
            header << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << k << ":="
                   << DoubleConvert(gradientVectors[gradientVecIndex][0] ) << "   "
                   << DoubleConvert(gradientVectors[gradientVecIndex][1] ) << "   "
                   << DoubleConvert(gradientVectors[gradientVecIndex][2] )
                   << std::endl;
            ++gradientVecIndex;
            }
          }
        }
      // write data in the same file is .nrrd was chosen
      header << std::endl;;
      if( nrrdFormat )
        {
        unsigned long nVoxels = dmImage->GetBufferedRegion().GetNumberOfPixels();
        header.write( reinterpret_cast<char *>(dmImage->GetBufferPointer() ),
                      nVoxels * sizeof(short) );
        }
      header.close();
      }
    else
      {
      // write out in FSL format
      if( WriteBValues<float>(bValues, outputFSLBValFilename) != EXIT_SUCCESS )
        {
        std::cerr << "Failed to write " << outputFSLBValFilename
                  << std::endl;
        return EXIT_FAILURE;
        }
      if( WriteBVectors<double>(gradientVectors, outputFSLBVecFilename) != EXIT_SUCCESS )
        {
        std::cerr << "Failed to write " << outputFSLBVecFilename
                  << std::endl;
        return EXIT_FAILURE;
        }
      }

    if( writeProtocolGradientsFile == true )
      {
      //////////////////////////////////////////////
      // writeProtocolGradientsFile write protocolGradientsFile file
      // This part follows a DWI NRRD file in NRRD format 5.
      // There should be a better way using itkNRRDImageIO.

      std::ofstream protocolGradientsFile;
      // std::string protocolGradientsFileFileName = outputDir + "/" + outputFileName;

      const std::string protocolGradientsFileName = outputVolumeHeaderName + ".txt";
      protocolGradientsFile.open( protocolGradientsFileName.c_str() );
      protocolGradientsFile << "ImageOrientationPatient (0020|0032): "
                            << DoubleConvert(LPSDirCos[0][0]) << "\\"
                            << DoubleConvert(LPSDirCos[1][0])
                            << "\\" << DoubleConvert(LPSDirCos[2][0]) << "\\"
                            << DoubleConvert(LPSDirCos[0][1])
                            << "\\" << DoubleConvert(LPSDirCos[1][1]) << "\\"
                            << DoubleConvert(LPSDirCos[2][1]) << "\\"
                            << std::endl;
      protocolGradientsFile << "==================================" << std::endl;
      protocolGradientsFile << "Direction Cosines: " << std::endl << LPSDirCos << std::endl;
      protocolGradientsFile << "==================================" << std::endl;
      protocolGradientsFile << "MeasurementFrame: " << std::endl << MeasurementFrame << std::endl;
      protocolGradientsFile << "==================================" << std::endl;
      for( unsigned int k = 0; k < nUsableVolumes; ++k )
        {
        float scaleFactor = 0;
        if( maxBvalue > 0 )
          {
          scaleFactor = sqrt( bValues[k] / maxBvalue );
          }
        protocolGradientsFile << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << k << "=["
                              << DoubleConvert(DiffusionVectors[k - nBaseline][0] * scaleFactor) << ";"
                              << DoubleConvert(DiffusionVectors[k - nBaseline][1] * scaleFactor) << ";"
                              << DoubleConvert(DiffusionVectors[k - nBaseline][2] * scaleFactor) << "]" << std::endl;
        }
      protocolGradientsFile << "==================================" << std::endl;
      for( unsigned int k = 0; k < nUsableVolumes; ++k )
        {
        float scaleFactor = 0;
        if( maxBvalue > 0 )
          {
          scaleFactor = sqrt( bValues[k] / maxBvalue );
          }
        const vnl_vector_fixed<double, 3u> ProtocolGradient = InverseMeasurementFrame * DiffusionVectors[k - nBaseline];
        protocolGradientsFile << "Protocol_gradient_" << std::setw(4) << std::setfill('0') << k << "=["
                              << DoubleConvert(ProtocolGradient[0] * scaleFactor) << ";"
                              << DoubleConvert(ProtocolGradient[1] * scaleFactor) << ";"
                              << DoubleConvert(ProtocolGradient[2] * scaleFactor) << "]" << std::endl;
        }
      protocolGradientsFile << "==================================" << std::endl;
      protocolGradientsFile.close();
      }
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception caught " << excp << std::endl;
    FreeHeaders(allHeaders);
    return EXIT_FAILURE;
    }

  FreeHeaders(allHeaders);
  return EXIT_SUCCESS;
}
