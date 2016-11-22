/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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
#include "StringContains.h"
#include "DWIConvertUtils.h"
#include "itkNumberToString.h"

/** the real computation goes on in DWIConverter classes, of which
 * there is one for each manufacturer we encounter.
 */
#include "DWIConverterFactory.h"
#include <BRAINSCommonLib.h>

#include "dcmtk/oflog/helpers/loglog.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmjpeg/djdecode.h"
#include "dcmtk/dcmjpls/djdecode.h"
#include "dcmtk/dcmdata/dcrledrg.h"

std::string
MakeFileComment(bool useBMatrixGradientDirections, bool useIdentityMeaseurementFrame, double smallGradientThreshold,
  const std::string& version)
{
  std::__1::stringstream commentSection;
  {

  commentSection << "#" << std::__1::endl << "#" << std::__1::endl;
  commentSection << "# This file was created by DWIConvert version " << version << std::__1::endl
         << "# https://github.com/BRAINSia/BRAINSTools" << std::__1::endl
         << "# part of the BRAINSTools package." << std::__1::endl
         << "# Command line options:" << std::__1::endl
         << "# --smallGradientThreshold " << smallGradientThreshold << std::__1::endl;
  if (useIdentityMeaseurementFrame) {
    commentSection << "# --useIdentityMeasurementFrame" << std::__1::endl;
  }
  if (useBMatrixGradientDirections) {
    commentSection << "# --useBMatrixGradientDirections" << std::__1::endl;
  }
  }
  return commentSection.str();
}

static void ManualWriteNRRDFile(const std::string& gradientVectorFile, bool useIdentityMeaseurementFrame,
  bool nrrdSingleFileFormat, const std::string& outputVolumeHeaderName, const std::string& outputVolumeDataName,
  DWIConverter* converter, const DWIMetaDataDictionaryValidator::GradientTableType& gradientVectors,
  const std::string commentSection)
{
  itk::NumberToString<double> DoubleConvert;
  std::__1::ofstream header;
  // std::string headerFileName = outputDir + "/" + outputFileName;

  const double maxBvalue = converter->GetMaxBValue();
  header.open(outputVolumeHeaderName.c_str(), std::__1::ios_base::out | std::__1::ios_base::binary);
  header << "NRRD0005" << std::__1::endl
         << std::__1::setprecision(17) << std::__1::scientific;


  header << commentSection;


  // stamp with DWIConvert branding

  if (!nrrdSingleFileFormat) {
    header << "content: exists(" << itksys::SystemTools::GetFilenameName(outputVolumeDataName) << ",0)"
           << std::__1::endl;
  }
  header << "type: short" << std::__1::endl;
  header << "dimension: 4" << std::__1::endl;
  header << "space: " << converter->GetNRRDSpaceDefinition() << "" << std::__1::endl;

  const DWIConverter::RotationMatrixType& NRRDSpaceDirection = converter->GetNRRDSpaceDirection();
  header << "sizes: " << converter->GetCols()
         << " " << converter->GetRows()
         << " " << converter->GetSlicesPerVolume()
         << " " << converter->GetNVolume() << std::__1::endl;
  header << "thicknesses:  NaN  NaN " << DoubleConvert(converter->GetSpacing()[2]) << " NaN" << std::__1::endl;
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
         << ") none" << std::__1::endl;
  header << "centerings: cell cell cell ???" << std::__1::endl;
  header << "kinds: space space space list" << std::__1::endl;

  header << "endian: little" << std::__1::endl;
  header << "encoding: raw" << std::__1::endl;
  header << "space units: \"mm\" \"mm\" \"mm\"" << std::__1::endl;

  const DWIConverter::VolumeType::PointType ImageOrigin = converter->GetOrigin();
  header << "space origin: "
         << "(" << DoubleConvert(ImageOrigin[0])
         << "," << DoubleConvert(ImageOrigin[1])
         << "," << DoubleConvert(ImageOrigin[2]) << ") " << std::__1::endl;
  if (!nrrdSingleFileFormat) {
    header << "data file: " << itksys::SystemTools::GetFilenameName(outputVolumeDataName) << std::__1::endl;
  }

  DWIConverter::RotationMatrixType MeasurementFrame = converter->GetMeasurementFrame();
  if (useIdentityMeaseurementFrame) {
    MeasurementFrame.SetIdentity();
  }
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
           << std::__1::endl;
  }

  header << "modality:=DWMRI" << std::__1::endl;
  // this is the norminal BValue, i.e. the largest one.
  header << "DWMRI_b-value:=" << DoubleConvert(maxBvalue) << std::__1::endl;

  //  the following three lines are for older NRRD format, where
  //  baseline images are always in the begining.
  //  header << "DWMRI_gradient_0000:=0  0  0" << std::endl;
  //  header << "DWMRI_NEX_0000:=" << nBaseline << std::endl;
  //  need to check
  if (gradientVectorFile!="") {
    for (unsigned int imageCount = 0; imageCount<converter->GetNVolume(); ++imageCount) {
      header << "DWMRI_gradient_" << std::__1::setw(4) << std::__1::setfill('0') << imageCount << ":="
             << DoubleConvert(gradientVectors[imageCount][0]) << "   "
             << DoubleConvert(gradientVectors[imageCount][1]) << "   "
             << DoubleConvert(gradientVectors[imageCount][2])
             << std::__1::endl;
    }
  }
  else {
    unsigned int gradientVecIndex = 0;
    for (unsigned int k = 0; k<gradientVectors.size(); ++k) {
      header << "DWMRI_gradient_" << std::__1::setw(4) << std::__1::setfill('0') << k << ":="
             << DoubleConvert(gradientVectors[gradientVecIndex][0]) << "   "
             << DoubleConvert(gradientVectors[gradientVecIndex][1]) << "   "
             << DoubleConvert(gradientVectors[gradientVecIndex][2])
             << std::__1::endl;
      ++gradientVecIndex;
    }
  }
  // write data in the same file is .nrrd was chosen
  header << std::__1::endl;;
  if (nrrdSingleFileFormat) {
    unsigned long nVoxels = converter->GetDiffusionVolume()->GetBufferedRegion().GetNumberOfPixels();
    header.write(reinterpret_cast<char*>(converter->GetDiffusionVolume()->GetBufferPointer()),
      nVoxels*sizeof(short));
  }
  else {
    // if we're writing out NRRD, and the split header/data NRRD
    // format is used, write out the image as a raw volume.
    itk::ImageFileWriter<DWIConverter::VolumeType>::Pointer
      rawWriter = itk::ImageFileWriter<DWIConverter::VolumeType>::New();
    itk::RawImageIO<DWIConverter::PixelValueType,3>::Pointer rawIO
      = itk::RawImageIO<DWIConverter::PixelValueType,3>::New();
    rawWriter->SetImageIO(rawIO);
    rawIO->SetByteOrderToLittleEndian();
    rawWriter->SetFileName(outputVolumeDataName.c_str());
    rawWriter->SetInput(converter->GetDiffusionVolume());
    try {
      rawWriter->Update();
    }
    catch (itk::ExceptionObject& excp) {
      std::__1::cerr << "Exception thrown while writing the series to"
                     << outputVolumeDataName << " " << excp << std::__1::endl;
      std::__1::cerr << excp << std::__1::endl;
      delete converter;
    }
  }
  header.close();
}


/** the DICOM datasets are read as 3D volumes, but they need to be
 *  written as 4D volumes for image types other than NRRD.
 */
int
Write4DVolume( DWIConverter::VolumeType::Pointer img, int nVolumes, const std::string & fname )
{
  typedef itk::Image<DWIConverter::PixelValueType, 4> Volume4DType;

  DWIConverter::VolumeType::SizeType      size3D(img->GetLargestPossibleRegion().GetSize() );
  DWIConverter::VolumeType::DirectionType direction3D(img->GetDirection() );
  DWIConverter::VolumeType::SpacingType   spacing3D(img->GetSpacing() );
  DWIConverter::VolumeType::PointType     origin3D(img->GetOrigin() );

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
  img4D->SetMetaDataDictionary(img->GetMetaDataDictionary());
  size_t bytecount = img4D->GetLargestPossibleRegion().GetNumberOfPixels();
  bytecount *= sizeof(DWIConverter::PixelValueType);
  memcpy(img4D->GetBufferPointer(), img->GetBufferPointer(), bytecount);
//#define DEBUG_WRITE4DVOLUME
#ifdef DEBUG_WRITE4DVOLUME
    {
   {
      //Set the qform and sfrom codes for the MetaDataDictionary.
      itk::MetaDataDictionary & thisDic = img->GetMetaDataDictionary();
      itk::EncapsulateMetaData< std::string >( thisDic, "qform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
      itk::EncapsulateMetaData< std::string >( thisDic, "sform_code_name", "NIFTI_XFORM_UNKNOWN" );
   }
    itk::ImageFileWriter<DWIConverter::VolumeType>::Pointer writer = itk::ImageFileWriter<DWIConverter::VolumeType>::New();
    writer->SetFileName( "/tmp/dwi3dconvert.nii.gz");
    writer->SetInput( img );
    writer->Update();
    }
#endif

   {
      //Set the qform and sfrom codes for the MetaDataDictionary.
      itk::MetaDataDictionary & thisDic = img4D->GetMetaDataDictionary();
      itk::EncapsulateMetaData< std::string >( thisDic, "qform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
      itk::EncapsulateMetaData< std::string >( thisDic, "sform_code_name", "NIFTI_XFORM_UNKNOWN" );
   }
  itk::ImageFileWriter<Volume4DType>::Pointer imgWriter = itk::ImageFileWriter<Volume4DType>::New();

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



int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const std::string version = commandLine.getVersion();
  BRAINSRegisterAlternateIO();
  dcmtk::log4cplus::helpers::LogLog::getLogLog()->setQuietMode(true);

  // register DCMTK codecs, otherwise they will not be available when
  // `itkDCMTKSeriesFileNames` is used to build a list of filenames,
  // so reading series with JPEG transfer syntax will fail.
  DJDecoderRegistration::registerCodecs();
  DcmRLEDecoderRegistration::registerCodecs();

  // just need one instance to do double to string conversions
  itk::NumberToString<double> DoubleConvert;

  if( outputVolume == "" )
  {
    std::cerr << "Missing output volume name" << std::endl;
    return EXIT_FAILURE;
  }

  // build a NRRD file out of FSL output, which is two text files
  // for gradients and b values plus a NIfTI file for the gradient volumes.
  if( conversionMode == "FSLToNrrd" )
    {
    std::vector<std::string> pathElements;
    pathElements.push_back(outputDirectory);
    pathElements.push_back("/");
    pathElements.push_back( outputVolume );
    std::string fullPathOutputVolume = itksys::SystemTools::JoinPath(pathElements);

    return FSLToNrrd(inputVolume, fullPathOutputVolume,fslNIFTIFile,
                     inputBValues, inputBVectors, transpose, allowLossyConversion);
    }
  // make FSL file set from a NRRD file.
  else if( conversionMode == "NrrdToFSL" )
    {
    return NrrdToFSL(inputVolume, outputVolume,
                     outputBValues, outputBVectors, allowLossyConversion);
    }
  //else //DicomToNrrd OR DicomToFSL

  bool nrrdSingleFileFormat(true);
  // check for required parameters
  if( inputDicomDirectory == "" )
    {
    std::cerr << "Missing DICOM input directory path" << std::endl;
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
      nrrdSingleFileFormat = false;
      }
    }
  else
    {
    // FSL output of gradients & BValues
    size_t extensionPos;
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

  // use the factor to instantiate a converter object based on the vender.
  DWIConverterFactory converterFactory(inputDicomDirectory,
                                       useBMatrixGradientDirections,
                                       smallGradientThreshold);
  DWIConverter *converter;
  try
    {
    converter = converterFactory.New();
    }
  catch( itk::ExceptionObject &excp)
    {
    std::cerr << "Exception creating converter " << excp << std::endl;
    return EXIT_FAILURE;
    }


  // read Dicom directory
  try
    {
    converter->LoadDicomDirectory();
    }
  catch( itk::ExceptionObject &excp)
    {
    std::cerr << "Exception creating converter " << excp << std::endl;
    delete converter;
    return EXIT_FAILURE;
    }
  // this is a punt, it will still write out the volume image
  // even if we don't know how to extract gradients.
  if(converterFactory.GetVendor() == "GENERIC")
    {
    std::cerr << "Can't extract DWI data from files created by vendor "
              << converterFactory.GetVendor() << std::endl;
    DWIConverter::VolumeType::Pointer vol = converter->GetDiffusionVolume();
    WriteVolume<DWIConverter::VolumeType>( vol, outputVolumeHeaderName );
    delete converter;
    return EXIT_SUCCESS;
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
    return EXIT_FAILURE;
    }

  if( conversionMode == "DicomToFSL" || fMRIOutput)
    {
    // FSLOutput requires a NIfTI file
    // copy the computed reference frame to the image so that ITK
    // writes the correct stuff out.
    const itk::Matrix<double, 3, 3> NIfTIDirCos = converter->GetLPSDirCos();
        /* // HACK
    for( unsigned i = 0; i < 3; ++i )
      {
      NIfTIDirCos[i][2] *= -1.0;
      }
         */
    converter->GetDiffusionVolume()->SetDirection(NIfTIDirCos);
    converter->GetDiffusionVolume()->SetSpacing(converter->GetSpacing());

    DWIConverter::VolumeType::PointType origin = converter->GetOrigin();
    converter->GetDiffusionVolume()->SetOrigin(origin);
    // write the image */
    if( Write4DVolume(converter->GetDiffusionVolume(), converter->GetNVolume(), outputVolumeHeaderName) != EXIT_SUCCESS )
      {
      delete converter;
      return EXIT_FAILURE;
      }
    else if(fMRIOutput) // skip writing out GVec/BValue files
      {
      delete converter;
      return EXIT_SUCCESS;
      }
    }




  // construct vector of gradients
  DWIMetaDataDictionaryValidator::GradientTableType gradientVectors;
  //TODO:  This should modify the gradient table IMMEDIATELY AFTER initialization.
  //       It should probably be part of the converter construction process!
  if( gradientVectorFile != "" )
    {
      gradientVectors= converter->readOverwriteGradientVectorFile(gradientVectorFile);
    }
  else
    {
      //TODO:  useIdentityMeasurementFrame should be part of the converter construction process as well
      gradientVectors = converter->computeBvalueScaledDiffusionTensors(useIdentityMeaseurementFrame);
    }

  //////////////////////////////////////////////
  // write header file
  // This part follows a DWI NRRD file in NRRD format 5.
  // There should be a better way using itkNRRDImageIO.
  if( conversionMode != "DicomToFSL" )
    {
#if 0 // Can not write ITK images directly for DWI images, meta data is lost.
#if 1

    Write4DVolume(converter->GetDiffusionVolume(),BvalueScaledDiffusionVectors.size(),outputVolumeHeaderName);
#else
      DWIMetaDataDictionaryValidator myDict;
      myDict.SetMeasurementFrame(converter->GetMeasurementFrame());
      myDict.SetBValue(maxBvalue);
      myDict.SetGradientTable(gradientVectors);
      converter->GetDiffusionVolume()->SetMetaDataDictionary(myDict.GetMetaDataDictionary());

      itk::ImageFileWriter<DWIConverter::VolumeType>::Pointer writer = itk::ImageFileWriter<DWIConverter::VolumeType>::New();
      writer->SetFileName( outputVolumeHeaderName );
      writer->SetInput( converter->GetDiffusionVolume() );
      writer->Update();
#endif

#else
    const std::string commentSection = MakeFileComment(useBMatrixGradientDirections, useIdentityMeaseurementFrame,
      smallGradientThreshold, version);

      ManualWriteNRRDFile(gradientVectorFile, useIdentityMeaseurementFrame, nrrdSingleFileFormat,
        outputVolumeHeaderName,
        outputVolumeDataName, converter, gradientVectors, commentSection);
  #endif
    }
  else
    {
    // write out in FSL format
    if( WriteBValues<double>(converter->GetBValues(), outputFSLBValFilename) != EXIT_SUCCESS )
      {
      std::cerr << "Failed to write " << outputFSLBValFilename
                << std::endl;
      return EXIT_FAILURE;
      }
    if( WriteBVectors(gradientVectors, outputFSLBVecFilename) != EXIT_SUCCESS )
      {
      std::cerr << "Failed to write " << outputFSLBVecFilename
                << std::endl;
      return EXIT_FAILURE;
      }
    }

  if( writeProtocolGradientsFile == true )
    {
      std::cerr << "ERROR: DEPRECATED IMPLEMENTED" << std::endl;
      return EXIT_FAILURE;
#if 0
    itk::NumberToString<double> DoubleConvert;
    //////////////////////////////////////////////
    // writeProtocolGradientsFile write protocolGradientsFile file
    // This part follows a DWI NRRD file in NRRD format 5.
    // There should be a better way using itkNRRDImageIO.

    std::ofstream protocolGradientsFile;
    // std::string protocolGradientsFileFileName = outputDir + "/" + outputFileName;

    const std::string protocolGradientsFileName = outputVolumeHeaderName + ".txt";
    DWIConverter::RotationMatrixType LPSDirCos = converter->GetLPSDirCos();
    DWIConverter::RotationMatrixType MeasurementFrame = converter->GetMeasurementFrame();

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
    for( unsigned int k = 0; k < converter->GetNVolume(); ++k )
      {
      protocolGradientsFile << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << k << "=["
                            << DoubleConvert(BvalueScaledDiffusionVectors[k][0]) << ";"
                            << DoubleConvert(BvalueScaledDiffusionVectors[k][1]) << ";"
                            << DoubleConvert(BvalueScaledDiffusionVectors[k][2]) << "]" << std::endl;
      }
    protocolGradientsFile << "==================================" << std::endl;
    for( unsigned int k = 0; k < converter->GetNVolume(); ++k )
      {
      const vnl_vector_fixed<double, 3u> ProtocolGradient = InverseMeasurementFrame * BvalueScaledDiffusionVectors[k];
      protocolGradientsFile << "Protocol_gradient_" << std::setw(4) << std::setfill('0') << k << "=["
                            << DoubleConvert(ProtocolGradient[0]) << ";"
                            << DoubleConvert(ProtocolGradient[1]) << ";"
                            << DoubleConvert(ProtocolGradient[2]) << "]" << std::endl;
      }
    protocolGradientsFile << "==================================" << std::endl;
    protocolGradientsFile.close();
#endif
    }
  return EXIT_SUCCESS;
}
