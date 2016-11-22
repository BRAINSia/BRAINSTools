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


static bool has_valid_nifti_extension( std::string outputVolumeHeaderName )
{
  const size_t NUMEXT=2;
  const char * const extList [NUMEXT] = {".nii.gz", ".nii"};
  for(size_t i = 0 ; i < NUMEXT; ++i)
  {
    const size_t extensionPos = outputVolumeHeaderName.find(extList[i]);
    if( extensionPos != std::string::npos )
    {
      return true;
    }
  }
  {
    std::cerr << "FSL Format output chosen, "
              << "but output Volume not a recognized "
              << "NIfTI filename " << outputVolumeHeaderName
              << std::endl;
    exit(1);
  }
  return false;
}

static DWIConverter * CreateDicomConverter(
  const std::string inputDicomDirectory,
  const bool useBMatrixGradientDirections,
  const double smallGradientThreshold)
{
// check for required parameters
  if( inputDicomDirectory == "" )
  {
    std::cerr << "Missing DICOM input directory path" << std::endl;
    return nullptr;
  }

// use the factor to instantiate a converter object based on the vender.
  DWIConverterFactory converterFactory(inputDicomDirectory,
    useBMatrixGradientDirections,
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


// read Dicom directory
  try
  {
    converter->LoadDicomDirectory();
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
#if 0  //I don't know why this is useful
    DWIConverter::VolumeType::Pointer vol = converter->GetDiffusionVolume();
    WriteVolume<DWIConverter::VolumeType>( vol, outputVolumeHeaderName );
#endif
    delete converter;
    exit(EXIT_SUCCESS);
  }
  return converter;
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

  // decide whether the output is a single file or
  // header/raw pair
  std::string outputVolumeDataName;
  std::string outputFSLBValFilename;
  std::string outputFSLBVecFilename;

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

  const size_t extensionPos = outputVolumeHeaderName.find(".nhdr");
  const bool nrrdSingleFileFormat = ( extensionPos != std::string::npos ) ? false : true;

  DWIConverter *converter;
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
  else if( conversionMode == "DicomToNrrd" || conversionMode == "DicomToFSL")
  {
    if (conversionMode == "DicomToNrrd")
    {
      if( extensionPos != std::string::npos )
      {
        outputVolumeDataName = outputVolumeHeaderName.substr(0, extensionPos);
        outputVolumeDataName += ".raw";
      }
    }
    else if (conversionMode == "DicomToFSL")
    {
      // FSL output of gradients & BValues
      has_valid_nifti_extension(outputVolumeHeaderName);
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
    converter = CreateDicomConverter(inputDicomDirectory,useBMatrixGradientDirections,smallGradientThreshold);
  }
  else
  {
    std::cerr << "Invalid conversion mode" << std::endl;
    exit(-1);
  }

  //Overwrite gradient directions
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
    //TODO:  Refactor Write4DVolume with only outputVolumeHeaderName
    if( converter->Write4DVolume(converter->GetDiffusionVolume(), converter->GetNVolume(), outputVolumeHeaderName) !=
      EXIT_SUCCESS )
    {
      delete converter;
      return EXIT_FAILURE;
    }
    else if(fMRIOutput) // skip writing out GVec/BValue files
    {
      delete converter;
      return EXIT_SUCCESS;
    }
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
  else
  {
  //////////////////////////////////////////////
  // write header file
  // This part follows a DWI NRRD file in NRRD format 5.
  // There should be a better way using itkNRRDImageIO.
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
    //TODO: Move this internal to converter class
    const std::string commentSection = converter->MakeFileComment(useBMatrixGradientDirections,
      useIdentityMeaseurementFrame, smallGradientThreshold, version);

    converter->ManualWriteNRRDFile(gradientVectorFile, useIdentityMeaseurementFrame, nrrdSingleFileFormat,
      outputVolumeHeaderName,
      outputVolumeDataName, gradientVectors, commentSection);
  #endif
    }

  return EXIT_SUCCESS;
}
