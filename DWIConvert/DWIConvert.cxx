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


#include "DWIConverter.h"

/** the real computation goes on in DWIConverter classes, of which
 * there is one for each manufacturer we encounter.
 */
#include "DWIConverterFactory.h"
#include "FSLDWIConverter.h"
#include "NRRDDWIConverter.h"

#include <BRAINSCommonLib.h>

#include "dcmtk/oflog/helpers/loglog.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmjpeg/djdecode.h"
#include "dcmtk/dcmjpls/djdecode.h"
#include "dcmtk/dcmdata/dcrledrg.h"


#undef HAVE_SSTREAM
#include "DWIConvertCLP.h"

static DWIConverter * CreateDicomConverter(
  const std::string inputDicomDirectory,
  const bool useBMatrixGradientDirections,
  const bool transpose,
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

// read Dicom directory
  try
  {
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

  if(fMRIOutput)
  {
    std::cerr << "Deprecated feature no longer supported: --fMRIOutput" << std::endl;
    return EXIT_FAILURE;
  }
  if( gradientVectorFile != "" ) {
    std::cerr << "Deprecated feature no longer supported: --gradientVectorFile" << std::endl;
    return EXIT_FAILURE;
  }

  if( outputVolume == "" )
  {
    std::cerr << "Missing output volume name" << std::endl;
    return EXIT_FAILURE;
  }

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

    DWIConverter::FileNamesContainer filesList;
    filesList.clear();
    filesList.push_back(inputVolume);

    std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
    FSLDWIConverter * FSLconverter = new FSLDWIConverter(filesList, inputBValues, inputBVectors,transpose);
    try
    {
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
  else if( conversionMode == "NrrdToFSL" )
    {
      std::vector<std::string> pathElements;
      pathElements.push_back(outputDirectory);
      pathElements.push_back("/");
      pathElements.push_back( outputVolume );
      std::string fullPathOutputVolume = itksys::SystemTools::JoinPath(pathElements);

      DWIConverter::FileNamesContainer filesList;
      filesList.clear();
      filesList.push_back(inputVolume);

      std::cout << "INPUT VOLUME: " << filesList[0] << std::endl;
      NRRDDWIConverter * NRRDconverter = new NRRDDWIConverter(filesList, transpose);
      useIdentityMeaseurementFrame = true; //Only true is valid for writing FSL
      NRRDconverter->SetUseIdentityMeaseurementFrame(useIdentityMeaseurementFrame);
      try
      {
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
  else if( conversionMode == "DicomToNrrd" || conversionMode == "DicomToFSL")
  {
    converter = CreateDicomConverter(inputDicomDirectory,useBMatrixGradientDirections, transpose,
      smallGradientThreshold);
    if (conversionMode == "DicomToFSL")
    {
      useIdentityMeaseurementFrame = true; //Only true is valid for writing FSL
    }
    converter->SetUseIdentityMeaseurementFrame(useIdentityMeaseurementFrame);
  }
  else
  {
    std::cerr << "Invalid conversion mode" << std::endl;
    exit(-1);
  }

#if 0 //This should use the bvec and bval file formats
  // NEED TO ADD --forceGradientOverwrite, and then read bvec and bval files
  // A test needs to be written for this case
  // ^^^^^^^^^^^^^^^^^^^^^^^ Done Reading Above this line
  //Overwrite gradient directions
  if( gradientVectorFile != "" )
  {
    converter->readOverwriteGradientVectorFile(gradientVectorFile);
  }
#endif
  //^^^^^^^^^^^^^^^^^^^^^^^^^Done modifying above this line vvvvvvvvvvvvvvvvvvvvv Write outputs
  converter->UpdateBVectorOrientation();



  std::string outputVolumeHeaderName(outputVolume);
  { // concatenate with outputDirectory
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
  }

  if( conversionMode == "DicomToFSL"  || conversionMode == "NrrdToFSL" )
  {
    useIdentityMeaseurementFrame = true; //Only true is valid for writing FSL
    converter->SetUseIdentityMeaseurementFrame(useIdentityMeaseurementFrame);
    // write the image */
    converter->WriteFSLFormattedFileSet(outputVolumeHeaderName, outputBValues, outputBVectors);
  }
  else
  {
  //////////////////////////////////////////////
  // write header file
  // This part follows a DWI NRRD file in NRRD format 5.
  // There should be a better way using itkNRRDImageIO.
    const std::string commentSection = converter->MakeFileComment(version,useBMatrixGradientDirections,
      useIdentityMeaseurementFrame, smallGradientThreshold, conversionMode );

    converter->ManualWriteNRRDFile(outputVolumeHeaderName, commentSection);
    std::cout << "Wrote file: " << outputVolumeHeaderName << std::endl;
  }

  return EXIT_SUCCESS;
}
