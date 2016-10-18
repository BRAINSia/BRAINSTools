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
#include "djdecode.h"
#include "StringContains.h"
#include "DWIConvertUtils.h"
#include "itkNumberToString.h"

/** the real computation goes on in DWIConverter classes, of which
 * there is one for each manufacturer we encounter.
 */
#include "DWIConverterFactory.h"
#include <BRAINSCommonLib.h>

#include "helpers/loglog.h"
#include "helpers/lloguser.h"

/** the DICOM datasets are read as 3D volumes, but they need to be
 *  written as 4D volumes for image types other than NRRD.
 */
int
Write4DVolume( DWIConverter::VolumeType::Pointer & img, int nVolumes, const std::string & fname )
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
#if DEBUG_WRITE4DVOLUME
    {
    itk::ImageFileWriter<DWIConverter::VolumeType>::Pointer writer = itk::ImageFileWriter<DWIConverter::VolumeType>::New();
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

double
ComputeMaxBvalue(const std::vector<double> &bValues)
{
  double maxBvalue(0.0);
  for( unsigned int k = 0; k < bValues.size(); ++k )
    {
    if( bValues[k] > maxBvalue )
      {
      maxBvalue = bValues[k];
      }
    }
  return maxBvalue;
}

DWIMetaDataDictionaryValidator::GradientTableType
computeScaledDiffusionVectors( const DWIMetaDataDictionaryValidator::GradientTableType &UnitNormDiffusionVectors,
                               const std::vector<double> &bValues,
                               const double maxBvalue)
{
  DWIMetaDataDictionaryValidator::GradientTableType BvalueScaledDiffusionVectors;
  for( unsigned int k = 0; k < UnitNormDiffusionVectors.size(); ++k )
    {
    vnl_vector_fixed<double,3> vec(3);
    float scaleFactor = 0;
    if( maxBvalue > 0 )
      {
      scaleFactor = sqrt( bValues[k] / maxBvalue );
      }
    std::cout << "Scale Factor for Multiple BValues: " << k << " -- sqrt( " << bValues[k] << " / " << maxBvalue << " ) = "
    << scaleFactor << std::endl;
    for( unsigned ind = 0; ind < 3; ++ind )
      {
      vec[ind] = UnitNormDiffusionVectors[k][ind] * scaleFactor;
      }
    BvalueScaledDiffusionVectors.push_back(vec);
    }
  return BvalueScaledDiffusionVectors;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const std::string version = commandLine.getVersion();
  BRAINSRegisterAlternateIO();
  dcmtk::log4cplus::helpers::LogLog::getLogLog()->setQuietMode(true);
  // just need one instance to do double to string conversions
  itk::NumberToString<double> DoubleConvert;

  // build a NRRD file out of FSL output, which is two text files
  // for gradients and b values plus a NIfTI file for the gradient volumes.
  if( conversionMode == "FSLToNrrd" )
    {
    return FSLToNrrd(inputVolume, outputVolume,fslNIFTIFile,
                     inputBValues, inputBVectors, transpose, allowLossyConversion);
    }
  // make FSL file set from a NRRD file.
  if( conversionMode == "NrrdToFSL" )
    {
    return NrrdToFSL(inputVolume, outputVolume,
                     outputBValues, outputBVectors, allowLossyConversion);
    }

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

  DWIConverter::VolumeType::Pointer dmImage = converter->GetDiffusionVolume();
  const DWIMetaDataDictionaryValidator::GradientTableType &UnitNormDiffusionVectors = converter->GetDiffusionVectors();
  const std::vector<double> &bValues = converter->GetBValues();
  const double maxBvalue = ComputeMaxBvalue(bValues);
  const DWIMetaDataDictionaryValidator::GradientTableType &BvalueScaledDiffusionVectors =
    computeScaledDiffusionVectors(UnitNormDiffusionVectors, bValues, maxBvalue);

  if( conversionMode != "DicomToFSL" && !fMRIOutput)
    {
    // if we're writing out NRRD, and the split header/data NRRD
    // format is used, write out the image as a raw volume.
    if( !nrrdFormat )
      {
      itk::ImageFileWriter<DWIConverter::VolumeType>::Pointer
        rawWriter = itk::ImageFileWriter<DWIConverter::VolumeType>::New();
      itk::RawImageIO<DWIConverter::PixelValueType, 3>::Pointer rawIO
        = itk::RawImageIO<DWIConverter::PixelValueType, 3>::New();
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
        std::cerr << "Exception thrown while writing the series to"
                  << outputVolumeDataName << " " << excp << std::endl;
        std::cerr << excp << std::endl;
        delete converter;
        return EXIT_FAILURE;
        }
      }
    }
  else
    {
    // FSLOutput requires a NIfT file
    // copy the computed reference frame to the image so that ITK
    // writes the correct stuff out.
    itk::Matrix<double, 3, 3> NIfTIDirCos = converter->GetLPSDirCos();
    for( unsigned i = 0; i < 3; ++i )
      {
      NIfTIDirCos[i][2] *= -1.0;
      }
    dmImage->SetDirection(NIfTIDirCos);
    dmImage->SetSpacing(converter->GetSpacing());

    DWIConverter::VolumeType::PointType origin = converter->GetOrigin();
    dmImage->SetOrigin(origin);
    // write the image */
    if( Write4DVolume(dmImage, converter->GetNVolume(), outputVolumeHeaderName) != EXIT_SUCCESS )
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

  const vnl_matrix_fixed<double, 3, 3> InverseMeasurementFrame =
    converter->GetMeasurementFrame().GetInverse();

  // construct vector of gradients
  DWIMetaDataDictionaryValidator::GradientTableType gradientVectors;
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
    if( numGradients != converter->GetNVolume() )
      {
      std::cerr << "number of Gradients doesn't match number of volumes" << std::endl;
      return EXIT_FAILURE;
      }
    for( unsigned int imageCount = 0; !gradientFile.eof(); ++imageCount )
      {
      DWIMetaDataDictionaryValidator::GradientDirectionType vec;
      for( unsigned i = 0; !gradientFile.eof() &&  i < 3; ++i )
        {
        gradientFile >> vec[i];
        }
      gradientVectors.push_back(vec);
      }
    }
  else
    {
    // grab the diffusion vectors.
    for( unsigned int k = 0; k < BvalueScaledDiffusionVectors.size(); ++k )
      {
      DWIMetaDataDictionaryValidator::GradientDirectionType vec;
      if( useIdentityMeaseurementFrame )
        {
        vnl_vector_fixed<double,3> RotatedScaledDiffusionVectors =
          InverseMeasurementFrame * (BvalueScaledDiffusionVectors[k]);
        for( unsigned ind = 0; ind < 3; ++ind )
          {
          vec[ind] = RotatedScaledDiffusionVectors[ind];
          }
        }
      else
        {
        for( unsigned ind = 0; ind < 3; ++ind )
          {
          vec[ind] = BvalueScaledDiffusionVectors[k][ind];
          }
        }
      gradientVectors.push_back(vec);
      }
    }

  //////////////////////////////////////////////
  // write header file
  // This part follows a DWI NRRD file in NRRD format 5.
  // There should be a better way using itkNRRDImageIO.
  if( conversionMode != "DicomToFSL" )
    {
#if 0 // Can not write ITK images directly for DWI images, meta data is lost.
#if 1

    Write4DVolume(dmImage,BvalueScaledDiffusionVectors.size(),outputVolumeHeaderName);
#else
      DWIMetaDataDictionaryValidator myDict;
      myDict.SetMeasurementFrame(converter->GetMeasurementFrame());
      myDict.SetBValue(maxBvalue);
      myDict.SetGradientTable(gradientVectors);
      dmImage->SetMetaDataDictionary(myDict.GetMetaDataDictionary());

      itk::ImageFileWriter<DWIConverter::VolumeType>::Pointer writer = itk::ImageFileWriter<DWIConverter::VolumeType>::New();
      writer->SetFileName( outputVolumeHeaderName );
      writer->SetInput( dmImage );
      writer->Update();
#endif

#else
    std::ofstream header;
    // std::string headerFileName = outputDir + "/" + outputFileName;

    header.open(outputVolumeHeaderName.c_str(), std::ios::out | std::ios::binary);
    header << "NRRD0005" << std::endl
           << std::setprecision(17) << std::scientific;

    // stamp with DWIConvert branding
    header << "# This file was created by DWIConvert version " << version << std::endl
           << "# https://github.com/BRAINSia/BRAINSTools" << std::endl
           << "# part of the BRAINSTools package." << std::endl
           << "# Command line options:" << std::endl
           << "# --smallGradientThreshold " << smallGradientThreshold << std::endl;
    if(useIdentityMeaseurementFrame)
      {
      header << "# --useIdentityMeasurementFrame" << std::endl;
      }
    if(useBMatrixGradientDirections)
      {
      header << "# --useBMatrixGradientDirections" << std::endl;
      }
    header << "#" << std::endl << "#" << std::endl;

    if( !nrrdFormat )
      {
      header << "content: exists(" << itksys::SystemTools::GetFilenameName(outputVolumeDataName) << ",0)"
             << std::endl;
      }
    header << "type: short" << std::endl;
    header << "dimension: 4" << std::endl;
    header << "space: " << converter->GetNRRDSpaceDefinition() << "" << std::endl;

    DWIConverter::RotationMatrixType NRRDSpaceDirection =
      converter->GetNRRDSpaceDirection();
    header << "sizes: " << converter->GetCols()
           << " " << converter->GetRows()
           << " " << converter->GetSlicesPerVolume()
           << " " << converter->GetNVolume() << std::endl;
    header << "thicknesses:  NaN  NaN " << DoubleConvert(converter->GetSpacing()[2]) << " NaN" << std::endl;
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

    DWIConverter::VolumeType::PointType
      ImageOrigin = converter->GetOrigin();
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
      DWIConverter::RotationMatrixType MeasurementFrame =
        converter->GetMeasurementFrame();
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
      for( unsigned int imageCount = 0; imageCount < converter->GetNVolume(); ++imageCount )
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
      unsigned int gradientVecIndex = 0;
      for( unsigned int k = 0; k < gradientVectors.size(); ++k )
        {
        header << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << k << ":="
               << DoubleConvert(gradientVectors[gradientVecIndex][0] ) << "   "
               << DoubleConvert(gradientVectors[gradientVecIndex][1] ) << "   "
               << DoubleConvert(gradientVectors[gradientVecIndex][2] )
               << std::endl;
        ++gradientVecIndex;
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
  #endif
    }
  else
    {
    // write out in FSL format
    if( WriteBValues<double>(bValues, outputFSLBValFilename) != EXIT_SUCCESS )
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
    }
  return EXIT_SUCCESS;
}
