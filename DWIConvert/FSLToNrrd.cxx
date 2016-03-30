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
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "DWIConvertUtils.h"
#include "itksys/SystemTools.hxx"
#include "vnl/vnl_math.h"

#include <itkExtractImageFilter.h>
#include <itkComposeImageFilter.h>
#include "DWIMetaDataDictionaryValidator.h"

typedef short                               PixelValueType;
typedef itk::Image<PixelValueType, 4>       Volume4DType;
typedef itk::Image<PixelValueType, 3>       Volume3DType;
typedef itk::VectorImage<PixelValueType, 3> VectorVolumeType;

int
FSLToNrrd(const std::string & inputVolume,
          const std::string & outputVolume,
          const std::string & fslNIFTIFile,
          const std::string & inputBValues,
          const std::string & inputBVectors,
          bool transpose,
          bool allowLossyConversion
         )
{
  if( (CheckArg<std::string>("Input Volume", inputVolume, "") == EXIT_FAILURE &&
       CheckArg<std::string>("Input Volume", fslNIFTIFile, "") == EXIT_FAILURE) ||
      CheckArg<std::string>("Output Volume", outputVolume, "") == EXIT_FAILURE)
    {
    return EXIT_FAILURE;
    }

  Volume4DType::Pointer inputVol;

  // string to use as template if no bval or bvec filename is given.
  std::string inputVolumeNameTemplate = inputVolume;
  if(fslNIFTIFile.size() > 0)
    {
    if( ReadVolume<Volume4DType>(inputVol, fslNIFTIFile, allowLossyConversion) != EXIT_SUCCESS )
      {
      return EXIT_FAILURE;
      }
    inputVolumeNameTemplate = fslNIFTIFile;
    }
  else if( inputVolume.size() == 0 || ReadVolume<Volume4DType>(inputVol, inputVolume, allowLossyConversion) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }

  std::string _inputBValues = inputBValues;
  if( CheckArg<std::string>("B Values", inputBValues, "") == EXIT_FAILURE )
    {
    _inputBValues = itksys::SystemTools::GetFilenameWithoutExtension(inputVolumeNameTemplate) +
      ".bval";
    }
  std::string _inputBVectors = inputBVectors;
  if( CheckArg<std::string>("B Vectors", inputBVectors, "") == EXIT_FAILURE )
    {
    _inputBVectors = itksys::SystemTools::GetFilenameWithoutExtension(inputVolumeNameTemplate) +
      ".bvec";
    }

  std::vector<double>               BVals;
  std::vector<std::vector<double> > BVecs;
  unsigned int                      bValCount = 0;
  unsigned int                      bVecCount = 0;
  double                            maxBValue(0.0);
  if( ReadBVals(BVals, bValCount, _inputBValues, maxBValue) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  if( ReadBVecs(BVecs, bVecCount, _inputBVectors,transpose) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  if( bValCount != bVecCount )
    {
    std::cerr << "Mismatch between count of B Vectors ("
              << bVecCount << ") and B Values ("
              << bValCount << ")" << std::endl;
    return EXIT_FAILURE;
    }

  // As suggested by Martin Styner:
  // The implicit normalization is implementing by dividing
  // all the gradients (or B-matrices) by the maximal
  // gradient (or B-matrix) magnitude. The magnitude of a
  // gradient direction vector is determined by the usual
  // L^2 norm, and the magnitude of a B-matrix is via the
  // [Frobenius Norm]. It is after this magnitude rescaling
  // that the nominal b-value (given via
  // "DWMRI_b-value:=b") applies... If the nominal b-value
  // is 1000 (via DWMRI_b-value:=1000), then to represent a
  // DWI with b=500, use a gradient vector (or B-matrix) who's norm is sqrt(1/2) the norm for the b=1000 DWI. "
  //
  // So, since all the b-vecs from FSL have norm 1, all the
  // gradients with maximal b-value (which is your Nrrd
  // DWMRI_b-value) should have norm 1 (i.e. for those you
  // can simply copy over the b-vec info). For all other
  // gradients (i.e. those with b-values below the maximal
  // b-value), you need to scale the coordinates of the
  // b-vector by sqrt(this-b-value/max-b-value).
  std::vector<double>::const_iterator bValIt = BVals.begin(),
    bValsEnd = BVals.end();
  std::vector<std::vector<double> >::iterator bVecIt = BVecs.begin(),
    bVecsEnd = BVecs.end();

  for(; bVecIt != bVecsEnd && bValIt != bValsEnd; ++bVecIt, ++bValIt)
    {
    if((*bValIt) == maxBValue)
      {
      continue;
      }
    double scale = std::sqrt((*bValIt) / maxBValue);
    std::vector<double> &cur = *bVecIt;
    for(unsigned int i = 0; i < 3; ++i)
      {
      cur[i] *= scale;
      }
    }

  Volume4DType::SizeType inputSize =
    inputVol->GetLargestPossibleRegion().GetSize();

  Volume4DType::IndexType inputIndex =
    inputVol->GetLargestPossibleRegion().GetIndex();

  const unsigned int volumeCount = inputSize[3];
  if( volumeCount != bValCount )
    {
    std::cerr << "Mismatch between BVector count ("
              << bVecCount << ") and image volume count ("
              << volumeCount << ")" << std::endl;
    return EXIT_SUCCESS;
    }

  // convert from image series to vector voxels
  Volume4DType::SpacingType   inputSpacing = inputVol->GetSpacing();
  std::cout << "Spacing :" << inputSpacing << std::endl;

  ////////
  // "inputVol" is read as a 4D image. Here we convert that to a VectorImageType:
  //
  typedef itk::ExtractImageFilter< Volume4DType, Volume3DType > ExtractFilterType;

  typedef itk::ComposeImageFilter<Volume3DType, VectorVolumeType> ComposeImageFilterType;
  ComposeImageFilterType::Pointer composer= ComposeImageFilterType::New();

  for( size_t componentNumber = 0; componentNumber < inputSize[3]; ++componentNumber )
     {
     Volume4DType::SizeType extractSize = inputSize;
     extractSize[3] = 0;
     Volume4DType::IndexType extractIndex = inputIndex;
     extractIndex[3] = componentNumber;
     Volume4DType::RegionType extractRegion(extractIndex, extractSize);

     ExtractFilterType::Pointer extracter = ExtractFilterType::New();
     extracter->SetExtractionRegion( extractRegion );
     extracter->SetInput( inputVol );
     extracter->SetDirectionCollapseToIdentity();
     extracter->Update();

     composer->SetInput(componentNumber,extracter->GetOutput());
     }
  composer->Update();
  VectorVolumeType::Pointer nrrdVolume = composer->GetOutput();

  const unsigned int nrrdNumOfComponents = nrrdVolume->GetNumberOfComponentsPerPixel();
  std::cout << "Number of components in converted Nrrd volume: " << nrrdNumOfComponents << std::endl;
  if( nrrdNumOfComponents != bVecCount )
    {
    std::cerr << "Mismatch between count of B Vectors ("
    << bVecCount << ") and number of components in converted vector image ("
    << nrrdNumOfComponents << ")" << std::endl;
    return EXIT_FAILURE;
    }
  ////////

  // Define nrrd volume metaData
  DWIMetaDataDictionaryValidator nrrdVolumeValidator;

  /* Fields that need to be set:
   - thickness
   - centerings
   - modality
   - measurement frame
   - b-value
   - gradients

   NOTE: "centerings" and "thickness" should be created based on "volume interleaved".
         If the image is "pixel interleave" (like vectorImage in ITK), the NrrdIO
         will automatically handle the correct permutation.
   */
  // centerings (optional)
  std::vector<std::string> tempCenterings(4,std::string("cell"));
  tempCenterings[3] = "???";
  nrrdVolumeValidator.SetCenterings(tempCenterings);

  // thickness (optional)
  std::vector<double> tempThickness(4,std::numeric_limits<double>::quiet_NaN());
  tempThickness[2] = inputSpacing[2];
  nrrdVolumeValidator.SetThicknesses(tempThickness);

  // modality
  std::string tempModality("DWMRI"); //The only valid DWI modality
  nrrdVolumeValidator.SetModality(tempModality);

  // measurement frame -> it is identity
  std::vector<std::vector<double> > msrFrame(3);
  for( unsigned int saxi = 0; saxi < 3; saxi++ )
    {
    msrFrame[saxi].resize(3);
    for( unsigned int saxj = 0; saxj < 3; saxj++ )
      {
      msrFrame[saxi][saxj] = 0.0;
      }
    }
  msrFrame[0][0] = 1.0; msrFrame[1][1] = 1.0; msrFrame[2][2] = 1.0;
  nrrdVolumeValidator.SetMeasurementFrame(msrFrame);

  // b-value
  nrrdVolumeValidator.SetBValue( maxBValue );

  // Gradient directions
  DWIMetaDataDictionaryValidator::GradientTableType gradientTable( bVecCount );
  MyArrayWrapper<double, 3> BVec_fixedSize;
  for( unsigned int i = 0; i < bVecCount; ++i )
    {
    // convert std::vector to MyArrayWrapper that is a vector with fixed size
    std::copy( BVecs[i].begin(), BVecs[i].begin()+3, BVec_fixedSize.begin() );
    gradientTable[i] = BVec_fixedSize;
    }
  nrrdVolumeValidator.SetGradientTable( gradientTable );

  // Add metaDataDictionary to Nrrd volume
  nrrdVolume->SetMetaDataDictionary(nrrdVolumeValidator.GetMetaDataDictionary());
  // Write Nrrd volume to disk
  typedef itk::ImageFileWriter<VectorVolumeType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( nrrdVolume );
  nrrdWriter->SetFileName( outputVolume );
  nrrdWriter->Update();

  return EXIT_SUCCESS;
}
