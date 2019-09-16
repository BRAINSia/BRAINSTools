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
#include "BRAINSCutApplyModel.h"

#include "itkLabelGeometryImageFilter.h"
#include "itkSimilarityIndexImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "SimilarityIndexCLP.h"

/*
 * This is to analyse the performance of the BRAINSCut result
 * It takes in manual volume and BRAINSCut continuous volume
 * and then apply different threshold with
 * BRAINSCut's post processing method
 */

inline LabelImagePointerType
ThresholdLabelImageToOneValue(LabelImagePointerType inputMaskVolume);

inline LabelImagePointerType
ReadBinaryImageByFilename(std::string filename);

inline WorkingImagePointer
ReadWorkingImageByFilename(std::string filename);

inline float
GetVolume(LabelImagePointerType image);

void
printToScreen(float manualVolume, float annVolume, float SI, float threshold);

void
printHeader();

int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  BRAINSCutApplyModel BRAINSCutPostProcessing;

  if (inputManualVolume == "")
  {
    std::cout << " inputManualVolume is necessary" << std::endl;
    exit(EXIT_FAILURE);
  }
  /* read continuous image */
  LabelImagePointerType manualVolume = ReadBinaryImageByFilename(inputManualVolume);
  manualVolume = ThresholdLabelImageToOneValue(manualVolume);

  /* temporary file to be compared */
  LabelImagePointerType annThresholdVolume;

  /* compute manual volume */
  float floatManualVolume = GetVolume(manualVolume);

  /* set up similarity index computation */
  using SimilarityIndexFilterType = itk::SimilarityIndexImageFilter<LabelImageType, LabelImageType>;
  SimilarityIndexFilterType::Pointer similarityIndexFilter = SimilarityIndexFilterType::New();

  similarityIndexFilter->SetInput1(manualVolume);

  printHeader();
  /** iterate through the threshold */
  for (float threshold = 0.0F; threshold <= 1.00F; threshold += thresholdInterval)
  {
    /* similarity index */
    annThresholdVolume = BRAINSCutPostProcessing.PostProcessingANN(ANNContinuousVolume, threshold);
    similarityIndexFilter->SetInput2(annThresholdVolume);
    similarityIndexFilter->Update();

    printToScreen(
      floatManualVolume, GetVolume(annThresholdVolume), similarityIndexFilter->GetSimilarityIndex(), threshold);
  }

  return 0;
}

inline LabelImagePointerType
ThresholdLabelImageToOneValue(LabelImagePointerType inputMaskVolume)
{
  using ThresholdType = itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>;
  ThresholdType::Pointer thresholder = ThresholdType::New();

  thresholder->SetInput(inputMaskVolume);
  thresholder->SetInsideValue(1);
  thresholder->SetOutsideValue(0);
  thresholder->SetLowerThreshold(1);
  thresholder->Update();

  LabelImagePointerType outputMask = thresholder->GetOutput();
  return outputMask;
}

inline WorkingImagePointer
ReadWorkingImageByFilename(std::string filename)
{
  using WorkingImageReaderType = itk::ImageFileReader<WorkingImageType>;
  WorkingImageReaderType::Pointer reader = WorkingImageReaderType::New();

  reader->SetFileName(filename);
  reader->Update();

  WorkingImagePointer image = reader->GetOutput();
  return image;
}

inline LabelImagePointerType
ReadBinaryImageByFilename(std::string filename)
{
  using BinaryImageReaderType = itk::ImageFileReader<LabelImageType>;
  BinaryImageReaderType::Pointer reader = BinaryImageReaderType::New();

  reader->SetFileName(filename);
  reader->Update();

  LabelImagePointerType image = reader->GetOutput();
  return image;
}

inline float
GetVolume(LabelImagePointerType image)
{
  unsigned char labelValue = 1;

  using MeasureFilterType = itk::LabelStatisticsImageFilter<LabelImageType, LabelImageType>;

  MeasureFilterType::Pointer manualVolumeMeasrueFilter = MeasureFilterType::New();

  manualVolumeMeasrueFilter->SetInput(image);
  manualVolumeMeasrueFilter->SetLabelInput(image);
  manualVolumeMeasrueFilter->Update();

  float count = manualVolumeMeasrueFilter->GetCount(labelValue);

  LabelImageType::SpacingType spacing = image->GetSpacing();

  float volumeOfOneVoxel = 1.0F;
  for (unsigned int i = 0; i < DIMENSION; i++)
  {
    volumeOfOneVoxel *= spacing[i];
  }

  return count * volumeOfOneVoxel;
}

void
printToScreen(float manualVolume, float annVolume, float SI, float threshold)
{
  std::cout << threshold << ", " << manualVolume << ", " << annVolume << ", " << SI << std::endl;
}

void
printHeader()
{
  std::cout << "threshold, manual, ann, SI" << std::endl;
}
