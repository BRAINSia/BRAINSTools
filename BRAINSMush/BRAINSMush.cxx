/*=========================================================================
Program: BRAINS (Brain Research: Analysis of Images, Networks, and Systems)

Copyright (c) Iowa Mental Health Clinical Research Center. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0.txt

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Input Example:

./BRAINSMush --inputFirstVolume T1.nii.gz --inputSecondVolume T2.nii.gz --inputMaskVolume brainMask.nii.gz[optional]
 --seed "128,128,128" [optional]  --outputVolume
mush_2.nii.gz [optional] --outputMask mask_2.nii.gz [optional] --outputWeightsFile weights.txt [optional]
--boundingBoxSize "90,60,75" [optional] --boundingBoxStart "83,113,80" [optional]

Minimal Input Example:
./BRAINSMush --inputFirstVolume T1.nii.gz --inputSecondVolume T2.nii.gz
=========================================================================*/

#include <string>
#include <iostream>
#include <cstdlib>

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkAmoebaOptimizer.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"

#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkBinaryBallStructuringElement.h"
// #include "itkLevenbergMarquardtOptimizer.h"
#include "itkMixtureStatisticCostFunction.h"

#include "BRAINSMushCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>
#include <itkIO.h>


// a simple print macro for use when debugging
#define PR(x) std::cout << #x " = " << (x) << "\n";


namespace BRAINSMush
{
constexpr int Dimension = 3;
}

namespace
{
using InputPixelType = unsigned char;
using PixelType = float;
using ImageType = itk::Image<PixelType, 3>;
using MaskPixelType = signed short;
using MaskImageType = itk::Image<MaskPixelType, 3>;
// using MaskIndexType = MaskImageType::IndexType;

using MaskImageWriterType = itk::ImageFileWriter<MaskImageType>;

using ReaderType = itk::ImageFileReader<ImageType>;
// using MaskReaderType = itk::ImageFileReader<MaskImageType>;

using ConstIteratorType = itk::ImageRegionConstIterator<ImageType>;
using MaskIteratorType = itk::ImageRegionIterator<MaskImageType>;
using ConstMaskIteratorType = itk::ImageRegionConstIterator<MaskImageType>;

using StructuringElementType = itk::BinaryBallStructuringElement<InputPixelType, BRAINSMush::Dimension>;
} // namespace


ImageType::Pointer
MixtureOptimizer(ImageType::Pointer &     firstImage,
                 ImageType::Pointer &     secondImage,
                 MaskImageType::Pointer & maskImage,
                 std::string              outputWeightsFile)
{
  using MixtureStatisticCostFunctionType = itk::MixtureStatisticCostFunction<ImageType, ImageType>;
  MixtureStatisticCostFunctionType::Pointer twoByTwoCostFunction = MixtureStatisticCostFunctionType::New();
  twoByTwoCostFunction->SetFirstImage(firstImage);
  twoByTwoCostFunction->SetSecondImage(secondImage);
  twoByTwoCostFunction->SetImageMask(maskImage);

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Initialized MixtureStatisticCostFunction! " << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  using OptimizerType = itk::AmoebaOptimizer;
  OptimizerType::Pointer twoByTwoOptimizer = OptimizerType::New();
  twoByTwoOptimizer->SetCostFunction(twoByTwoCostFunction.GetPointer());
  OptimizerType::ParametersType initialParameters(1);
  initialParameters[0] = 1.0;
  twoByTwoOptimizer->SetInitialPosition(initialParameters);
  twoByTwoOptimizer->SetFunctionConvergenceTolerance(.001);
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Updating Optimizer... " << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  try
  {
    twoByTwoOptimizer->StartOptimization();
  }
  catch (const itk::ExceptionObject & exp)
  {
    std::cerr << "LMO FAIL" << std::endl;
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }
  // End of Sending to Optimizer
  /* ------------------------------------------------------------------------------------
   */

  OptimizerType::ParametersType optimalParameters = twoByTwoOptimizer->GetCurrentPosition();
  OptimizerType::MeasureType    optimalMeasures = twoByTwoOptimizer->GetValue();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Obtained Output from Optimizer!" << optimalParameters << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  /* ------------------------------------------------------------------------------------
   * Save the optimization:  one line with the image coeficients in order,
   * and one line with the error measure ito desired mean and variance.  Print
   * it all out, too.
   */
  const double secondImageBlendValue = optimalParameters[0];

  std::cout << "First Weight:  " << secondImageBlendValue << std::endl;
  std::cout << "Optimality of Variance:  " << optimalMeasures << std::endl << std::endl;

  // write a text file named outputWeightsFile

  std::ofstream to(outputWeightsFile.c_str());
  if (to.is_open())
  {
    to << secondImageBlendValue << "  " << std::endl;
    to << optimalMeasures << "  " << optimalMeasures << std::endl;
    to.close();
  }
  else
  {
    std::cout << "Can't open file for writing! --- " << outputWeightsFile << std::endl;
  }

  /* ------------------------------------------------------------------------------------
   * declare and compute mixtureImage.
   */

  ImageType::Pointer mixtureImage = ImageType::New();
  mixtureImage->SetRegions(firstImage->GetLargestPossibleRegion());
  mixtureImage->SetSpacing(firstImage->GetSpacing());
  mixtureImage->SetOrigin(firstImage->GetOrigin());
  mixtureImage->SetDirection(firstImage->GetDirection());
  mixtureImage->Allocate();

  ConstIteratorType firstIt(firstImage, firstImage->GetRequestedRegion());
  ConstIteratorType secondIt(secondImage, secondImage->GetRequestedRegion());

  using MixtureIteratorType = itk::ImageRegionIterator<ImageType>;
  MixtureIteratorType mixtureIt(mixtureImage, mixtureImage->GetRequestedRegion());
  for (mixtureIt.GoToBegin(), firstIt.GoToBegin(), secondIt.GoToBegin(); !mixtureIt.IsAtEnd();
       ++mixtureIt, ++firstIt, ++secondIt)
  {
    const auto & firstValue = firstIt.Get();
    const auto & secondValue = secondIt.Get();
    const double mixtureValue = 0.5 * (firstValue + secondImageBlendValue * secondValue);
    mixtureIt.Set(mixtureValue);
  }
  return mixtureImage;
}

void
GenerateBrainVolume(ImageType::Pointer &     firstImage,
                    ImageType::Pointer &     secondImage,
                    MaskImageType::Pointer & maskImage,
                    std::string              inputMaskVolume,
                    double                   lowerThresholdFactor,
                    double                   upperThresholdFactor,
                    std::vector<int>         boundingBoxSize,
                    std::vector<int>         boundingBoxStart,
                    // std::vector<int> seed,
                    std::string outputVolume,
                    //  std::string outputMask,
                    std::string              outputWeightsFile,
                    MaskImageType::Pointer & resultImage)
{
  /* ------------------------------------------------------------------------------------
   * Send to Optimizer
   */

  ImageType::Pointer mixtureImage = MixtureOptimizer(firstImage, secondImage, maskImage, outputWeightsFile);

  /* ------------------------------------------------------------------------------------
   * Generate brain volume mask (adapted but heavily modified from proc
   * MushPiece in brainsAutoWorkupPhase2.tcl)
   */

  /* ------------------------------------------------------------------------------------
   * Write out MUSH Image
   */
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->UseCompressionOn();
  writer->SetInput(mixtureImage);
  writer->SetFileName(outputVolume);
  try
  {
    std::cout << "Writing mixture image: " << outputVolume << std::endl;
    writer->Update();
  }
  catch (const itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  /* ------------------------------------------------------------------------------------
   * If region of interest mask is supplied, then use it to generate an initial
   * brain mask and perform the thresholding calculations
   * Otherwise, use the initializer region to compute same
   */

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Calculating thresholds..." << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  // using MaskImageWriterType = itk::ImageFileWriter<MaskImageType>;
  MaskImageWriterType::Pointer maskWriter = MaskImageWriterType::New();
  maskWriter->UseCompressionOn();

  double mean;
  double upper;
  double lower;

  if (inputMaskVolume == "no_mask_exists")
  {
    /* ------------------------------------------------------------------------------------
     * Draw a cuboid box of fixed dimensions inside the brain and use it to
     * determine mean and bounds
     */

    int    voxelCount = 0;
    double signalTotal = 0.0;

    MaskImageType::RegionType            regionOfInterest;
    MaskImageType::RegionType::SizeType  regionOfInterestSize;
    MaskImageType::RegionType::IndexType regionOfInterestStart;

    regionOfInterestSize[0] = boundingBoxSize[0];
    regionOfInterestSize[1] = boundingBoxSize[1];
    regionOfInterestSize[2] = boundingBoxSize[2];

    regionOfInterest.SetSize(regionOfInterestSize);

    regionOfInterestStart[0] = boundingBoxStart[0];
    regionOfInterestStart[1] = boundingBoxStart[1];
    regionOfInterestStart[2] = boundingBoxStart[2];

    regionOfInterest.SetIndex(regionOfInterestStart);

    ConstIteratorType     volumeIt(mixtureImage, regionOfInterest);
    ConstMaskIteratorType labelIt(maskImage, regionOfInterest);
    for (volumeIt.GoToBegin(), labelIt.GoToBegin(); !volumeIt.IsAtEnd() && !labelIt.IsAtEnd(); ++volumeIt, ++labelIt)
    {
      MaskPixelType labelValue = labelIt.Get();
      if (labelValue == 1)
      {
        PixelType signalValue = volumeIt.Get();
        signalTotal += signalValue;
      }
      ++voxelCount;
    }
    PR(voxelCount);
    PR(signalTotal);
    mean = signalTotal / voxelCount;
    // these definitions use magic numbers obtained through manual thresholding
    // and experimentation
    lower = (mean / lowerThresholdFactor);
    upper = (mean / upperThresholdFactor);
  }
  else
  {
    /* ------------------------------------------------------------------------------------
     * Binary erosion to generate initial brain mask
     */

    /* ------------------------------------------------------------------------------------
     * Perform binary threshold on image
     */
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Performing Initial Binary Threshold..." << std::endl;
    std::cout << "---------------------------------------------------" << std::endl << std::endl;
    using BinaryThresholdMaskFilterType = itk::BinaryThresholdImageFilter<MaskImageType, MaskImageType>;
    BinaryThresholdMaskFilterType::Pointer threshToBrainCoreMask = BinaryThresholdMaskFilterType::New();
    threshToBrainCoreMask->SetInput(maskImage);
    threshToBrainCoreMask->SetLowerThreshold(1);
    threshToBrainCoreMask->SetUpperThreshold(1);
    threshToBrainCoreMask->SetInsideValue(1);
    threshToBrainCoreMask->SetOutsideValue(0);

    try
    {
      threshToBrainCoreMask->Update();
    }
    catch (const itk::ExceptionObject & excp)
    {
      std::cerr << "Exception caught ! " << std::endl;
      std::cerr << excp << std::endl;
      return;
    }

    using binaryErodeFilterType = itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType>;

    int erosionValue = 7;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Beginning initial erosion..." << std::endl;
    std::cout << "Eroding by: " << erosionValue << std::endl;
    std::cout << "---------------------------------------------------" << std::endl << std::endl;

    binaryErodeFilterType::Pointer initialMaskImage = binaryErodeFilterType::New();

    StructuringElementType structuringElement;
    structuringElement.SetRadius(erosionValue);
    structuringElement.CreateStructuringElement();
    initialMaskImage->SetKernel(structuringElement);
    initialMaskImage->SetInput(threshToBrainCoreMask->GetOutput());

    // Templating requires different writer to output mask images (type short
    // for initialMaskImage vs type float for mixtureImage)

    try
    {
      initialMaskImage->Update();
    }
    catch (const itk::ExceptionObject & exp)
    {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << exp << std::endl;
    }

    /* ------------------------------------------------------------------------------------
     * Obtain mean of image; calculate lower and upper bounds
     */

    using LabelFilterType = itk::LabelStatisticsImageFilter<ImageType, MaskImageType>;
    LabelFilterType::Pointer labelFilter = LabelFilterType::New();
    labelFilter->SetInput(mixtureImage);
    labelFilter->SetLabelInput(initialMaskImage->GetOutput());

    try
    {
      labelFilter->Update();
    }
    catch (const itk::ExceptionObject & excp)
    {
      std::cerr << "Exception caught ! " << std::endl;
      std::cerr << excp << std::endl;
      return;
    }

    // using StatisticRealType = LabelFilterType::RealType;
    mean = labelFilter->GetMean(1);

    // these definitions use magic numbers obtained through manual thresholding
    // and experimentation
    lower = (mean / lowerThresholdFactor);
    upper = (mean / upperThresholdFactor);
  }

  std::cout << "MushROI Mean:   " << mean << std::endl
            << "Lower Bound:    " << lower << std::endl
            << "Upper Bound:    " << upper << std::endl
            << std::endl;

  /* ------------------------------------------------------------------------------------
   * Perform binary threshold on image
   */
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Performing Initial Binary Threshold..." << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;
  using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, MaskImageType>;
  BinaryThresholdFilterType::Pointer threshToHeadMask = BinaryThresholdFilterType::New();
  threshToHeadMask->SetInput(mixtureImage);
  auto smaller = lower < upper ? lower : upper;
  auto larger = lower < upper ? upper : lower;
  threshToHeadMask->SetLowerThreshold(static_cast<float>(smaller));
  threshToHeadMask->SetUpperThreshold(static_cast<float>(larger));
  threshToHeadMask->SetInsideValue(1);
  threshToHeadMask->SetOutsideValue(0);

  try
  {
    threshToHeadMask->Update();
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Exception caught ! " << std::endl;
    std::cerr << excp << std::endl;
    // return EXIT_FAILURE;
  }

  using BinaryDilateFilterType = itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType>;

  using BinaryImageErodeFilterType = itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType>;

  double ClosingSize = 6;

  const ImageType::SpacingType & spacing = mixtureImage->GetSpacing();

  // Compute minimum object size as the number of voxels in the structuring
  // element, an ellipsoidal ball.
  const double FourThirdsPi = 3.141592653589793238459 * 1.333333333333333333333;
  double       ClosingElementVolume = FourThirdsPi * ClosingSize * ClosingSize * ClosingSize;
  double       VoxelVolume = spacing[0] * spacing[1] * spacing[2];
  int          MinimumObjectSize = static_cast<int>(ClosingElementVolume / VoxelVolume);

  // Define binary erosion and dilation structuring element
  StructuringElementType           ball;
  StructuringElementType::SizeType ballSize;
  for (int d = 0; d < 3; d++)
  {
    ballSize[d] = static_cast<int>((0.5 * ClosingSize) / spacing[d]);
  }
  ball.SetRadius(ballSize);
  ball.CreateStructuringElement();

  /* ------------------------------------------------------------------------------------
   * Binary erosion
   */
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Eroding largest filled region..." << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;
  BinaryImageErodeFilterType::Pointer binaryErodeFilter = BinaryImageErodeFilterType::New();
  binaryErodeFilter->SetErodeValue(1);
  binaryErodeFilter->SetKernel(ball);

  /* ------------------------------------------------------------------------------------
   * Update the ITK pipeline
   */
  try
  {
    MaskImageType::Pointer thresholdOutput = threshToHeadMask->GetOutput();
    thresholdOutput->DisconnectPipeline();

    binaryErodeFilter->SetInput(thresholdOutput);
    binaryErodeFilter->Update();
  }
  catch (const itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  /* ------------------------------------------------------------------------------------
   * Perform binary threshold on image
   */
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Performing Special Binary Threshold..." << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;
  using BinaryThresholdMaskFilterType = itk::BinaryThresholdImageFilter<MaskImageType, MaskImageType>;
  BinaryThresholdMaskFilterType::Pointer threshToBrainCoreMask = BinaryThresholdMaskFilterType::New();
  threshToBrainCoreMask->SetInput(binaryErodeFilter->GetOutput());
  threshToBrainCoreMask->SetLowerThreshold(1);
  threshToBrainCoreMask->SetUpperThreshold(1);
  threshToBrainCoreMask->SetInsideValue(1);
  threshToBrainCoreMask->SetOutsideValue(0);

  try
  {
    threshToBrainCoreMask->Update();
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Exception caught ! " << std::endl;
    std::cerr << excp << std::endl;
    return;
  }

  /* ------------------------------------------------------------------------------------
   * Obtain Largest region filled mask
   */
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Obtaining Largest Filled Region..." << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  using CCImageType = itk::Image<unsigned int, 3>;
  using ConnectedComponentFilterType = itk::ConnectedComponentImageFilter<MaskImageType, CCImageType>;
  using RelabelComponentFilterType = itk::RelabelComponentImageFilter<CCImageType, MaskImageType>;
  using ThresholdFilterType = itk::ThresholdImageFilter<MaskImageType>;

  ConnectedComponentFilterType::Pointer filter = ConnectedComponentFilterType::New();
  RelabelComponentFilterType::Pointer   relabel = RelabelComponentFilterType::New();

  ThresholdFilterType::Pointer LargestFilter = ThresholdFilterType::New();
  LargestFilter->SetOutsideValue(0);
  LargestFilter->ThresholdAbove(1);
  LargestFilter->ThresholdBelow(1);

  /* ------------------------------------------------------------------------------------
   * Update the ITK pipeline
   */

  try
  {
    {
      MaskImageType::Pointer threshToBrainCoreOutput = threshToBrainCoreMask->GetOutput();
      threshToBrainCoreOutput->DisconnectPipeline();

      filter->SetInput(threshToBrainCoreOutput); // ConnectedComponentFilter
      filter->Update();
    }

    {
      CCImageType::Pointer filterOutput = filter->GetOutput();
      filterOutput->DisconnectPipeline();

      relabel->SetInput(filterOutput);

      if (MinimumObjectSize > 0)
      {
        relabel->SetMinimumObjectSize(MinimumObjectSize);
        std::cerr << "MinimumObjectSize: " << MinimumObjectSize << std::endl;
      }

      relabel->Update();

      unsigned int numObjects = relabel->GetNumberOfObjects();
      std::cout << "Removed " << numObjects - 1 << " smaller objects." << std::endl << std::endl;
    }

    {
      MaskImageType::Pointer relabelOutput = relabel->GetOutput();
      relabelOutput->DisconnectPipeline();

      LargestFilter->SetInput(relabelOutput);
      LargestFilter->Update();
    }
  }
  catch (const itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  /* ------------------------------------------------------------------------------------
   * Binary dilation
   */
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Dilating largest filled region..." << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;
  BinaryDilateFilterType::Pointer binaryDilateFilter = BinaryDilateFilterType::New();
  binaryDilateFilter->SetDilateValue(1);
  binaryDilateFilter->SetKernel(ball);

  /* ------------------------------------------------------------------------------------
   * Update the ITK pipeline
   */
  try
  {
    MaskImageType::Pointer largestFilterOutput = LargestFilter->GetOutput();
    largestFilterOutput->DisconnectPipeline();

    binaryDilateFilter->SetInput(largestFilterOutput);
    binaryDilateFilter->Update();
  }
  catch (const itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
    return;
  }

  /* ------------------------------------------------------------------------------------
   * Perform binary threshold on image
   */
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Performing Special Binary Threshold..." << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;
  BinaryThresholdMaskFilterType::Pointer threshToClosureMask = BinaryThresholdMaskFilterType::New();
  threshToClosureMask->SetInput(binaryDilateFilter->GetOutput());
  threshToClosureMask->SetLowerThreshold(1);
  threshToClosureMask->SetUpperThreshold(1);
  threshToClosureMask->SetInsideValue(1);
  threshToClosureMask->SetOutsideValue(0);

  try
  {
    threshToClosureMask->Update();
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Exception caught ! " << std::endl;
    std::cerr << excp << std::endl;
    return;
  }

  MaskImageType::Pointer dilatedOutput = threshToClosureMask->GetOutput();
  dilatedOutput->DisconnectPipeline();

  using LFFMaskFilterType = itk::LargestForegroundFilledMaskImageFilter<MaskImageType>;
  LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
  LFF->SetInput(dilatedOutput);
  LFF->SetOtsuPercentileThreshold(0);
  LFF->SetClosingSize(5);
  LFF->Update();
  resultImage = LFF->GetOutput();
}


static ImageType::Pointer
LoadImage(const std::string & imageName)
{
  return itkUtil::ScaleAndCast<ImageType, ImageType>(itk::ReadImage<ImageType>(imageName), 0, 2048);
}

static MaskImageType::Pointer
LoadMaskImage(const std::string & imageName)
{
  //  typename itk::BinaryThresholdImageFilter<MaskImageType,MaskImageType>::Pointer thf =
  //    itk::BinaryThresholdImageFilter<MaskImageType,MaskImageType>::New();
  //  thf->SetInput(itk::ReadImage<MaskImageType>(imageName));
  //  thf->SetLowerThreshold(0);
  //  thf->SetUpperThreshold(itk::NumericTraits<typename MaskImageType::PixelType>::max());
  //  thf->Update();
  //  return thf->GetOutput();
  return itk::ReadImage<MaskImageType>(imageName);
}

static MaskImageType::Pointer
GenerateInitializerRegion(ImageType::Pointer & referenceImage,
                          std::vector<int>     boundingBoxSize,
                          std::vector<int>     boundingBoxStart)
{
  MaskImageType::Pointer initializeMask = MaskImageType::New();

  initializeMask->SetRegions(referenceImage->GetLargestPossibleRegion());
  initializeMask->CopyInformation(referenceImage);
  initializeMask->Allocate();

  MaskImageType::RegionType            regionOfInterest;
  MaskImageType::RegionType::IndexType regionStart;
  MaskImageType::RegionType::SizeType  regionSize;

  regionStart[0] = boundingBoxStart[0];
  regionStart[1] = boundingBoxStart[1];
  regionStart[2] = boundingBoxStart[2];

  regionSize[0] = boundingBoxSize[0];
  regionSize[1] = boundingBoxSize[1];
  regionSize[2] = boundingBoxSize[2];

  regionOfInterest.SetSize(regionSize);
  regionOfInterest.SetIndex(regionStart);

  MaskIteratorType initializeIt(initializeMask, regionOfInterest);
  // sets a binary cuboid mask
  for (initializeIt.GoToBegin(); !initializeIt.IsAtEnd(); ++initializeIt)
  {
    initializeIt.Set(1);
  }

  return initializeMask;
}


int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  bool debug = true;
  if (debug)
  {
    std::cout << "First Mixture Component Image: " << inputFirstVolume << std::endl;
    std::cout << "Second Mixture Component Image: " << inputSecondVolume << std::endl;
    std::cout << "Region Of Interest Image Mask: " << inputMaskVolume << std::endl;

    std::cout << "Seed Point: {" << seed[0] << ", " << seed[1] << ", " << seed[2] << "}" << std::endl;

    std::cout << "Bounding Box Size: {" << boundingBoxSize[0] << ", " << boundingBoxSize[1] << ", "
              << boundingBoxSize[2] << "}" << std::endl;

    std::cout << "Bounding Box Start: {" << boundingBoxStart[0] << ", " << boundingBoxStart[1] << ", "
              << boundingBoxStart[2] << "}" << std::endl;

    std::cout << "Output Image Name: " << outputVolume << std::endl;
    std::cout << "Output Mask Name: " << outputMask << std::endl;
    std::cout << "Output Weights File: " << outputWeightsFile << std::endl;
    std::cout << "Preliminary Lower Threshold Factor: " << lowerThresholdFactorPre << std::endl;
    std::cout << "Preliminary Upper Threshold Factor: " << upperThresholdFactorPre << std::endl << std::endl;
    std::cout << "Main Lower Threshold Factor: " << lowerThresholdFactor << std::endl;
    std::cout << "Main Upper Threshold Factor: " << upperThresholdFactor << std::endl << std::endl;
  }
  /* ------------------------------------------------------------------------------------
   * Load Images
   */
  ImageType::Pointer firstImage = LoadImage(inputFirstVolume);
  ImageType::Pointer secondImage = LoadImage(inputSecondVolume);

  /* ------------------------------------------------------------------------------------
   * Load or automatically generate the MaskImage to define the region on which
   * to optimize mixture image uniformity.
   * In the case where no mask was given on the command line, the method calls
   * for GenerateBrainVolume to be run twice.
   * Twice.
   */
  MaskImageType::Pointer maskImage = MaskImageType::New();

  if (inputMaskVolume == "no_mask_exists") // "no_mask_exists" is the default
                                           // when no mask is specified on
                                           // command line
  {
    MaskImageType::Pointer boxImage = GenerateInitializerRegion(firstImage, boundingBoxSize, boundingBoxStart);
    GenerateBrainVolume(firstImage,
                        secondImage,
                        boxImage,
                        inputMaskVolume,
                        lowerThresholdFactorPre,
                        upperThresholdFactorPre,
                        boundingBoxSize,
                        boundingBoxStart,
                        outputVolume,
                        outputWeightsFile,
                        maskImage);
    inputMaskVolume = "The_mask_was_generated"; // used as an anti-sentinel
  }
  else
  {
    maskImage = LoadMaskImage(inputMaskVolume);
  }

  // Use the MaskImage from above to define the region on which to optimize
  // mixture image uniformity.
  MaskImageType::Pointer resultImage = MaskImageType::New();
  GenerateBrainVolume(firstImage,
                      secondImage,
                      maskImage,
                      inputMaskVolume,
                      lowerThresholdFactor,
                      upperThresholdFactor,
                      boundingBoxSize,
                      boundingBoxStart,
                      outputVolume,
                      outputWeightsFile,
                      resultImage);

  MaskImageWriterType::Pointer maskWriter = MaskImageWriterType::New();
  maskWriter->SetInput(resultImage);
  maskWriter->SetFileName(outputMask);
  try
  {
    std::cout << "Writing mush mask : " << resultImage << std::endl;
    maskWriter->Update();
  }
  catch (const itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  return 0;

  // End of Output
}
