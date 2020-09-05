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
--desiredMean 10000 [optional] --desiredVariance 0 [optional] --seed "128,128,128" [optional]  --outputVolume
mush_2.nii.gz [optional] --outputMask mask_2.nii.gz [optional] --outputWeightsFile weights.txt [optional]
--boundingBoxSize "90,60,75" [optional] --boundingBoxStart "83,113,80" [optional]

Minimal Input Example:
./BRAINSMush --inputFirstVolume T1.nii.gz --inputSecondVolume T2.nii.gz
=========================================================================*/

#include "BRAINSMush.h"
#include "BRAINSMushCLP.h"
#include "BRAINSThreadControl.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkThresholdImageFilter.h"

#include "itkLargestForegroundFilledMaskImageFilter.h"

#include <fstream>
#include <string>
#include <BRAINSCommonLib.h>

#define PR(x)                                                                                                          \
  std::cout << #x " = " << x << "\n"; // a simple print macro for
                                      // use when debugging

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
    std::cout << "Desired Mean: " << desiredMean << std::endl;
    std::cout << "Desired Variance: " << desiredVariance << std::endl;

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
  ImageType::Pointer firstImage = LoadImage(inputFirstVolume.c_str());

  ImageType::Pointer secondImage = LoadImage(inputSecondVolume.c_str());

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
                        desiredMean,
                        desiredVariance,
                        lowerThresholdFactorPre,
                        upperThresholdFactorPre,
                        boundingBoxSize,
                        boundingBoxStart,
                        //      seed,
                        outputVolume,
                        //  outputMask,
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
                      desiredMean,
                      desiredVariance,
                      lowerThresholdFactor,
                      upperThresholdFactor,
                      boundingBoxSize,
                      boundingBoxStart,
                      //  seed,
                      outputVolume,
                      //    outputMask,
                      outputWeightsFile,
                      resultImage);

  MaskImageWriterType::Pointer maskWriter = MaskImageWriterType::New();
  maskWriter->SetInput(resultImage);
  maskWriter->SetFileName(outputMask);
  try
  {
    maskWriter->Update();
  }
  catch (itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  return 0;

  // End of Output
}

void
GenerateBrainVolume(ImageType::Pointer &     firstImage,
                    ImageType::Pointer &     secondImage,
                    MaskImageType::Pointer & maskImage,
                    std::string              inputMaskVolume,
                    double                   desiredMean,
                    double                   desiredVariance,
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
  ImageType::Pointer mixtureImage =
    MixtureOptimizer(firstImage, secondImage, maskImage, desiredMean, desiredVariance, outputWeightsFile);

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
    writer->Update();
  }
  catch (itk::ExceptionObject & exp)
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
      voxelCount++;
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
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Exception caught ! " << std::endl;
      std::cerr << excp << std::endl;
      // return EXIT_FAILURE;
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
    catch (itk::ExceptionObject & exp)
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
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Exception caught ! " << std::endl;
      std::cerr << excp << std::endl;
      // return EXIT_FAILURE;
    }

    // using StatisticRealType = LabelFilterType::RealType;
    mean = labelFilter->GetMean(1);

    // these definitions use magic numbers obtained through manual thresholding
    // and experimentation
    lower = (mean / lowerThresholdFactor);
    upper = (mean / upperThresholdFactor);
  }

  std::cout << "MushROI Mean:   " << mean
            << std::endl
            // << "MushROI StdDev: " <<   sigma << std::endl
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
  threshToHeadMask->SetLowerThreshold(lower);
  threshToHeadMask->SetUpperThreshold(upper);
  threshToHeadMask->SetInsideValue(1);
  threshToHeadMask->SetOutsideValue(0);

  try
  {
    threshToHeadMask->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Exception caught ! " << std::endl;
    std::cerr << excp << std::endl;
    // return EXIT_FAILURE;
  }

  using BinaryDilateFilterType = itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType>;

  using BinaryImageErodeFilterType = itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType>;

  double ClosingSize = 6;

  const ImageType::SpacingType & spacing = mixtureImage->GetSpacing();

  // Compute minumum object size as the number of voxels in the structuring
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
  catch (itk::ExceptionObject & exp)
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
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Exception caught ! " << std::endl;
    std::cerr << excp << std::endl;
    // return EXIT_FAILURE;
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
  catch (itk::ExceptionObject & exp)
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
  catch (itk::ExceptionObject & exp)
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
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Exception caught ! " << std::endl;
    std::cerr << excp << std::endl;
    // return EXIT_FAILURE;
  }

  MaskImageType::Pointer dilatedOutput = threshToClosureMask->GetOutput();
  dilatedOutput->DisconnectPipeline();

  // resultImage = FindLargestForgroundFilledMask<MaskImageType>( dilatedOutput,
  //                                                               0,
  //                                                               5 );
  using LFFMaskFilterType = itk::LargestForegroundFilledMaskImageFilter<MaskImageType>;
  LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
  LFF->SetInput(dilatedOutput);
  LFF->SetOtsuPercentileThreshold(0);
  LFF->SetClosingSize(5);
  LFF->Update();
  resultImage = LFF->GetOutput();
}

ImageType::Pointer
LoadImage(std::string imageName)
{
  ReaderType::Pointer loadImageReader = ReaderType::New();

  loadImageReader->SetFileName(imageName);
  try
  {
    loadImageReader->Update();
  }
  catch (itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  ImageType::Pointer image;

  return image = loadImageReader->GetOutput();
}

MaskImageType::Pointer
LoadMaskImage(std::string imageName)
{
  MaskReaderType::Pointer loadImageReader = MaskReaderType::New();

  loadImageReader->SetFileName(imageName);
  try
  {
    loadImageReader->Update();
  }
  catch (itk::ExceptionObject & exp)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }

  MaskImageType::Pointer image;

  return image = loadImageReader->GetOutput();
}

MaskImageType::Pointer
GenerateInitializerRegion(ImageType::Pointer & referenceImage,
                          std::vector<int>     boundingBoxSize,
                          std::vector<int>     boundingBoxStart)
{
  MaskImageType::Pointer initializeMask = MaskImageType::New();

  initializeMask->SetRegions(referenceImage->GetLargestPossibleRegion());
  // initializeMask->SetSpacing(referenceImage->GetSpacing());
  // initializeMask->SetOrigin(referenceImage->GetOrigin());
  // initializeMask->SetDirection(referenceImage->GetDirection());
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

ImageType::Pointer
MixtureOptimizer(ImageType::Pointer &     firstImage,
                 ImageType::Pointer &     secondImage,
                 MaskImageType::Pointer & maskImage,
                 double                   desiredMean,
                 double                   desiredVariance,
                 std::string              outputWeightsFile)
{
  using MixtureStatisticCostFunctionType = itk::MixtureStatisticCostFunction<ImageType, ImageType>;
  MixtureStatisticCostFunctionType::Pointer twoByTwoCostFunction = MixtureStatisticCostFunctionType::New();
  twoByTwoCostFunction->SetDesiredMean(desiredMean);
  twoByTwoCostFunction->SetDesiredVariance(desiredVariance);
  twoByTwoCostFunction->SetFirstImage(firstImage);
  twoByTwoCostFunction->SetSecondImage(secondImage);
  twoByTwoCostFunction->SetImageMask(maskImage);
  twoByTwoCostFunction->Initialize(1);

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Initialized MixtureStatisticCostFunction! " << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  double firstMean = twoByTwoCostFunction->GetSumOfFirstMaskVoxels() / twoByTwoCostFunction->GetNumberOfMaskVoxels();
  double secondMean = twoByTwoCostFunction->GetSumOfSecondMaskVoxels() / twoByTwoCostFunction->GetNumberOfMaskVoxels();
  double jointFactor = 1.0 / (firstMean + secondMean);

  using LevenbergMarquardtOptimizerType = itk::LevenbergMarquardtOptimizer;
  LevenbergMarquardtOptimizerType::Pointer twoByTwoOptimizer = LevenbergMarquardtOptimizerType::New();
  twoByTwoOptimizer->SetUseCostFunctionGradient(false);
  twoByTwoOptimizer->SetCostFunction(twoByTwoCostFunction);
  LevenbergMarquardtOptimizerType::ParametersType initialParameters(2);
  initialParameters[0] = firstMean * jointFactor;
  initialParameters[1] = secondMean * jointFactor;
  twoByTwoOptimizer->SetInitialPosition(initialParameters);

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Updating Levenberg-Marquardt Optimizer... " << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  try
  {
    twoByTwoOptimizer->StartOptimization();
  }
  catch (itk::ExceptionObject & exp)
  {
    std::cerr << "LMO FAIL" << std::endl;
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << exp << std::endl;
  }
  // End of Sending to Optimizer
  /* ------------------------------------------------------------------------------------
   */

  LevenbergMarquardtOptimizerType::ParametersType optimalParameters = twoByTwoOptimizer->GetCurrentPosition();
  LevenbergMarquardtOptimizerType::MeasureType    optimalMeasures = twoByTwoOptimizer->GetValue();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Obtained Output from Levenberg-Marquardt Optimizer!" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl << std::endl;

  /* ------------------------------------------------------------------------------------
   * Save the optimization:  one line with the image coeficients in order,
   * and one line with the error measure ito desired mean and variance.  Print
   * it all out, too.
   */
  const double firstWeight = optimalParameters[0];
  const double secondWeight = optimalParameters[1];

  std::cout << "First Weight:  " << firstWeight << std::endl;
  std::cout << "Second Weight:  " << secondWeight << std::endl;
  std::cout << "Optimality of Mean:  " << optimalMeasures[0] << std::endl;
  std::cout << "Optimality of Variance:  " << optimalMeasures[1] << std::endl << std::endl;

  // write a text file named outputWeightsFile

  std::ofstream to(outputWeightsFile.c_str());
  if (to.is_open())
  {
    to << firstWeight << "  " << secondWeight << std::endl;
    to << optimalMeasures[0] << "  " << optimalMeasures[1] << std::endl;
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
    PixelType firstValue = firstIt.Get();
    PixelType secondValue = secondIt.Get();

    double mixtureValue = firstWeight * firstValue + secondWeight * secondValue;

    mixtureIt.Set(mixtureValue);
  }

  return mixtureImage;
}
