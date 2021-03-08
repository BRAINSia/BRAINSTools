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
#include <sstream>
#include "itkMedianImageFilter.h"
#include "itkExtractImageFilter.h"
#include "BRAINSCommonLib.h"
#include "BRAINSThreadControl.h"
#include "BRAINSFitHelper.h"
#include "BRAINSFitCLP.h"

#include "BRAINSToolsVersion.h"


// This program was modified from
// Insight/Examples/Registration/ImageRegistration8.cxx
// and is an improved replacement for the old (and defective)

using BRAINSFitPixelType = float;
using FixedVolumeType = itk::Image<BRAINSFitPixelType, Dimension>;
using MovingVolumeType = itk::Image<BRAINSFitPixelType, Dimension>;

using InputImageType = itk::Image<BRAINSFitPixelType, MaxInputDimension>;
using FixedVolumeReaderType = itk::ImageFileReader<InputImageType>;
using MovingVolumeReaderType = itk::ImageFileReader<InputImageType>;
using AffineTransformPointer = itk::AffineTransform<double, 3>::Pointer;

template <typename ImageType>
typename ImageType::Pointer
ExtractImage(typename InputImageType::Pointer & inputImage, unsigned int InputImageTimeIndex)
{
  using ExtractImageFilterType = typename itk::ExtractImageFilter<InputImageType, ImageType>;
  typename ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
  extractImageFilter->SetDirectionCollapseToSubmatrix();

  // fixedVolumeReader->GetOutput();
  InputImageType::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
  InputImageType::SizeType   inputSize = inputRegion.GetSize();
  inputSize[3] = 0;
  inputRegion.SetSize(inputSize);

  InputImageType::IndexType inputIndex = inputRegion.GetIndex();
  inputIndex[0] = 0;
  inputIndex[1] = 0;
  inputIndex[2] = 0;
  inputIndex[3] = InputImageTimeIndex;
  inputRegion.SetIndex(inputIndex);
  extractImageFilter->SetExtractionRegion(inputRegion);
  extractImageFilter->SetInput(inputImage);

  try
  {
    extractImageFilter->Update();
  }
  catch (...)
  {
    std::cout << "Error while extracting a time indexed fixed image." << std::endl;
    throw;
  }

  typename ImageType::Pointer extractImage = extractImageFilter->GetOutput();
  //  std::cerr << "Extract fixed image origin" << extractImage->GetOrigin() << std::endl;

  return extractImage;
}

template <typename ImageType>
typename ImageType::Pointer
DoMedian(typename ImageType::Pointer & input, typename ImageType::SizeType indexRadius)
{
  using MedianFilterType = typename itk::MedianImageFilter<ImageType, ImageType>;
  typename MedianFilterType::Pointer medianFilter = MedianFilterType::New();

  medianFilter->SetRadius(indexRadius);
  medianFilter->SetInput(input);
  medianFilter->Update();
  typename ImageType::Pointer result = medianFilter->GetOutput();
  return result;
}

#ifdef USE_DebugImageViewer
/*************************
 * Have a global variable to
 * add debugging information.
 */
DebugImageViewerClient DebugImageDisplaySender;
#endif

int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  if (printVersionInfo)
  {
    std::cout << BRAINSTools::Version::ExtendedVersionString() << std::endl;
    return EXIT_SUCCESS;
  }

  BRAINSRegisterAlternateIO();

#ifdef USE_DebugImageViewer
  if (UseDebugImageViewer)
  {
    DebugImageDisplaySender.SetEnabled(UseDebugImageViewer);
  }
#endif
  using GenericTransformType = itk::Transform<double, 3, 3>;

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if (debugLevel > 1)
  {
    std::cout << "Number Of Threads used: " << numberOfThreads << std::endl;
    std::cout << " " << itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads() << std::endl;
  }

  std::string localInitializeTransformMode = initializeTransformMode;

  std::vector<int> BSplineGridSize(3, 0);
  {
    // Verify that spline grid size has enough parameters
    // If it has more than 3 parameters, only the first 3 parameters are used for three dimensions.
    // However, if it has less than 3 parameters, the first parameter is used for all three dimensions.
    if (splineGridSize.size() >= 3)
    {
      for (unsigned int sgs = 0; sgs < BSplineGridSize.size(); ++sgs)
      {
        // Verify that the spline grid sizes are greater than 3
        if (splineGridSize[sgs] < 3)
        {
          std::cout << "splineGridSize[" << sgs << "]= " << splineGridSize[sgs]
                    << " is invalid.  There must be at least 3 divisions in each dimension of the image." << std::endl;
          return EXIT_FAILURE;
        }
        BSplineGridSize[sgs] = splineGridSize[sgs];
      }
    }
    else
    {
      for (int & sgs : BSplineGridSize)
      {
        if (splineGridSize[0] < 3)
        {
          std::cout << "splineGridSize = " << splineGridSize[0]
                    << " is invalid.  There must be at least 3 divisions in each dimension of the image." << std::endl;
          return EXIT_FAILURE;
        }
        sgs = splineGridSize[0];
      }
    }
  }

  std::vector<std::string> localTransformType;
  // See if the individual boolean registration options are being used.  If any
  // of these are set, then transformType is not used.
  if ((useRigid == true) || (useScaleVersor3D == true) || (useScaleSkewVersor3D == true) || (useAffine == true) ||
      (useBSpline == true) || (useSyN == true))
  {
    localTransformType.resize(0); // Set to zero length
    if (useRigid == true)
    {
      localTransformType.emplace_back("Rigid");
    }
    if (useScaleVersor3D == true)
    {
      localTransformType.emplace_back("ScaleVersor3D");
    }
    if (useScaleSkewVersor3D == true)
    {
      localTransformType.emplace_back("ScaleSkewVersor3D");
    }
    if (useAffine == true)
    {
      localTransformType.emplace_back("Affine");
    }
    if (useBSpline == true)
    {
      localTransformType.emplace_back("BSpline");
    }
    if (useSyN == true)
    {
      localTransformType.emplace_back("SyN");
    }
    if (useComposite)
    {
      localTransformType.emplace_back("Composite3D");
    }
  }
  else if (!transformType.empty())
  {
    localTransformType = transformType;
  }
  else if ((!initialTransform.empty()) || (initializeTransformMode != "Off"))
  {
    // Only do the initialization phase;
  }
  else
  {
    std::cerr << "ERROR: No registration phases specified to perform!" << std::endl;
    return EXIT_FAILURE;
  }

  // In order to make the Slicer interface work, a few alternate command line
  // options need to available
  std::string localOutputTransform;
  if ((!linearTransform.empty() && !bsplineTransform.empty()) ||
      (!linearTransform.empty() && !outputTransform.empty()) || (!outputTransform.empty() && !bsplineTransform.empty()))
  {
    std::cout << "Error:  user can only specify one output transform type." << std::endl;
    return EXIT_FAILURE;
  }
  if (!linearTransform.empty())
  {
    localOutputTransform = linearTransform;
    if ((!localTransformType.empty()) && ((localTransformType[localTransformType.size() - 1] == "BSpline") ||
                                          (localTransformType[localTransformType.size() - 1] == "SyN")))
    {
      std::cout << "Error:  Linear transforms can not be used for BSpline or SyN registration!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (!bsplineTransform.empty())
  {
    localOutputTransform = bsplineTransform;
    if ((!localTransformType.empty()) && (localTransformType[localTransformType.size() - 1] != "BSpline"))
    {
      std::cout << "Error:  BSpline registrations require output transform to be of type BSpline!" << std::endl;
      return EXIT_FAILURE;
    }
    else if (localTransformType.empty())
    {
      std::cout << "Error:  Initializer only registrations require output transform to be of type Linear!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (!outputTransform.empty())
  {
    localOutputTransform = outputTransform;
  }

  if (localOutputTransform.empty() && strippedOutputTransform.empty() && outputVolume.empty())
  {
    std::cout << "Error:  user requested neither localOutputTransform,"
              << " nor strippedOutputTransform,"
              << " nor outputVolume." << std::endl;
    return 2;
  }

  if (numberOfIterations.size() != localTransformType.size())
  {
    if (numberOfIterations.size() != 1)
    {
      std::cerr << "The numberOfIterations array must match the length of the transformType"
                << "length, or have a single value that is used for all registration phases." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      // replicate throughout
      numberOfIterations.resize(localTransformType.size());
      const int numberOf = numberOfIterations[0];
      for (unsigned int i = 1; i < localTransformType.size(); i++)
      {
        numberOfIterations[i] = numberOf;
      }
    }
  }

  // Need to ensure that the order of transforms is from smallest to largest.
  if (localInitializeTransformMode == "Off")
  {
    try
    {
      itk::ValidateTransformRankOrdering(localTransformType);
    }
    catch (itk::ExceptionObject & err)
    {
      std::cout << "Exception Object caught: " << std::endl;
      std::cout << err << std::endl;
      throw;
    }
  }

  // Extracting a timeIndex cube from the fixed image goes here....
  // Also MedianFilter
  FixedVolumeType::Pointer  extractFixedVolume;
  MovingVolumeType::Pointer extractMovingVolume;
  InputImageType::Pointer   OriginalFixedVolume(itkUtil::ReadImage<InputImageType>(fixedVolume));

  std::cout << "Original Fixed image origin" << OriginalFixedVolume->GetOrigin() << std::endl;
  // fixedVolumeTimeIndex lets lets us treat 3D as 4D.
  /***********************
   * Acquire Fixed Image Index
   **********************/
  extractFixedVolume = ExtractImage<FixedVolumeType>(OriginalFixedVolume, fixedVolumeTimeIndex);
  // Extracting a timeIndex cube from the moving image goes here....

  InputImageType::Pointer OriginalMovingVolume(itkUtil::ReadImage<InputImageType>(movingVolume));
  // This default lets us treat 3D as 4D.
  // const unsigned int movingVolumeTimeIndex;

  /***********************
   * Acquire Moving Image Index
   **********************/
  extractMovingVolume = ExtractImage<MovingVolumeType>(OriginalMovingVolume, movingVolumeTimeIndex);
  // Multimodal registration input setting
  FixedVolumeType::Pointer  extractFixedVolume2 = nullptr;
  MovingVolumeType::Pointer extractMovingVolume2 = nullptr;
  if (!fixedVolume2.empty() && !movingVolume2.empty())
  {
    InputImageType::Pointer OriginalFixedVolume2(itkUtil::ReadImage<InputImageType>(fixedVolume2));
    std::cout << "Second Fixed image original origin" << OriginalFixedVolume2->GetOrigin() << std::endl;
    extractFixedVolume2 = ExtractImage<FixedVolumeType>(OriginalFixedVolume2, fixedVolumeTimeIndex);

    InputImageType::Pointer OriginalMovingVolume2(itkUtil::ReadImage<InputImageType>(movingVolume2));
    std::cout << "Second Moving image original origin" << OriginalMovingVolume2->GetOrigin() << std::endl;
    extractMovingVolume2 = ExtractImage<MovingVolumeType>(OriginalMovingVolume2, movingVolumeTimeIndex);
  }

#ifdef USE_DebugImageViewer
  if (DebugImageDisplaySender.Enabled())
  {
    DebugImageDisplaySender.SendImage<itk::Image<float, 3>>(extractFixedVolume, 0);
    DebugImageDisplaySender.SendImage<itk::Image<float, 3>>(extractMovingVolume, 1);
  }
#endif

  // get median filter radius.
  // const unsigned int MedianFilterRadius =
  // command.GetValueAsInt(MedianFilterRadiusText, IntegerText);
  // Median Filter images if requested.
  if ((medianFilterSize[0] > 0) || (medianFilterSize[1] > 0) || (medianFilterSize[2] > 0))
  {
    FixedVolumeType::SizeType indexRadius;
    indexRadius[0] = static_cast<long int>(medianFilterSize[0]); // radius
                                                                 // along x
    indexRadius[1] = static_cast<long int>(medianFilterSize[1]); // radius
                                                                 // along y
    indexRadius[2] = static_cast<long int>(medianFilterSize[2]); // radius
                                                                 // along z
    // DEBUG
    std::cout << "Median radius  " << indexRadius << std::endl;
    std::cout << "Fixed Image size     " << extractFixedVolume->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "Moving Image size     " << extractMovingVolume->GetLargestPossibleRegion().GetSize() << std::endl;
    extractFixedVolume = DoMedian<FixedVolumeType>(extractFixedVolume, indexRadius);
    extractMovingVolume = DoMedian<MovingVolumeType>(extractMovingVolume, indexRadius);
  }

  // If masks are associated with the images, then read them into the correct
  // orientation.
  // if they've been defined assign the masks...
  ImageMaskPointer fixedMask = nullptr;
  ImageMaskPointer movingMask = nullptr;
  if (maskProcessingMode == "NOMASK")
  {
    if (!fixedBinaryVolume.empty() || !movingBinaryVolume.empty())
    {
      std::cout
        << "ERROR:  Can not specify mask file names when the default of NOMASK is used for the maskProcessingMode"
        << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (maskProcessingMode == "ROIAUTO")
  {
    if (!fixedBinaryVolume.empty() || !movingBinaryVolume.empty())
    {
      std::cout << "ERROR:  Can not specify mask file names when ROIAUTO is used for the maskProcessingMode"
                << std::endl;
      return EXIT_FAILURE;
    }
    {
      using ROIAutoType = itk::BRAINSROIAutoImageFilter<FixedVolumeType, itk::Image<unsigned char, 3>>;
      ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
      ROIFilter->SetInput(extractFixedVolume);
      ROIFilter->SetClosingSize(ROIAutoClosingSize);
      ROIFilter->SetDilateSize(ROIAutoDilateSize);
      ROIFilter->Update();
      fixedMask = ROIFilter->GetSpatialObjectROI();
    }
    {
      using ROIAutoType = itk::BRAINSROIAutoImageFilter<MovingVolumeType, itk::Image<unsigned char, 3>>;
      ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
      ROIFilter->SetInput(extractMovingVolume);
      ROIFilter->SetClosingSize(ROIAutoClosingSize);
      ROIFilter->SetDilateSize(ROIAutoDilateSize);
      ROIFilter->Update();
      movingMask = ROIFilter->GetSpatialObjectROI();
    }
  }
  else if (maskProcessingMode == "ROI")
  {
    if (fixedBinaryVolume.empty() && movingBinaryVolume.empty())
    {
      std::cout << "ERROR:  Must specify mask file names when ROI is used for the maskProcessingMode" << std::endl;
      return EXIT_FAILURE;
    }
    if (!fixedBinaryVolume.empty())
    {
      fixedMask = ReadImageMask<SpatialObjectType, Dimension>(fixedBinaryVolume, extractFixedVolume.GetPointer());
    }
    else
    {
      fixedMask = nullptr;
    }

    if (!movingBinaryVolume.empty())
    {
      movingMask = ReadImageMask<SpatialObjectType, Dimension>(movingBinaryVolume, extractMovingVolume.GetPointer());
    }
    else
    {
      movingMask = nullptr;
    }
  }

  using CompositeTransformType = itk::CompositeTransform<double, 3>;
  CompositeTransformType::Pointer currentGenericTransform;
  if (!initialTransform.empty())
  {
    currentGenericTransform = CompositeTransformType::New();
    GenericTransformType::Pointer movingInitialTransform = itk::ReadTransformFromDisk(initialTransform);
    currentGenericTransform->AddTransform(movingInitialTransform);
  }

  FixedVolumeType::Pointer resampledImage;
  /*
   *  Everything prior to this point is preprocessing
   *  Start Processing
   *
   */

  {
    using HelperType = itk::BRAINSFitHelper;
    HelperType::Pointer myHelper = HelperType::New();
    myHelper->SetTransformType(localTransformType);
    myHelper->SetFixedVolume(extractFixedVolume);
    myHelper->SetMovingVolume(extractMovingVolume);
    if (extractFixedVolume2.IsNotNull() && extractMovingVolume2.IsNotNull())
    {
      myHelper->SetFixedVolume2(extractFixedVolume2);
      myHelper->SetMovingVolume2(extractMovingVolume2);
    }
    myHelper->SetHistogramMatch(histogramMatch);
    myHelper->SetRemoveIntensityOutliers(removeIntensityOutliers);
    myHelper->SetNumberOfMatchPoints(numberOfMatchPoints);
    myHelper->SetFixedBinaryVolume(fixedMask);
    myHelper->SetMovingBinaryVolume(movingMask);
    myHelper->SetOutputFixedVolumeROI(outputFixedVolumeROI);
    myHelper->SetOutputMovingVolumeROI(outputMovingVolumeROI);
    if (numberOfSamples > 0)
    {
      const unsigned long numberOfAllSamples = extractFixedVolume->GetBufferedRegion().GetNumberOfPixels();
      samplingPercentage = static_cast<double>(numberOfSamples) / numberOfAllSamples;
      std::cout << "WARNING --numberOfSamples is deprecated, please use --samplingPercentage instead " << std::endl;
      std::cout << "WARNING: Replacing command line --samplingPercentage " << samplingPercentage << std::endl;
    }
    myHelper->SetSamplingPercentage(samplingPercentage);
    myHelper->SetNumberOfHistogramBins(numberOfHistogramBins);
    myHelper->SetNumberOfIterations(numberOfIterations);
    myHelper->SetMaximumStepLength(maximumStepLength);
    myHelper->SetMinimumStepLength(minimumStepLength);
    myHelper->SetRelaxationFactor(relaxationFactor);
    myHelper->SetTranslationScale(translationScale);
    myHelper->SetReproportionScale(reproportionScale);
    myHelper->SetSkewScale(skewScale);
    myHelper->SetBackgroundFillValue(backgroundFillValue);
    myHelper->SetInitializeTransformMode(localInitializeTransformMode);
    myHelper->SetMaskInferiorCutOffFromCenter(maskInferiorCutOffFromCenter);
    myHelper->SetCurrentGenericTransform(currentGenericTransform);
    myHelper->SetSplineGridSize(BSplineGridSize);
    myHelper->SetCostFunctionConvergenceFactor(costFunctionConvergenceFactor);
    myHelper->SetProjectedGradientTolerance(projectedGradientTolerance);
    myHelper->SetMaxBSplineDisplacement(maxBSplineDisplacement);
    myHelper->SetDisplayDeformedImage(UseDebugImageViewer);
    myHelper->SetPromptUserAfterDisplay(PromptAfterImageSend);
    myHelper->SetDebugLevel(debugLevel);
    myHelper->SetCostMetricName(costMetric);
    myHelper->SetUseROIBSpline(useROIBSpline);
    myHelper->SetSamplingStrategy(metricSamplingStrategy);
    myHelper->SetInitializeRegistrationByCurrentGenericTransform(initializeRegistrationByCurrentGenericTransform);
    myHelper->SetMaximumNumberOfEvaluations(maximumNumberOfEvaluations);
    myHelper->SetMaximumNumberOfCorrections(maximumNumberOfCorrections);
    myHelper->SetWriteOutputTransformInFloat(writeOutputTransformInFloat);

    // HACK: create a flag for normalization
    bool NormalizeInputImages = false;
    myHelper->SetNormalizeInputImages(NormalizeInputImages);

    if (debugLevel > 7)
    {
      myHelper->PrintCommandLine(true, "BF");
    }
    try
    {
      myHelper->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cerr << "Exception during registration: " << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
    currentGenericTransform = myHelper->GetCurrentGenericTransform();

    std::string currentGenericTransformFileType;
    if (currentGenericTransform.IsNotNull())
    {
      currentGenericTransformFileType = std::string(currentGenericTransform->GetNameOfClass());
    }
    if (currentGenericTransformFileType != "CompositeTransform")
    {
      itkGenericExceptionMacro(<< "ERROR: Output transform is null.");
    }
    CompositeTransformType::Pointer outputComposite =
      static_cast<CompositeTransformType *>(currentGenericTransform.GetPointer());

    MovingVolumeType::ConstPointer preprocessedMovingVolume = myHelper->GetPreprocessedMovingVolume();
    if (NormalizeInputImages)
    {
      preprocessedMovingVolume = extractMovingVolume; // The resampled image should not be normalized
                                                      // because it my be casted to zero later.
    }

    {
      resampledImage = GenericTransformImage<MovingVolumeType>(preprocessedMovingVolume,
                                                               extractFixedVolume,
                                                               currentGenericTransform.GetPointer(),
                                                               backgroundFillValue,
                                                               interpolationMode,
                                                               false);
    }

    // If --logFileReport myReport.csv is specified on the command line, then write out this simple CSV file.
    if (!logFileReport.empty())
    {
      const double      finalMetricValue = myHelper->GetFinalMetricValue();
      std::stringstream myLogFileReportStream; // ( logFileReport );
      myLogFileReportStream << "#MetricName,MetricValue,FixedImageName,FixedMaskName,MovingImageName,MovingMaskName"
                            << std::endl;
      myLogFileReportStream << costMetric << ",";
      myLogFileReportStream << finalMetricValue << ",";
      myLogFileReportStream << fixedVolume << ",";
      myLogFileReportStream << fixedBinaryVolume << ",";
      myLogFileReportStream << movingVolume << ",";
      myLogFileReportStream << movingBinaryVolume << std::endl;

      std::ofstream LogScript;
      LogScript.open(logFileReport.c_str());
      if (!LogScript.is_open())
      {
        std::cerr << "Error: Can't write log file report file " << logFileReport << std::endl;
        std::cerr.flush();
        return EXIT_FAILURE;
      }
      LogScript << myLogFileReportStream.str();
      LogScript.close();
    }
  }
  /*
   *  At this point we can save the resampled image.
   */

  if (!outputVolume.empty())
  {
    //      std::cout << "=========== resampledImage :\n" <<
    // resampledImage->GetDirection() << std::endl;
    // Set in PARSEARGS const bool scaleOutputValues=false;//INFO: Make this a
    // command line parameter
    if (outputVolumePixelType == "float")
    {
      // itkUtil::WriteCastImage<itk::Image<float,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      using WriteOutImageType = itk::Image<float, FixedVolumeType::ImageDimension>;
      WriteOutImageType::Pointer CastImage =
        (scaleOutputValues == true) ? (itkUtil::PreserveCast<FixedVolumeType, WriteOutImageType>(resampledImage))
                                    : (itkUtil::TypeCast<FixedVolumeType, WriteOutImageType>(resampledImage));
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
    }
    else if (outputVolumePixelType == "short")
    {
      // itkUtil::WriteCastImage<itk::Image<signed short,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      using WriteOutImageType = itk::Image<signed short, FixedVolumeType::ImageDimension>;
      WriteOutImageType::Pointer CastImage =
        (scaleOutputValues == true) ? (itkUtil::PreserveCast<FixedVolumeType, WriteOutImageType>(resampledImage))
                                    : (itkUtil::TypeCast<FixedVolumeType, WriteOutImageType>(resampledImage));
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
    }
    else if (outputVolumePixelType == "ushort")
    {
      // itkUtil::WriteCastImage<itk::Image<unsigned short,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      using WriteOutImageType = itk::Image<unsigned short, FixedVolumeType::ImageDimension>;
      WriteOutImageType::Pointer CastImage =
        (scaleOutputValues == true) ? (itkUtil::PreserveCast<FixedVolumeType, WriteOutImageType>(resampledImage))
                                    : (itkUtil::TypeCast<FixedVolumeType, WriteOutImageType>(resampledImage));
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
    }
    else if (outputVolumePixelType == "int")
    {
      // itkUtil::WriteCastImage<itk::Image<signed int,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      using WriteOutImageType = itk::Image<signed int, FixedVolumeType::ImageDimension>;
      WriteOutImageType::Pointer CastImage =
        (scaleOutputValues == true) ? (itkUtil::PreserveCast<FixedVolumeType, WriteOutImageType>(resampledImage))
                                    : (itkUtil::TypeCast<FixedVolumeType, WriteOutImageType>(resampledImage));
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
    }
    else if (outputVolumePixelType == "uint")
    {
      // itkUtil::WriteCastImage<itk::Image<unsigned int,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      using WriteOutImageType = itk::Image<unsigned int, FixedVolumeType::ImageDimension>;
      WriteOutImageType::Pointer CastImage =
        (scaleOutputValues == true) ? (itkUtil::PreserveCast<FixedVolumeType, WriteOutImageType>(resampledImage))
                                    : (itkUtil::TypeCast<FixedVolumeType, WriteOutImageType>(resampledImage));
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
    }
    else if (outputVolumePixelType == "uchar")
    {
      // itkUtil::WriteCastImage<itk::Image<unsigned char,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      using WriteOutImageType = itk::Image<unsigned char, FixedVolumeType::ImageDimension>;
      WriteOutImageType::Pointer CastImage =
        (scaleOutputValues == true) ? (itkUtil::PreserveCast<FixedVolumeType, WriteOutImageType>(resampledImage))
                                    : (itkUtil::TypeCast<FixedVolumeType, WriteOutImageType>(resampledImage));
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
    }
  }
  if (writeOutputTransformInFloat)
  {
    std::cout << "Write the output transform in single (float) precision..." << std::endl;
    itk::WriteBothTransformsToDisk<double, float>(
      currentGenericTransform.GetPointer(), localOutputTransform, strippedOutputTransform);
  }
  else
  {
    itk::WriteBothTransformsToDisk<double, double>(
      currentGenericTransform.GetPointer(), localOutputTransform, strippedOutputTransform);
  }


  return 0;
}
