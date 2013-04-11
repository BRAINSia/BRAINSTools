/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile$
Language:  C++
Date:      $Date: 2007-08-31 11:20:20 -0500 (Fri, 31 Aug 2007) $
Version:   $Revision: 10358 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <sstream>
#include "itkMedianImageFilter.h"
#include "itkExtractImageFilter.h"
#include "BRAINSCommonLib.h"
#include "BRAINSThreadControl.h"
#include "BRAINSFitHelper.h"
#include "BRAINSFitCLP.h"

// This program was modified from
// Insight/Examples/Registration/ImageRegistration8.cxx
// and is an improved replacement for the old (and defective)

typedef float                                     BRAINSFitPixelType;
typedef itk::Image<BRAINSFitPixelType, Dimension> FixedVolumeType;
typedef itk::Image<BRAINSFitPixelType, Dimension> MovingVolumeType;

typedef itk::Image<BRAINSFitPixelType, MaxInputDimension> InputImageType;
typedef itk::ImageFileReader<InputImageType>              FixedVolumeReaderType;
typedef itk::ImageFileReader<InputImageType>              MovingVolumeReaderType;
typedef AffineTransformType::Pointer                      AffineTransformPointer;
// typedef itk::Vector<double, Dimension>           BRAINSFitVectorType;

template <class ImageType>
typename ImageType::Pointer ExtractImage(
  typename InputImageType::Pointer & inputImage,
  unsigned int InputImageTimeIndex)
{
  typedef typename itk::ExtractImageFilter<InputImageType, ImageType> ExtractImageFilterType;
  typename ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
#if  ITK_VERSION_MAJOR >= 4
  extractImageFilter->SetDirectionCollapseToSubmatrix();
#endif

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
  catch( ... )
    {
    std::cout << "Error while extracting a time indexed fixed image."
              << std::endl;
    throw;
    }

  typename ImageType::Pointer extractImage = extractImageFilter->GetOutput();
  //  std::cerr << "Extract fixed image origin" << extractImage->GetOrigin() <<
  // std::endl;

  return extractImage;
}

template <class ImageType>
typename ImageType::Pointer DoMedian(typename ImageType::Pointer & input,
                                     typename ImageType::SizeType indexRadius)
{
  typedef typename itk::MedianImageFilter<ImageType,
                                          ImageType> MedianFilterType;
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

int main(int argc, char *argv[])
{
  PARSE_ARGS;

#ifdef USE_DebugImageViewer
  if( UseDebugImageViewer )
    {
    DebugImageDisplaySender.SetEnabled(UseDebugImageViewer);
    }
#endif

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  std::string localInitializeTransformMode = initializeTransformMode;

  // std::vector<int> zeroOrigin(3, 0);

  // Verify that the spline grid sizes are greater than 3
    {
//    bool validGridSize = true;
    for( unsigned int sgs = 0; sgs < splineGridSize.size(); sgs++ )
      {
      if( splineGridSize[sgs] < 3 )
        {
//        validGridSize = false;
        std::cout << "splineGridSize[" << sgs << "]= " << splineGridSize[sgs]
                  << " is invalid.  There must be at lest 3 divisions in each dimension of the image." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  std::vector<std::string> localTransformType;
  // See if the individual boolean registration options are being used.  If any
  // of these are set, then transformType is not used.
  if( ( useRigid == true ) || ( useScaleVersor3D == true ) || ( useScaleSkewVersor3D == true )
      || ( useAffine == true ) || ( useBSpline == true ) || ( useSyN == true ) )
    {
    localTransformType.resize(0); // Set to zero length
    if( useRigid == true )
      {
      localTransformType.push_back("Rigid");
      }
    if( useScaleVersor3D == true )
      {
      localTransformType.push_back("ScaleVersor3D");
      }
    if( useScaleSkewVersor3D == true )
      {
      localTransformType.push_back("ScaleSkewVersor3D");
      }
    if( useAffine == true )
      {
      localTransformType.push_back("Affine");
      }
    if( useBSpline == true )
      {
      localTransformType.push_back("BSpline");
      }
    if( useSyN == true )
      {
      localTransformType.push_back("SyN");
      }
    if( useComposite )
      {
      localTransformType.push_back("Composite3D");
      }
    }
  else if( transformType.size() > 0 )
    {
    localTransformType = transformType;
    }
  else if( (initialTransform.size() > 0) || (initializeTransformMode != "Off") )
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
  if( ( linearTransform.size() > 0 && bsplineTransform.size() > 0 )
      || ( linearTransform.size() > 0 && outputTransform.size() > 0 )
      || ( outputTransform.size() > 0 && bsplineTransform.size() > 0 ) )
    {
    std::cout << "Error:  user can only specify one output transform type." << std::endl;
    return EXIT_FAILURE;
    }
  if( linearTransform.size() > 0 )
    {
    localOutputTransform = linearTransform;
    if( ( !localTransformType.empty() ) &&
        ( (localTransformType[localTransformType.size() - 1] == "BSpline") ||
          (localTransformType[localTransformType.size() - 1] == "SyN") ) )
      {
      std::cout << "Error:  Linear transforms can not be used for BSpline or SyN registration!" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if( !bsplineTransform.empty() )
    {
    localOutputTransform = bsplineTransform;
    if( ( !localTransformType.empty() ) && ( localTransformType[localTransformType.size() - 1] != "BSpline" ) )
      {
      std::cout << "Error:  BSpline registrations require output transform to be of type BSpline!" << std::endl;
      return EXIT_FAILURE;
      }
    else if( localTransformType.empty() )
      {
      std::cout << "Error:  Initializer only registrations require output transform to be of type Linear!" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if( !outputTransform.empty() )
    {
    localOutputTransform = outputTransform;
    }

  if( localOutputTransform.empty()
      && strippedOutputTransform.empty()
      && outputVolume.empty() )
    {
    std::cout << "Error:  user requested neither localOutputTransform,"
              << " nor strippedOutputTransform,"
              << " nor outputVolume." << std::endl;
    return 2;
    }

  if( numberOfIterations.size() != localTransformType.size() )
    {
    if( numberOfIterations.size() != 1 )
      {
      std::cerr << "The numberOfIterations array must match the length of the transformType"
                << "length, or have a single value that is used for all registration phases." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      // replicate throughout
      numberOfIterations.resize( localTransformType.size() );
      const int numberOf = numberOfIterations[0];
      for( unsigned int i = 1; i < localTransformType.size(); i++ )
        {
        numberOfIterations[i] = numberOf;
        }
      }
    }
  if( minimumStepLength.size() != localTransformType.size() )
    {
    if( minimumStepLength.size() != 1 )
      {
      std::cerr << "The minimumStepLength array must match the localTransformType length" << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      // replicate throughout
      const double stepSize = minimumStepLength[0];
      for( unsigned int i = 1; i < localTransformType.size(); i++ )
        {
        minimumStepLength.push_back(stepSize);
        }
      }
    }

  // Need to ensure that the order of transforms is from smallest to largest.
  try
    {
    itk::ValidateTransformRankOrdering(localTransformType);
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    throw;
    }

  // Extracting a timeIndex cube from the fixed image goes here....
  // Also MedianFilter
  FixedVolumeType::Pointer  extractFixedVolume;
  MovingVolumeType::Pointer extractMovingVolume;
  InputImageType::Pointer
    OriginalFixedVolume( itkUtil::ReadImage<InputImageType>(fixedVolume) );

  std::cout << "Original Fixed image origin"
            << OriginalFixedVolume->GetOrigin() << std::endl;
  // fixedVolumeTimeIndex lets lets us treat 3D as 4D.
  /***********************
   * Acquire Fixed Image Index
   **********************/
  extractFixedVolume = ExtractImage<FixedVolumeType>(OriginalFixedVolume,
                                                     fixedVolumeTimeIndex);
  // Extracting a timeIndex cube from the moving image goes here....

  InputImageType::Pointer OriginalMovingVolume(
    itkUtil::ReadImage<InputImageType>(movingVolume) );
  // This default lets us treat 3D as 4D.
  // const unsigned int movingVolumeTimeIndex;

  /***********************
   * Acquire Moving Image Index
   **********************/
  extractMovingVolume = ExtractImage<MovingVolumeType>(OriginalMovingVolume,
                                                       movingVolumeTimeIndex);

#ifdef USE_DebugImageViewer
  if( DebugImageDisplaySender.Enabled() )
    {
    DebugImageDisplaySender.SendImage<itk::Image<float, 3> >(extractFixedVolume, 0);
    DebugImageDisplaySender.SendImage<itk::Image<float, 3> >(extractMovingVolume, 1);
    }
#endif

  // get median filter radius.
  // const unsigned int MedianFilterRadius =
  // command.GetValueAsInt(MedianFilterRadiusText, IntegerText);
  // Median Filter images if requested.
  if( (medianFilterSize[0] > 0) || (medianFilterSize[1] > 0)
      || (medianFilterSize[2] > 0) )
    {
    FixedVolumeType::SizeType indexRadius;
    indexRadius[0] = static_cast<long int>( medianFilterSize[0] );   // radius
                                                                     // along x
    indexRadius[1] = static_cast<long int>( medianFilterSize[1] );   // radius
                                                                     // along y
    indexRadius[2] = static_cast<long int>( medianFilterSize[2] );   // radius
                                                                     // along z
    // DEBUG
    std::cout << "Median radius  " << indexRadius << std::endl;
    std::cout << "Fixed Image size     "  << extractFixedVolume->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "Moving Image size     " << extractMovingVolume->GetLargestPossibleRegion().GetSize() << std::endl;
    extractFixedVolume = DoMedian<FixedVolumeType>(extractFixedVolume,
                                                   indexRadius);
    extractMovingVolume = DoMedian<MovingVolumeType>(extractMovingVolume,
                                                     indexRadius);
    }

  // If masks are associated with the images, then read them into the correct
  // orientation.
  // if they've been defined assign the masks...
  ImageMaskPointer fixedMask = NULL;
  ImageMaskPointer movingMask = NULL;
  if( maskProcessingMode == "NOMASK" )
    {
    if( fixedBinaryVolume != "" || movingBinaryVolume != "" )
      {
      std::cout
        << "ERROR:  Can not specify mask file names when the default of NOMASK is used for the maskProcessingMode"
        << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if( maskProcessingMode == "ROIAUTO" )
    {
    if( fixedBinaryVolume != "" || movingBinaryVolume != "" )
      {
      std::cout
        << "ERROR:  Can not specify mask file names when ROIAUTO is used for the maskProcessingMode"
        << std::endl;
      return EXIT_FAILURE;
      }
      {
      typedef itk::BRAINSROIAutoImageFilter<FixedVolumeType, itk::Image<unsigned char, 3> > ROIAutoType;
      ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
      ROIFilter->SetInput(extractFixedVolume);
      ROIFilter->SetClosingSize(ROIAutoClosingSize);
      ROIFilter->SetDilateSize(ROIAutoDilateSize);
      ROIFilter->Update();
      fixedMask = ROIFilter->GetSpatialObjectROI();
      }
      {
      typedef itk::BRAINSROIAutoImageFilter<MovingVolumeType, itk::Image<unsigned char, 3> > ROIAutoType;
      ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
      ROIFilter->SetInput(extractMovingVolume);
      ROIFilter->SetClosingSize(ROIAutoClosingSize);
      ROIFilter->SetDilateSize(ROIAutoDilateSize);
      ROIFilter->Update();
      movingMask = ROIFilter->GetSpatialObjectROI();
      }
    }
  else if( maskProcessingMode == "ROI" )
    {
    if( fixedBinaryVolume == "" && movingBinaryVolume == "" )
      {
      std::cout
        <<
        "ERROR:  Must specify mask file names when ROI is used for the maskProcessingMode"
        << std::endl;
      return EXIT_FAILURE;
      }
    if( fixedBinaryVolume != "" )
      {
      fixedMask = ReadImageMask<SpatialObjectType, Dimension>(
          fixedBinaryVolume,
          extractFixedVolume.GetPointer() );
      }
    else
      {
      fixedMask = NULL;
      }

    if( movingBinaryVolume != "" )
      {
      movingMask = ReadImageMask<SpatialObjectType, Dimension>(
          movingBinaryVolume,
          extractMovingVolume.GetPointer() );
      }
    else
      {
      movingMask = NULL;
      }
    }
  GenericTransformType::Pointer currentGenericTransform;
  if( initialTransform != "" )
    {
    currentGenericTransform = itk::ReadTransformFromDisk(initialTransform);
    }

  FixedVolumeType::Pointer resampledImage;
  /*
   *  Everything prior to this point is preprocessing
   *  Start Processing
   *
   */
//  int actualIterations = 0;
//  int permittedIterations = 0;
//  int allLevelsIterations = 0;

    {
    typedef itk::BRAINSFitHelper HelperType;
    HelperType::Pointer myHelper = HelperType::New();
    myHelper->SetTransformType(localTransformType);
    myHelper->SetFixedVolume(extractFixedVolume);
    myHelper->SetForceMINumberOfThreads(forceMINumberOfThreads);
    myHelper->SetMovingVolume(extractMovingVolume);
    myHelper->SetHistogramMatch(histogramMatch);
    myHelper->SetRemoveIntensityOutliers(removeIntensityOutliers);
    myHelper->SetNumberOfMatchPoints(numberOfMatchPoints);
    myHelper->SetFixedBinaryVolume(fixedMask);
    myHelper->SetMovingBinaryVolume(movingMask);
    myHelper->SetOutputFixedVolumeROI(outputFixedVolumeROI);
    myHelper->SetOutputMovingVolumeROI(outputMovingVolumeROI);
    myHelper->SetPermitParameterVariation(permitParameterVariation);
    myHelper->SetNumberOfSamples(numberOfSamples);
    myHelper->SetNumberOfHistogramBins(numberOfHistogramBins);
    myHelper->SetNumberOfIterations(numberOfIterations);
    myHelper->SetMaximumStepLength(maximumStepLength);
    myHelper->SetMinimumStepLength(minimumStepLength);
    myHelper->SetRelaxationFactor(relaxationFactor);
    myHelper->SetTranslationScale(translationScale);
    myHelper->SetReproportionScale(reproportionScale);
    myHelper->SetSkewScale(skewScale);
    myHelper->SetUseExplicitPDFDerivativesMode(useExplicitPDFDerivativesMode);
    myHelper->SetUseCachingOfBSplineWeightsMode(useCachingOfBSplineWeightsMode);
    myHelper->SetBackgroundFillValue(backgroundFillValue);
    myHelper->SetInitializeTransformMode(localInitializeTransformMode);
    myHelper->SetMaskInferiorCutOffFromCenter(maskInferiorCutOffFromCenter);
    myHelper->SetCurrentGenericTransform(currentGenericTransform);
    myHelper->SetSplineGridSize(splineGridSize);
    myHelper->SetCostFunctionConvergenceFactor(costFunctionConvergenceFactor);
    myHelper->SetProjectedGradientTolerance(projectedGradientTolerance);
    myHelper->SetMaxBSplineDisplacement(maxBSplineDisplacement);
    myHelper->SetDisplayDeformedImage(UseDebugImageViewer);
    myHelper->SetPromptUserAfterDisplay(PromptAfterImageSend);
    myHelper->SetDebugLevel(debugLevel);
    myHelper->SetCostMetric(costMetric);
    myHelper->SetUseROIBSpline(useROIBSpline);
    if( debugLevel > 7 )
      {
      myHelper->PrintCommandLine(true, "BF");
      }
    myHelper->Update();
    currentGenericTransform = myHelper->GetCurrentGenericTransform();
    MovingVolumeType::ConstPointer preprocessedMovingVolume = myHelper->GetPreprocessedMovingVolume();
#if 0
    if( interpolationMode == "ResampleInPlace" )
      {
      % {
        VersorRigid3DTransformType::ConstPointer versor3D =
          dynamic_cast<const VersorRigid3DTransformType *>(currentGenericTransform.GetPointer() );
        if( versor3D.IsNotNull() )
          {
          FixedVolumeType::Pointer tempInPlaceResample = itk::SetRigidTransformInPlace<FixedVolumeType>(
              versor3D.GetPointer(), extractMovingVolume.GetPointer() );
          resampledImage = itkUtil::TypeCast<FixedVolumeType, MovingVolumeType>(tempInPlaceResample);
          }
        else
          {
          // This should be an exception thow instead of exit.
          std::cout << "could not convert to rigid versor type" << std::endl;
          return EXIT_FAILURE;
          }
        }
      else
        {
        // Remember:  the Data is Moving's, the shape is Fixed's.
        resampledImage = TransformResample<MovingVolumeType, FixedVolumeType>(
            preprocessedMovingVolume,
            extractFixedVolume,
            backgroundFillValue,
            GetInterpolatorFromString<MovingVolumeType>(interpolationMode),
            currentGenericTransform);
        }
      }
#else
      {
      typedef float                                                                     VectorComponentType;
      typedef itk::Vector<VectorComponentType, GenericTransformImageNS::SpaceDimension> VectorPixelType;
      typedef itk::Image<VectorPixelType,  GenericTransformImageNS::SpaceDimension>     DisplacementFieldType;
      resampledImage = GenericTransformImage<MovingVolumeType, FixedVolumeType, DisplacementFieldType>(
          preprocessedMovingVolume,
          extractFixedVolume,
          NULL,
          currentGenericTransform,
          backgroundFillValue,
          interpolationMode,
          false);
      }
#endif
//    actualIterations = myHelper->GetActualNumberOfIterations();
//    permittedIterations = myHelper->GetPermittedNumberOfIterations();
      /*
      allLevelsIterations=myHelper->GetAccumulatedNumberOfIterationsForAllLevels();
      */

      //If --logFileReport myReport.csv is specified on the command line, then write out this simple CSV file.
    if( logFileReport != "" )
      {
      const double finalMetricValue = myHelper->GetFinalMetricValue();
      std::stringstream myLogFileReportStream; // ( logFileReport );
      myLogFileReportStream << "#MetricName,MetricValue,FixedImageName,FixedMaskName,MovingImageName,MovingMaskName" << std::endl;
      myLogFileReportStream << costMetric << ",";
      myLogFileReportStream << finalMetricValue << ",";
      myLogFileReportStream << fixedVolume << ",";
      myLogFileReportStream << fixedBinaryVolume << ",";
      myLogFileReportStream << movingVolume << ",";
      myLogFileReportStream << movingBinaryVolume << std::endl;

      std::ofstream     LogScript;
      LogScript.open( logFileReport.c_str() );
      if( !LogScript.is_open() )
        {
        std::cerr << "Error: Can't write log file report file "
        << logFileReport << std::endl;
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

  if( outputVolume.size() > 0 )
    {
    //      std::cout << "=========== resampledImage :\n" <<
    // resampledImage->GetDirection() << std::endl;
    // Set in PARSEARGS const bool scaleOutputValues=false;//TODO: Make this a
    // command line parameter
    if( outputVolumePixelType == "float" )
      {
      // itkUtil::WriteCastImage<itk::Image<float,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      typedef itk::Image<float,
                         FixedVolumeType::ImageDimension> WriteOutImageType;
      WriteOutImageType::Pointer CastImage =
        ( scaleOutputValues == true ) ?
        ( itkUtil::PreserveCast<FixedVolumeType,
                                WriteOutImageType>(resampledImage) ) :
        ( itkUtil::TypeCast<FixedVolumeType,
                            WriteOutImageType>(resampledImage) );
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
      }
    else if( outputVolumePixelType == "short" )
      {
      // itkUtil::WriteCastImage<itk::Image<signed short,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      typedef itk::Image<signed short,
                         FixedVolumeType::ImageDimension> WriteOutImageType;
      WriteOutImageType::Pointer CastImage =
        ( scaleOutputValues == true ) ?
        ( itkUtil::PreserveCast<FixedVolumeType,
                                WriteOutImageType>(resampledImage) ) :
        ( itkUtil::TypeCast<FixedVolumeType,
                            WriteOutImageType>(resampledImage) );
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
      }
    else if( outputVolumePixelType == "ushort" )
      {
      // itkUtil::WriteCastImage<itk::Image<unsigned short,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      typedef itk::Image<unsigned short,
                         FixedVolumeType::ImageDimension> WriteOutImageType;
      WriteOutImageType::Pointer CastImage =
        ( scaleOutputValues == true ) ?
        ( itkUtil::PreserveCast<FixedVolumeType,
                                WriteOutImageType>(resampledImage) ) :
        ( itkUtil::TypeCast<FixedVolumeType,
                            WriteOutImageType>(resampledImage) );
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
      }
    else if( outputVolumePixelType == "int" )
      {
      // itkUtil::WriteCastImage<itk::Image<signed int,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      typedef itk::Image<signed int,
                         FixedVolumeType::ImageDimension> WriteOutImageType;
      WriteOutImageType::Pointer CastImage =
        ( scaleOutputValues == true ) ?
        ( itkUtil::PreserveCast<FixedVolumeType,
                                WriteOutImageType>(resampledImage) ) :
        ( itkUtil::TypeCast<FixedVolumeType,
                            WriteOutImageType>(resampledImage) );
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
      }
    else if( outputVolumePixelType == "uint" )
      {
      // itkUtil::WriteCastImage<itk::Image<unsigned int,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      typedef itk::Image<unsigned int,
                         FixedVolumeType::ImageDimension> WriteOutImageType;
      WriteOutImageType::Pointer CastImage =
        ( scaleOutputValues == true ) ?
        ( itkUtil::PreserveCast<FixedVolumeType,
                                WriteOutImageType>(resampledImage) ) :
        ( itkUtil::TypeCast<FixedVolumeType,
                            WriteOutImageType>(resampledImage) );
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
      }
#if 0
    else if( outputVolumePixelType == "char" )
      {
      // itkUtil::WriteCastImage<itk::Image<signed char,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      typedef itk::Image<char,
                         FixedVolumeType::ImageDimension> WriteOutImageType;
      WriteOutImageType::Pointer CastImage =
        ( scaleOutputValues == true ) ?
        ( itkUtil::PreserveCast<FixedVolumeType,
                                WriteOutImageType>(resampledImage) ) :
        ( itkUtil::TypeCast<FixedVolumeType,
                            WriteOutImageType>(resampledImage) );
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
      }
#endif
    else if( outputVolumePixelType == "uchar" )
      {
      // itkUtil::WriteCastImage<itk::Image<unsigned char,
      // FixedVolumeType::ImageDimension>,
      // FixedVolumeType>(resampledImage,outputVolume);
      typedef itk::Image<unsigned char,
                         FixedVolumeType::ImageDimension> WriteOutImageType;
      WriteOutImageType::Pointer CastImage =
        ( scaleOutputValues == true ) ?
        ( itkUtil::PreserveCast<FixedVolumeType,
                                WriteOutImageType>(resampledImage) ) :
        ( itkUtil::TypeCast<FixedVolumeType,
                            WriteOutImageType>(resampledImage) );
      itkUtil::WriteImage<WriteOutImageType>(CastImage, outputVolume);
      }
    }

#ifdef USE_ANTS
  if( localTransformType[localTransformType.size() - 1] == "SyN" )
    {
    CompositeTransformType::Pointer tempSyNCompositeTransform =
      dynamic_cast<CompositeTransformType *>( currentGenericTransform.GetPointer() );
    // write out transform actually computed, so skip the initial transform
    CompositeTransformType::TransformTypePointer tempSyNFinalTransform =
      tempSyNCompositeTransform.GetPointer();
// tempSyNCompositeTransform->GetNthTransform( 1 );

    if( tempSyNFinalTransform.IsNull() )
      {
      std::cout << "Error in type conversion" << __FILE__ << __LINE__ << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      itk::ants::WriteTransform<3>( tempSyNFinalTransform, localOutputTransform );
      std::cout << "SyN warped transform is written to the disk." << std::endl;
      }
    }
  else
#endif
    {
    /*const int write_status=*/
    itk::WriteBothTransformsToDisk(currentGenericTransform.GetPointer(),
                                   localOutputTransform, strippedOutputTransform);
    }

  return 0;
}
