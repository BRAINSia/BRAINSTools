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

#include "BRAINSFitUtils.h"
#include "BRAINSFitHelper.h"
#include "BRAINSFitHelperTemplate.h"

#include "genericRegistrationHelper.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkKullbackLeiblerCompareHistogramImageToImageMetric.h"
#include "itkHistogramImageToImageMetric.h"
#include "itkKappaStatisticImageToImageMetric.h"
#include "itkMeanReciprocalSquareDifferenceImageToImageMetric.h"
#include "itkJointHistogramMutualInformationImageToImageMetricv4.h"
#include "itkGradientDifferenceImageToImageMetric.h"
#include "itkCompareHistogramImageToImageMetric.h"
#include "itkCorrelationCoefficientHistogramImageToImageMetric.h"
#include "itkMatchCardinalityImageToImageMetric.h"
#include "itkMeanSquaresHistogramImageToImageMetric.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"

// A little dummy function to make it easy to stop the debugger.
void debug_catch(void)
{
  std::cout << "HERE" << __FILE__ << " " << __LINE__ << std::endl;

  return;
}

// convert spatial object to image
MaskImageType::ConstPointer
ExtractConstPointerToImageMaskFromImageSpatialObject( SpatialObjectType::ConstPointer inputSpatialObject )
{
  ImageMaskSpatialObjectType const * const temp =
    dynamic_cast<ImageMaskSpatialObjectType const *>( inputSpatialObject.GetPointer() );

  if( temp == NULL )
    {
    itkGenericExceptionMacro(<< "Invalid mask converstion attempted.");
    }
  ImageMaskSpatialObjectType::ConstPointer ImageMask( temp );
  const MaskImageType::ConstPointer        tempOutputVolumeROI = ImageMask->GetImage();
  return tempOutputVolumeROI;
}

// convert image to mask (spatial object)
SpatialObjectType::ConstPointer
ConvertMaskImageToSpatialMask( MaskImageType::ConstPointer inputImage )
{
  ImageMaskSpatialObjectType::Pointer mask = ImageMaskSpatialObjectType::New();
  mask->SetImage(inputImage);
  mask->ComputeObjectToWorldTransform();
  // return pointer to mask
  SpatialObjectType::Pointer p = dynamic_cast<SpatialObjectType *>( mask.GetPointer() );
  if( p.IsNull() )
    {
    itkGenericExceptionMacro(<< "Failed conversion to Mask");
    }

  SpatialObjectType::ConstPointer objectMask(p);
  return objectMask;
}

namespace itk
{
BRAINSFitHelper::BRAINSFitHelper() :
  m_FixedVolume(NULL),
  m_MovingVolume(NULL),
  m_PreprocessedMovingVolume(NULL),
  m_FixedBinaryVolume(NULL),
  m_MovingBinaryVolume(NULL),
  m_OutputFixedVolumeROI(""),
  m_OutputMovingVolumeROI(""),
  m_PermitParameterVariation(0),
  m_NumberOfSamples(0), // normally this variable should NOT be used; however, it is kept for backward compatiblity.
  m_SamplingPercentage(1), // instead or number of samples, sampling% should be used that is a number between 0 and 1.
  m_NumberOfHistogramBins(50),
  m_HistogramMatch(false),
  m_RemoveIntensityOutliers(0.00),
  m_NumberOfMatchPoints(10),
  m_NumberOfIterations(1, 1500),
  m_MaximumStepLength(0.2),
  m_MinimumStepLength(1, 0.005),
  m_RelaxationFactor(0.5),
  m_TranslationScale(1000.0),
  m_ReproportionScale(1.0),
  m_SkewScale(1.0),
  m_UseCachingOfBSplineWeightsMode("ON"),
  m_BackgroundFillValue(0.0),
  m_TransformType(1, "Rigid"),
  m_InitializeTransformMode("Off"),
  m_MaskInferiorCutOffFromCenter(1000),
  m_SplineGridSize(3, 10),
  m_CostFunctionConvergenceFactor(1e+9),
  m_ProjectedGradientTolerance(1e-5),
  m_MaxBSplineDisplacement(0.0),
  m_ActualNumberOfIterations(0),
  m_PermittedNumberOfIterations(0),
  // m_AccumulatedNumberOfIterationsForAllLevels(0),
  m_DebugLevel(0),
  m_CurrentGenericTransform(NULL),
  //m_GenericTransformList(0),
  m_DisplayDeformedImage(false),
  m_PromptUserAfterDisplay(false),
  m_FinalMetricValue(0.0),
  m_ObserveIterations(true),
  m_CostMetric("MMI"), // Default to Mattes Mutual Information Metric
  m_UseROIBSpline(false),
  m_Helper(NULL),
  m_SamplingStrategy(AffineRegistrationType::NONE),
  m_NormalizeInputImages(false),
  m_InitializeRegistrationByCurrentGenericTransform(true),
  m_ForceMINumberOfThreads(-1)
{
  m_SplineGridSize[0] = 14;
  m_SplineGridSize[1] = 10;
  m_SplineGridSize[2] = 12;
}

/*
This function returns a normalized image with values between 0 and 1.
HACK: parameters are hard coded but some of them should be passed by flags.
*/
template <typename ImageType>
typename ImageType::Pointer
NormalizeImage(typename ImageType::Pointer inputImage)
{
  typedef itk::Statistics::ImageToHistogramFilter<ImageType>   HistogramFilterType;
  typedef typename HistogramFilterType::InputBooleanObjectType InputBooleanObjectType;
  typedef typename HistogramFilterType::HistogramSizeType      HistogramSizeType;
  typedef typename HistogramFilterType::HistogramType          HistogramType;

  HistogramSizeType histogramSize( 1 );
  histogramSize[0] = 256;

  typename InputBooleanObjectType::Pointer autoMinMaxInputObject = InputBooleanObjectType::New();
  autoMinMaxInputObject->Set( true );

  typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  histogramFilter->SetInput( inputImage );
  histogramFilter->SetAutoMinimumMaximumInput( autoMinMaxInputObject );
  histogramFilter->SetHistogramSize( histogramSize );
  histogramFilter->SetMarginalScale( 10.0 );
  histogramFilter->Update();

  float lowerValue = histogramFilter->GetOutput()->Quantile( 0, 0 );
  float upperValue = histogramFilter->GetOutput()->Quantile( 0, 1 );

  typedef itk::IntensityWindowingImageFilter<ImageType, ImageType> IntensityWindowingImageFilterType;
  typename IntensityWindowingImageFilterType::Pointer windowingFilter = IntensityWindowingImageFilterType::New();
  windowingFilter->SetInput( inputImage );
  windowingFilter->SetWindowMinimum( lowerValue );
  windowingFilter->SetWindowMaximum( upperValue );
  windowingFilter->SetOutputMinimum( 0 );
  windowingFilter->SetOutputMaximum( 1 );
  windowingFilter->Update();

  typename ImageType::Pointer outputImage = NULL;
  outputImage = windowingFilter->GetOutput();
  outputImage->Update();
  outputImage->DisconnectPipeline();

  return outputImage;
}

void
BRAINSFitHelper::Update(void)
{
  // Do remove intensity outliers if requested
  if(  m_RemoveIntensityOutliers > vcl_numeric_limits<float>::epsilon() )
    {
    this->m_FixedVolume = ClampNoisyTailsOfImage<FixedImageType, FixedBinaryVolumeType>(
        m_RemoveIntensityOutliers, this->m_FixedVolume.GetPointer(), this->m_FixedBinaryVolume.GetPointer() );
    this->m_PreprocessedMovingVolume = ClampNoisyTailsOfImage<MovingImageType, MovingBinaryVolumeType>(
        m_RemoveIntensityOutliers, this->m_MovingVolume.GetPointer(), this->m_MovingBinaryVolume.GetPointer() );
      {
      if( this->m_DebugLevel > 9 )
        {
          {
          typedef itk::ImageFileWriter<FixedImageType> WriterType;
          WriterType::Pointer writer = WriterType::New();
          writer->UseCompressionOn();
          writer->SetFileName("DEBUGNormalizedFixedVolume.nii.gz");
          writer->SetInput(this->m_FixedVolume);
          try
            {
            writer->Update();
            }
          catch( itk::ExceptionObject & err )
            {
            std::cout << "Exception Object caught: " << std::endl;
            std::cout << err << std::endl;
            throw;
            }
          }
          {
          typedef itk::ImageFileWriter<MovingImageType> WriterType;
          WriterType::Pointer writer = WriterType::New();
          writer->UseCompressionOn();
          writer->SetFileName("DEBUGNormalizedMovingVolume.nii.gz");
          writer->SetInput(this->m_PreprocessedMovingVolume);
          try
            {
            writer->Update();
            }
          catch( itk::ExceptionObject & err )
            {
            std::cout << "Exception Object caught: " << std::endl;
            std::cout << err << std::endl;
            throw;
            }
          }
        }
      }
    }
  else
    {
    this->m_PreprocessedMovingVolume = this->m_MovingVolume;
    }

  // Do Histogram equalization on moving image if requested.
  if( m_HistogramMatch )
    {
    typedef itk::OtsuHistogramMatchingImageFilter<FixedImageType, MovingImageType> HistogramMatchingFilterType;
    HistogramMatchingFilterType::Pointer histogramfilter = HistogramMatchingFilterType::New();

    // TODO:  Regina:  Write various histogram matching specializations and
    // compare them.
    // histogramfilter->SetForegroundMode("Otsu"); // A simple Otsu threshold
    // for each image .... BUT BE CAREFUL, need to do some quantile checking for
    // degenerate images
    // histogramfilter->SetForegroundMode("Simple"); // A simple average value
    // of the image should be used
    // histogramfilter->SetForegroundMode("Quantile"); // Only values between
    // the 25th and 66th quantile should be used.
    // histogramfilter->SetForegroundMode("Masks"); // Foreground is
    // specifically defined by masks.

    histogramfilter->SetReferenceImage(this->m_FixedVolume);
    if( this->m_FixedBinaryVolume.IsNull() )
      {
      itkGenericExceptionMacro(<< "ERROR:  Histogram matching requires a fixed mask.");
      }
    histogramfilter->SetReferenceMask( m_FixedBinaryVolume.GetPointer() );
    histogramfilter->SetInput(this->m_PreprocessedMovingVolume);
    if( this->m_MovingBinaryVolume.IsNull() )
      {
      itkGenericExceptionMacro(<< "ERROR:  Histogram matching requires a moving mask.");
      }
    histogramfilter->SetSourceMask( m_MovingBinaryVolume.GetPointer() );
    histogramfilter->SetNumberOfHistogramLevels(this->m_NumberOfHistogramBins);
    histogramfilter->SetNumberOfMatchPoints(this->m_NumberOfMatchPoints);
    histogramfilter->Update();
    this->m_PreprocessedMovingVolume = histogramfilter->GetOutput();
    if( this->m_DebugLevel > 5 )
      {
      typedef itk::ImageFileWriter<MovingImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();
      writer->SetFileName("DEBUGHISTOGRAMMATCHEDMOVING.nii.gz");
      writer->SetInput(this->m_PreprocessedMovingVolume);
      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "Exception Object caught: " << std::endl;
        std::cout << err << std::endl;
        throw;
        }
      }
    }
  else
    {
    this->m_PreprocessedMovingVolume = this->m_MovingVolume;
    }

  if( m_NormalizeInputImages )
    {
    this->m_FixedVolume = NormalizeImage< FixedImageType >( this->m_FixedVolume );
    this->m_PreprocessedMovingVolume = NormalizeImage< MovingImageType >( this->m_PreprocessedMovingVolume );
    }

  const bool     gradientfilter = false;

  GenericMetricType::Pointer metric;
  if( this->m_CostMetric == "MMI" )
    {
    typedef COMMON_MMI_METRIC_TYPE<FixedImageType, MovingImageType, FixedImageType, double> MIMetricType;
    MIMetricType::Pointer mutualInformationMetric = MIMetricType::New();
    mutualInformationMetric = mutualInformationMetric;
    mutualInformationMetric->SetNumberOfHistogramBins( this->m_NumberOfHistogramBins );
    mutualInformationMetric->SetUseMovingImageGradientFilter( gradientfilter );
    mutualInformationMetric->SetUseFixedImageGradientFilter( gradientfilter );
    mutualInformationMetric->SetUseFixedSampledPointSet( false );
    metric = mutualInformationMetric;

    this->SetupRegistration< MIMetricType >(metric);
    this->RunRegistration< MIMetricType >();
    }
  else if( this->m_CostMetric == "MSE" )
    {
    typedef itk::MeanSquaresImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, double> MSEMetricType;
    MSEMetricType::Pointer meanSquareMetric = MSEMetricType::New();
    meanSquareMetric = meanSquareMetric;
    metric = meanSquareMetric;

    this->SetupRegistration< MSEMetricType >(metric);
    this->RunRegistration< MSEMetricType >();
    }
  else if( this->m_CostMetric == "NC" )
    {
    typedef itk::CorrelationImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, double> corrMetricType;
    corrMetricType::Pointer corrMetric = corrMetricType::New();
    metric = corrMetric;

    this->SetupRegistration< corrMetricType >(metric);
    this->RunRegistration< corrMetricType >();
    }
  else if( this->m_CostMetric == "MIH" )
    {
    typedef itk::JointHistogramMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType, FixedImageType, double> MutualInformationMetricType;
    MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
    mutualInformationMetric = mutualInformationMetric;
    mutualInformationMetric->SetNumberOfHistogramBins( this->m_NumberOfHistogramBins );
    mutualInformationMetric->SetUseMovingImageGradientFilter( gradientfilter );
    mutualInformationMetric->SetUseFixedImageGradientFilter( gradientfilter );
    mutualInformationMetric->SetUseFixedSampledPointSet( false );
    mutualInformationMetric->SetVarianceForJointPDFSmoothing( 1.0 );
    metric = mutualInformationMetric;

    this->SetupRegistration< MutualInformationMetricType >(metric);
    this->RunRegistration< MutualInformationMetricType >();
    }
  else
    {
    std::cout << "Metric \"" << this->m_CostMetric << "\" not valid!" << std::endl;
    }
}

void
BRAINSFitHelper::PrintSelf(std::ostream & os, Indent indent) const
{
  // Superclass::PrintSelf(os,indent);
  os << indent << "FixedVolume:\n"  <<   this->m_FixedVolume << std::endl;
  os << indent << "MovingVolume:\n" <<   this->m_MovingVolume << std::endl;
  os << indent << "PreprocessedMovingVolume:\n" <<   this->m_PreprocessedMovingVolume << std::endl;
  if( this->m_FixedBinaryVolume.IsNotNull() )
    {
    os << indent << "FixedBinaryVolume:\n" << this->m_FixedBinaryVolume << std::endl;
    }
  else
    {
    os << indent << "FixedBinaryVolume: IS NULL" << std::endl;
    }
  if( this->m_MovingBinaryVolume.IsNotNull() )
    {
    os << indent << "MovingBinaryVolume:\n" << this->m_MovingBinaryVolume << std::endl;
    }
  else
    {
    os << indent << "MovingBinaryVolume: IS NULL" << std::endl;
    }
  os << indent << "SamplingPercentage:      " << this->m_SamplingPercentage << std::endl;

  os << indent << "NumberOfIterations:    [";
  for( unsigned int q = 0; q < this->m_NumberOfIterations.size(); ++q )
    {
    os << this->m_NumberOfIterations[q] << " ";
    }
  os << "]" << std::endl;
  os << indent << "NumberOfHistogramBins:" << this->m_NumberOfHistogramBins << std::endl;
  os << indent << "MaximumStepLength:    " << this->m_MaximumStepLength << std::endl;
  os << indent << "MinimumStepLength:     [";
  for( unsigned int q = 0; q < this->m_MinimumStepLength.size(); ++q )
  {
      os << this->m_MinimumStepLength[q] << " ";
  }
  os << "]" << std::endl;
  os << indent << "TransformType:     [";
  for( unsigned int q = 0; q < this->m_TransformType.size(); ++q )
    {
    os << this->m_TransformType[q] << " ";
    }
  os << "]" << std::endl;

  os << indent << "RelaxationFactor:    " << this->m_RelaxationFactor << std::endl;
  os << indent << "TranslationScale:    " << this->m_TranslationScale << std::endl;
  os << indent << "ReproportionScale:   " << this->m_ReproportionScale << std::endl;
  os << indent << "SkewScale:           " << this->m_SkewScale << std::endl;
  os << indent << "UseCachingOfBSplineWeightsMode: " << this->m_UseCachingOfBSplineWeightsMode << std::endl;
  os << indent << "BackgroundFillValue:            " << this->m_BackgroundFillValue << std::endl;
  os << indent << "InitializeTransformMode:        " << this->m_InitializeTransformMode << std::endl;
  os << indent << "MaskInferiorCutOffFromCenter:   " << this->m_MaskInferiorCutOffFromCenter << std::endl;
  os << indent << "ActualNumberOfIterations:       " << this->m_ActualNumberOfIterations << std::endl;
  os << indent << "PermittedNumberOfIterations:       " << this->m_PermittedNumberOfIterations << std::endl;

  os << indent << "SplineGridSize:     [";
  for( unsigned int q = 0; q < this->m_SplineGridSize.size(); ++q )
    {
    os << this->m_SplineGridSize[q] << " ";
    }
  os << "]" << std::endl;

  os << indent << "PermitParameterVariation:     [";
  for( unsigned int q = 0; q < this->m_PermitParameterVariation.size(); ++q )
    {
    os << this->m_PermitParameterVariation[q] << " ";
    }
  os << "]" << std::endl;

  if( m_CurrentGenericTransform.IsNotNull() )
    {
    os << indent << "CurrentGenericTransform:\n" << this->m_CurrentGenericTransform << std::endl;
    }
  else
    {
    os << indent << "CurrentGenericTransform: IS NULL" << std::endl;
    }
  os << indent << "CostMetric:       " << this->m_CostMetric << std::endl;
}

void
BRAINSFitHelper::PrintCommandLine(const bool dumpTempVolumes, const std::string & suffix) const
{
  std::cout << "The equivalent command line to the current run would be:" << std::endl;

  const std::string fixedVolumeString("DEBUGFixedVolume_" + suffix + ".nii.gz");
  const std::string movingVolumeString("DEBUGMovingVolume_" + suffix + ".nii.gz");
  const std::string fixedBinaryVolumeString("DEBUGFixedBinaryVolume_" + suffix + ".nii.gz");
  const std::string movingBinaryVolumeString("DEBUGMovingBinaryVolume_" + suffix + ".nii.gz");

  std::ostringstream oss;

  oss << "BRAINSFit \\" << std::endl;
  if( dumpTempVolumes == true )
    {
      {
      typedef itk::ImageFileWriter<FixedImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();
      writer->SetFileName(fixedVolumeString);
      writer->SetInput(this->m_FixedVolume);
      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        oss << "Exception Object caught: " << std::endl;
        oss << err << std::endl;
        throw;
        }
      }
      {
      typedef itk::ImageFileWriter<MovingImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();
      writer->SetFileName(movingVolumeString);
      writer->SetInput(this->m_MovingVolume);
      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        oss << "Exception Object caught: " << std::endl;
        oss << err << std::endl;
        throw;
        }
      }
    }
  oss << "--costMetric " << this->m_CostMetric << " \\" << std::endl;
  oss << "--fixedVolume "  <<  fixedVolumeString   << "  \\" << std::endl;
  oss << "--movingVolume " <<  movingVolumeString  << "  \\" << std::endl;
  if( this->m_HistogramMatch )
    {
    oss << "--histogramMatch " <<  "  \\" << std::endl;
    }

    {
    if( this->m_FixedBinaryVolume.IsNotNull() )
      {
      oss << "--fixedBinaryVolume " << fixedBinaryVolumeString  << "  \\" << std::endl;
        {
          {
          const MaskImageType::ConstPointer tempOutputFixedVolumeROI =
            ExtractConstPointerToImageMaskFromImageSpatialObject(m_FixedBinaryVolume.GetPointer() );
          itkUtil::WriteConstImage<MaskImageType>(tempOutputFixedVolumeROI.GetPointer(), fixedBinaryVolumeString);
          }
        }
      }
    if( this->m_MovingBinaryVolume.IsNotNull() )
      {
      oss << "--movingBinaryVolume " << movingBinaryVolumeString  << "  \\" << std::endl;
        {
          {
          const MaskImageType::ConstPointer tempOutputMovingVolumeROI =
            ExtractConstPointerToImageMaskFromImageSpatialObject(m_MovingBinaryVolume.GetPointer() );
          itkUtil::WriteConstImage<MaskImageType>(tempOutputMovingVolumeROI.GetPointer(), movingBinaryVolumeString);
          }
        }
      }
    if( this->m_FixedBinaryVolume.IsNotNull()  || this->m_MovingBinaryVolume.IsNotNull() )
      {
      oss << "--maskProcessingMode ROI "   << "  \\" << std::endl;
      }
    }
  oss << "--samplingPercentage " << this->m_SamplingPercentage  << "  \\" << std::endl;

  oss << "--numberOfIterations ";
  for( unsigned int q = 0; q < this->m_NumberOfIterations.size(); ++q )
    {
    oss << this->m_NumberOfIterations[q];
    if( q < this->m_NumberOfIterations.size() - 1 )
      {
      oss << ",";
      }
    }
  oss << " \\" << std::endl;
  oss << "--numberOfHistogramBins " << this->m_NumberOfHistogramBins  << "  \\" << std::endl;
  oss << "--maximumStepLength " << this->m_MaximumStepLength  << "  \\" << std::endl;
  oss << "--minimumStepLength ";
  for( unsigned int q = 0; q < this->m_MinimumStepLength.size(); ++q )
  {
    oss << this->m_MinimumStepLength[q];
    if( q < this->m_MinimumStepLength.size() - 1 )
    {
      oss << ",";
    }
  }
  oss << " \\" << std::endl;
  oss << "--transformType ";
  for( unsigned int q = 0; q < this->m_TransformType.size(); ++q )
    {
    oss << this->m_TransformType[q];
    if( q < this->m_TransformType.size() - 1 )
      {
      oss << ",";
      }
    }
  oss << " \\" << std::endl;

  oss << "--relaxationFactor " << this->m_RelaxationFactor  << "  \\" << std::endl;
  oss << "--translationScale " << this->m_TranslationScale  << "  \\" << std::endl;
  oss << "--reproportionScale " << this->m_ReproportionScale  << "  \\" << std::endl;
  oss << "--skewScale " << this->m_SkewScale  << "  \\" << std::endl;
  oss << "--useCachingOfBSplineWeightsMode " << this->m_UseCachingOfBSplineWeightsMode  << "  \\" << std::endl;
  oss << "--maxBSplineDisplacement " << this->m_MaxBSplineDisplacement << " \\" << std::endl;
  oss << "--projectedGradientTolerance " << this->m_ProjectedGradientTolerance << " \\" << std::endl;
  oss << "--costFunctionConvergenceFactor " << this->m_CostFunctionConvergenceFactor << " \\" << std::endl;
  oss << "--backgroundFillValue " << this->m_BackgroundFillValue  << "  \\" << std::endl;
  oss << "--initializeTransformMode " << this->m_InitializeTransformMode  << "  \\" << std::endl;
  oss << "--maskInferiorCutOffFromCenter " << this->m_MaskInferiorCutOffFromCenter  << "  \\" << std::endl;
  oss << "--splineGridSize ";
  for( unsigned int q = 0; q < this->m_SplineGridSize.size(); ++q )
    {
    oss << this->m_SplineGridSize[q];
    if( q < this->m_SplineGridSize.size() - 1 )
      {
      oss << ",";
      }
    }
  oss << " \\" << std::endl;

  if( !this->m_PermitParameterVariation.empty() )
    {
    oss << "--permitParameterVariation ";
    for( unsigned int q = 0; q < this->m_PermitParameterVariation.size(); ++q )
      {
      oss << this->m_PermitParameterVariation[q];
      if( q < this->m_PermitParameterVariation.size() - 1 )
        {
        oss << ",";
        }
      }
    oss << " \\" << std::endl;
    }
  if( m_CurrentGenericTransform.IsNotNull() )
    {
    const std::string initialTransformString("DEBUGInitialTransform_" + suffix + ".h5");
    WriteBothTransformsToDisk(this->m_CurrentGenericTransform.GetPointer(), initialTransformString, "");
    oss << "--initialTransform " << initialTransformString  << "  \\" << std::endl;
    }
    {
    const std::string outputVolume("DEBUGOutputVolume_" + suffix + ".nii.gz");
    oss << "--outputVolume " << outputVolume  << "  \\" << std::endl;
    std::cout << oss.str() << std::endl;
    }
    {
    const std::string outputTransform("DEBUGOutputTransform" + suffix + ".h5");
    oss << "--outputTransform " << outputTransform  << "  \\" << std::endl;
    std::cout << oss.str() << std::endl;
    }
  oss << "--useROIBSpline " << this->m_UseROIBSpline << " \\" << std::endl;
  const std::string TesterScript("DEBUGScript" + suffix + ".sh");
  std::ofstream     myScript;
  myScript.open( TesterScript.c_str() );
  myScript << oss.str() << std::endl;
  myScript.close();
}

void
BRAINSFitHelper::GenerateData()
{
  this->Update();
}
} // end namespace itk
