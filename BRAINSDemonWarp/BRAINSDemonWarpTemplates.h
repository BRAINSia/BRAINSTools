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
#ifndef __BRAINSDemonWarpTemplates_h
#define __BRAINSDemonWarpTemplates_h

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include "itkImage.h"
#include "itkIndex.h"
#include "itkSize.h"
#include "itkExceptionObject.h"
#include "itkPDEDeformableRegistrationFilter.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"
#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkDiffeomorphicDemonsRegistrationWithMaskFilter.h"
#include "itkVectorDiffeomorphicDemonsRegistrationFilter.h"
#include "itkESMDemonsRegistrationWithMaskFunction.h"
#include "itkArray.h"

#include "BRAINSCommonLib.h"
#include "BRAINSTypes.h"

#include "BRAINSDemonWarpCommonLibWin32Header.h"
#include "GenericTransformImage.h"
#include "VBRAINSDemonWarp.h"
#include "BRAINSDemonWarp.h"

#include "itkSymmetricForcesDemonsRegistrationFunction.h"
#include "itkFastSymmetricForcesDemonsRegistrationFunction.h"

#include "itkIO.h"
#if defined( USE_DebugImageViewer )
#  include "DebugImageViewerClient.h"
extern DebugImageViewerClient DebugImageDisplaySender;
#endif

#include "VcommandIterationupdate.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"

#include "commandIterationupdate.h"
#include "ReadMask.h"
#include "itkImageMaskSpatialObject.h"

extern void
PrintDataTypeStrings( void );

extern int
CompareNoCase( const std::string & s, const std::string & s2 );

extern int
BRAINSResample( int argc, char * argv[] );

extern void
ProcessOutputType_uchar( struct BRAINSDemonWarpAppParameters & command );

extern void
ProcessOutputType_short( struct BRAINSDemonWarpAppParameters & command );

extern void
ProcessOutputType_ushort( struct BRAINSDemonWarpAppParameters & command );

extern void
ProcessOutputType_int( struct BRAINSDemonWarpAppParameters & command );

extern void
ProcessOutputType_uint( struct BRAINSDemonWarpAppParameters & command );

extern void
ProcessOutputType_float( struct BRAINSDemonWarpAppParameters & command );

extern void
ProcessOutputType_double( struct BRAINSDemonWarpAppParameters & command );

extern void
VectorProcessOutputType_uchar( struct BRAINSDemonWarpAppParameters & command );

extern void
VectorProcessOutputType_short( struct BRAINSDemonWarpAppParameters & command );

extern void
VectorProcessOutputType_ushort( struct BRAINSDemonWarpAppParameters & command );

extern void
VectorProcessOutputType_int( struct BRAINSDemonWarpAppParameters & command );

extern void
VectorProcessOutputType_uint( struct BRAINSDemonWarpAppParameters & command );

extern void
VectorProcessOutputType_float( struct BRAINSDemonWarpAppParameters & command );

extern void
VectorProcessOutputType_double( struct BRAINSDemonWarpAppParameters & command );

struct BRAINSDemonWarpAppParameters
{
  std::string     movingVolume;
  std::string     fixedVolume;
  std::string     outputVolume;
  std::string     outputDisplacementFieldVolume;
  std::string     inputPixelType;
  std::string     outputPixelType;
  std::string     outputDisplacementFieldPrefix;
  std::string     outputCheckerboardVolume;
  std::string     registrationFilterType;
  itk::Index< 3 > checkerboardPatternSubdivisions;
  bool            outputNormalized;
  bool            outputDebug;
  std::string     maskProcessingMode;
  std::string     fixedBinaryVolume;
  std::string     movingBinaryVolume;
  int             lowerThresholdForBOBF;
  int             upperThresholdForBOBF;
  int             backgroundFillValue;
  itk::Index< 3 > seedForBOBF;
  itk::Index< 3 > neighborhoodForBOBF;
  itk::Size< 3 >  medianFilterSize;
  /*/Not yet implemented
    bool forceCoronalZeroOrigin;
    std::string movingLandmarks;
    std::string fixedLandmarks;
    std::string initializeWithFourier;
    */
  std::string  initializeWithDisplacementField;
  std::string  initializeWithTransform;
  unsigned int numberOfBCHApproximationTerms;

  /** Smoothing sigma for the deformation field at each iteration.*/
  float smoothDisplacementFieldSigma;

  /** Maximum lengthof an update vector. */
  float maxStepLength;

  /** Type of gradient used for computing the demons force. */
  unsigned int gradientType;

  /** Smoothing sigma for the update field at each iteration. */
  float smoothingUp;

  /** Intensity_histogram_matching. */
  bool histogramMatch;

  /** ShrinkFactors type. */
  using ShrinkFactorsType = itk::FixedArray< unsigned int, 3 >;

  /** IterationArray type. */
  using IterationsArrayType = itk::Array< unsigned int >;
  unsigned long       numberOfHistogramLevels;
  unsigned long       numberOfMatchPoints;
  unsigned short      numberOfLevels;
  ShrinkFactorsType   theMovingImageShrinkFactors;
  ShrinkFactorsType   theFixedImageShrinkFactors;
  IterationsArrayType numberOfIterations;
  // VECTORPARAMS
  std::vector< std::string > vectorMovingVolume;
  std::vector< std::string > vectorFixedVolume;
  bool                       makeBOBF;
  using WeightFactorsType = itk::Array< float >;
  WeightFactorsType weightFactors;
  std::string       interpolationMode;
};

// This function calls the Thirion registration filter setting all the
// parameters.
template < typename InPixelType, typename OutPixelType >
void
ThirionFunction( const struct BRAINSDemonWarpAppParameters & command )
{
  constexpr int dims = 3;

  using ImageType = itk::Image< InPixelType, dims >;
  using TRealImage = itk::Image< float, dims >;
  using OutputImageType = itk::Image< OutPixelType, dims >;
  using TDisplacementField = itk::Image< itk::Vector< float, dims >, dims >;

  using MaskPixelType = unsigned char;
  using MaskImageType = itk::Image< MaskPixelType, dims >;
  using CastImageFilter = itk::CastImageFilter< TRealImage, MaskImageType >;

  using ImageMaskSpatialObjectType = itk::ImageMaskSpatialObject< dims >;

  //
  // If optional landmark files given, will use landmark registration to
  // generate
  // a deformation field to prime the thirion demons registration.

  typedef typename itk::BRAINSDemonWarp< ImageType, TRealImage, OutputImageType > AppType;
  typename AppType::Pointer                                                       app = AppType::New();

  // Set up the diffeomorphic demons filter with mask

  if ( command.outputDebug )
  {
    std::cout << command.registrationFilterType << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }

  typename CommandIterationUpdate< float, 3 >::Pointer observer;
  if ( command.outputDebug )
  {
    observer = CommandIterationUpdate< float, 3 >::New();
  }
  {
    // Set up the demons filter
    using BaseRegistrationFilterType =
      typename itk::PDEDeformableRegistrationFilter< TRealImage, TRealImage, TDisplacementField >;
    // BaseRegistrationFilterType::Pointer filter =
    //   BaseRegistrationFilterType::New();
    typename BaseRegistrationFilterType::Pointer filter;

    if ( command.registrationFilterType == "Demons" )
    {
      using ActualRegistrationFilterType =
        typename itk::DemonsRegistrationFilter< TRealImage, TRealImage, TDisplacementField >;
      ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
      // INFO:  Review this value setting with Insight Journal Diffeomorphic
      // implementation.
      // actualfilter->SetMaximumUpdateStepLength( command.maxStepLength );
      // NOTE: GRADIENT TYPE NOT AVAILABLE IN DemonsRegistrationFilter
      // using GradientType = ActualRegistrationFilterType::GradientType;
      // actualfilter->SetUseGradientType(
      // static_cast<GradientType>(command.gradientType) );
      // actualfilter->SetUseMovingImageGradient(true);
      filter = actualfilter;
    }
    else if ( command.registrationFilterType == "Diffeomorphic" )
    {
      using ActualRegistrationFilterType =
        typename itk::DiffeomorphicDemonsRegistrationWithMaskFilter< TRealImage, TRealImage, TDisplacementField >;
      typename ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();

      using GradientType = typename ActualRegistrationFilterType::GradientType;
      actualfilter->SetMaximumUpdateStepLength( command.maxStepLength );
      actualfilter->SetUseGradientType( static_cast< GradientType >( command.gradientType ) );
      // It would be preferable that this would be part of the "Application"
      // INFO:  Move this bit of data into the application portion.
      if ( command.maskProcessingMode == "ROIAUTO" )
      {
        if ( ( command.fixedBinaryVolume != "" ) || ( command.movingBinaryVolume != "" ) )
        {
          itkGenericExceptionMacro( << "ERROR:  Can not specify mask file names when ROIAUTO "
                                    << "is used for the maskProcessingMode" )
        }
        std::cout << "Diffeomorphic with autogenerated Mask!!!!!!!" << std::endl;
        typename TRealImage::Pointer movingBinaryVolumeImage;
        typename TRealImage::Pointer fixedBinaryVolumeImage;
        constexpr double             otsuPercentileThreshold = 0.01;
        constexpr int                closingSize = 7;
        // using LargeIntegerImage = itk::Image<signed long, dims>;

        typename TRealImage::Pointer fixedVolume = itkUtil::ReadImage< TRealImage >( command.fixedVolume.c_str() );
        //       fixedBinaryVolumeImage =
        // FindLargestForgroundFilledMask<TRealImage>(
        //       fixedVolume,
        //       otsuPercentileThreshold,
        //       closingSize);
        using LFFMaskFilterType = itk::LargestForegroundFilledMaskImageFilter< TRealImage >;
        LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
        LFF->SetInput( fixedVolume );
        LFF->SetOtsuPercentileThreshold( otsuPercentileThreshold );
        LFF->SetClosingSize( closingSize );
        LFF->Update();
        fixedBinaryVolumeImage = LFF->GetOutput();

        typename CastImageFilter::Pointer castFixedMaskImage = CastImageFilter::New();
        castFixedMaskImage->SetInput( fixedBinaryVolumeImage );
        castFixedMaskImage->Update();

        typename MaskImageType::Pointer fm = castFixedMaskImage->GetOutput();
        DebugOutput( MaskImageType, fm );

        // convert mask image to mask
        typename ImageMaskSpatialObjectType::Pointer fixedMask = ImageMaskSpatialObjectType::New();
        fixedMask->SetImage( castFixedMaskImage->GetOutput() );
        fixedMask->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()

        typename TRealImage::Pointer movingVolume = itkUtil::ReadImage< TRealImage >( command.movingVolume.c_str() );
        LFF->SetInput( movingVolume );
        LFF->SetOtsuPercentileThreshold( otsuPercentileThreshold );
        LFF->SetClosingSize( closingSize );
        LFF->Update();
        movingBinaryVolumeImage = LFF->GetOutput();

        typename CastImageFilter::Pointer castMovingMaskImage = CastImageFilter::New();
        castMovingMaskImage->SetInput( movingBinaryVolumeImage );
        castMovingMaskImage->Update();
        typename MaskImageType::Pointer mm = castMovingMaskImage->GetOutput();
        DebugOutput( MaskImageType, mm );

        // convert mask image to mask
        typename ImageMaskSpatialObjectType::Pointer movingMask = ImageMaskSpatialObjectType::New();
        movingMask->SetImage( castMovingMaskImage->GetOutput() );
        movingMask->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()

        actualfilter->SetFixedImageMask( dynamic_cast< SpatialObjectType * >( fixedMask.GetPointer() ) );
        actualfilter->SetMovingImageMask( dynamic_cast< SpatialObjectType * >( movingMask.GetPointer() ) );
      }
      else if ( command.maskProcessingMode == "ROI" )
      {
        if ( ( command.fixedBinaryVolume == "" ) || ( command.movingBinaryVolume == "" ) )
        {
          itkGenericExceptionMacro( << "ERROR:  Must specify mask file names"
                                       " when ROI is used for the maskProcessingMode" );
        }
        std::cout << "Diffeomorphic with Mask!!!!!!!" << std::endl;
        typename TRealImage::Pointer fixedVolume = itkUtil::ReadImage< TRealImage >( command.fixedVolume.c_str() );
        typename TRealImage::Pointer movingVolume = itkUtil::ReadImage< TRealImage >( command.movingVolume.c_str() );

        SpatialObjectType::Pointer fixedMask =
          ReadImageMask< SpatialObjectType, dims >( command.fixedBinaryVolume, fixedVolume );
        SpatialObjectType::Pointer movingMask =
          ReadImageMask< SpatialObjectType, dims >( command.movingBinaryVolume, movingVolume );
        actualfilter->SetFixedImageMask( fixedMask );
        actualfilter->SetMovingImageMask( movingMask );
      }
      filter = actualfilter;
    }
    else if ( command.registrationFilterType == "FastSymmetricForces" )
    {
      // s <- s + u (ITK basic implementation)
      using ActualRegistrationFilterType =
        typename itk::FastSymmetricForcesDemonsRegistrationFilter< TRealImage, TRealImage, TDisplacementField >;
      using GradientType = typename ActualRegistrationFilterType::GradientType;
      typename ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
      // INFO:  Review this value setting.
      actualfilter->SetMaximumUpdateStepLength( command.maxStepLength );
      actualfilter->SetUseGradientType( static_cast< GradientType >( command.gradientType ) );
      filter = actualfilter;
    }
    /*
    else if(command.registrationFilterType == "UseFirstOrderExpOn Diffeomorphic Registration")
    {
    //INFO:  Review this value setting with Insight Journal Diffeomorphic implementation.
    // s <- s o (Id + u) (Diffeomorphic demons)
    // This is simply a crude diffeomorphic demons
    // where the exponential is computed in 0 iteration
    using ActualRegistrationFilterType = typename itk::DiffeomorphicDemonsRegistrationFilter  < TRealImage, TRealImage,
    TDisplacementField>; using GradientType = typename ActualRegistrationFilterType::GradientType;
    ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
    //INFO:  HACK: Make sure that MaxLength and GradientTypes are set.
    actualfilter->SetMaximumUpdateStepLength( command.maxStepLength );
    actualfilter->SetUseGradientType( static_cast<GradientType>(command.gradientType) );
    actualfilter->UseFirstOrderExpOn();
    filter = actualfilter;
    }
    */
    else
    {
      std::cerr << "Unknown Registration Filter type: " << command.registrationFilterType << std::endl;
      std::cerr.flush();
      throw;
    }

    // INFO:  Review this value setting with Insight Journal Diffeomorphic
    // implementation.
    if ( command.smoothDisplacementFieldSigma > 0.1 )
    {
      if ( command.outputDebug )
      {
        std::cout << " Smoothing is on ....." << std::endl;
      }
      filter->SmoothDisplacementFieldOn();
      filter->SetStandardDeviations( command.smoothDisplacementFieldSigma );
    }
    else
    {
      filter->SmoothDisplacementFieldOff();
    }
    if ( command.smoothingUp > 0.1 )
    {
      if ( command.outputDebug )
      {
        std::cout << " Smoothing at update....." << std::endl;
      }
      filter->SmoothUpdateFieldOn();
      filter->SetUpdateFieldStandardDeviations( command.smoothingUp );
    }
    else
    {
      filter->SmoothUpdateFieldOff();
    }
    if ( command.outputDebug )
    {
      filter->AddObserver( itk::IterationEvent(), observer );
    }

    app->SetRegistrationFilter( filter );
  }
  /*NOT YET IMPLEMENTED
    if ( command.fixedLandmarks != "none"
      && command.fixedLandmarks != ""
      && command.movingLandmarks != "none"
      && command.movingLandmarks != "" )
      {
      app->SetMovingLandmarkFilename(command.movingLandmarks);
      app->SetFixedLandmarkFilename(command.fixedLandmarks);
    if ( command.initializeWithFourier != "" )
      {
      app->SetInitialCoefficientFilename( command.initializeWithFourier.c_str () );
      }
      }
      */
  if ( command.initializeWithDisplacementField != "" )
  {
    app->SetInitialDisplacementFieldFilename( command.initializeWithDisplacementField.c_str() );
  }
  if ( command.initializeWithTransform != "" )
  {
    app->SetInitialTransformFilename( command.initializeWithTransform.c_str() );
  }

  app->SetTheMovingImageFilename( command.movingVolume.c_str() );
  app->SetTheFixedImageFilename( command.fixedVolume.c_str() );
  if ( command.outputDebug )
  {
    typename TRealImage::Pointer fixedVolume = itkUtil::ReadImage< TRealImage >( command.fixedVolume.c_str() );
    typename TRealImage::Pointer movingVolume = itkUtil::ReadImage< TRealImage >( command.movingVolume.c_str() );
    observer->SetMovingImage( movingVolume );
    observer->SetFixedImage( fixedVolume );
  }
  app->SetWarpedImageName( command.outputVolume.c_str() );
  app->SetMedianFilterSize( command.medianFilterSize );

  // Set the other optional arguments if specified by the user.
  if ( command.outputDisplacementFieldPrefix != "" )
  {
    app->SetDisplacementBaseName( command.outputDisplacementFieldPrefix.c_str() );
  }
  if ( command.outputDisplacementFieldVolume != "" )
  {
    app->SetDisplacementFieldOutputName( command.outputDisplacementFieldVolume.c_str() );
  }

  if ( command.outputCheckerboardVolume != "" )
  {
    app->SetCheckerBoardFilename( command.outputCheckerboardVolume.c_str() );
    unsigned int array[3] = { static_cast< unsigned int >( command.checkerboardPatternSubdivisions[0] ),
                              static_cast< unsigned int >( command.checkerboardPatternSubdivisions[1] ),
                              static_cast< unsigned int >( command.checkerboardPatternSubdivisions[2] ) };
    app->SetCheckerBoardPattern( array );
  }

  if ( command.outputNormalized )
  {
    std::string normalize = "ON"; // INFO:  SetOutNormalized should be a
                                  // boolean
                                  // not a string.
    app->SetOutNormalized( normalize.c_str() );
  }

  if ( command.outputDebug )
  {
    bool debug = true;
    app->SetOutDebug( debug ); // INFO:  SetOutDebug should be a boolean not a
                               // string.
  }

  app->SetTheMovingImageShrinkFactors( command.theMovingImageShrinkFactors );
  app->SetTheFixedImageShrinkFactors( command.theFixedImageShrinkFactors );

  app->SetUseHistogramMatching( command.histogramMatch );
  if ( app->GetUseHistogramMatching() )
  {
    if ( command.outputDebug )
    {
      std::cout << " Use Histogram Matching....." << std::endl;
    }
    app->SetNumberOfHistogramLevels( command.numberOfHistogramLevels );
    app->SetNumberOfMatchPoints( command.numberOfMatchPoints );
  }

  app->SetNumberOfLevels( command.numberOfLevels );
  app->SetNumberOfIterations( command.numberOfIterations );
  app->SetInterpolationMode( command.interpolationMode );

  if ( ( command.maskProcessingMode == "NOMASK" ) &&
       ( ( command.fixedBinaryVolume != "" ) || ( command.movingBinaryVolume != "" ) ) )
  {
    itkGenericExceptionMacro( << "ERROR:  Can not specify mask file names when "
                              << "the default of NOMASK is used for the maskProcessingMode" );
  }
  // If making BOBF option is specified Initialize its parameters
  if ( command.maskProcessingMode == "BOBF" )
  {
    if ( ( command.fixedBinaryVolume == "" ) || ( command.movingBinaryVolume == "" ) )
    {
      itkGenericExceptionMacro( << "Error: If BOBF option is set for maskProcessingMode"
                                << " then the fixed mask name and moving mask file "
                                << "name should be specified" );
    }

    app->SetFixedBinaryVolume( command.fixedBinaryVolume.c_str() );
    app->SetMovingBinaryVolume( command.movingBinaryVolume.c_str() );
    app->SetLower( command.lowerThresholdForBOBF );
    app->SetUpper( command.upperThresholdForBOBF );
    typename ImageType::SizeType radius;
    radius[0] = command.neighborhoodForBOBF[0]; // Radius along X
    radius[1] = command.neighborhoodForBOBF[1]; // Radius along Y
    radius[2] = command.neighborhoodForBOBF[2]; // Radius along Z
    app->SetRadius( radius );
    typename ImageType::IndexType seed;
    seed[0] = command.seedForBOBF[0]; // Seed in X dimension;
    seed[1] = command.seedForBOBF[1]; // Seed in Y dimension;
    seed[2] = command.seedForBOBF[2]; // Seed in Z dimension;
    app->SetSeed( seed );
  }
  if ( command.outputDebug )
  {
    std::cout << "Setting Default PixelValue: " << command.backgroundFillValue << "." << std::endl;
  }
  app->SetDefaultPixelValue( command.backgroundFillValue );
  if ( command.outputDebug )
  {
    std::cout << "Running Thirion Registration" << std::endl;
  }
  try
  {
    app->Execute();
  }
  catch ( itk::ExceptionObject & err )
  {
    std::cout << "Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
  }
  catch ( ... )
  {
    std::cout << "Caught a non-ITK exception " << __FILE__ << " " << __LINE__ << std::endl;
  }

  return;
}

// This function calls the Thirion registration filter setting all the
// parameters.
template < typename InPixelType, typename OutPixelType >
void
ProcessAppType( const struct BRAINSDemonWarpAppParameters & command )
{
  ThirionFunction< InPixelType, OutPixelType >( command );
}

// This function processes the output data type.
template < typename PixelType >
void
ProcessOutputType( struct BRAINSDemonWarpAppParameters & command )
{
  if ( command.outputPixelType != "" )
  {
    // process the string for the data type
    if ( CompareNoCase( command.outputPixelType, std::string( "uchar" ) ) == 0 )
    {
      ProcessAppType< PixelType, unsigned char >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "short" ) ) == 0 )
    {
      ProcessAppType< PixelType, short >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "ushort" ) ) == 0 )
    {
      ProcessAppType< PixelType, unsigned short >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "int" ) ) == 0 )
    {
      ProcessAppType< PixelType, int >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "float" ) ) == 0 )
    {
      ProcessAppType< PixelType, float >( command );
    }
#ifdef _USE_UNCOMMON_TYPES // This is commented out because it causes too many
                           // segments in one object file for the intel
                           // compiler
    else if ( CompareNoCase( command.outputPixelType, std::string( "uint" ) ) == 0 )
    {
      ProcessAppType< PixelType, unsigned int >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "double" ) ) == 0 )
    {
      ProcessAppType< PixelType, double >( command );
    }
#endif
    else
    {
      std::cout << "Error. Invalid data type for -outtype!  Use one of these:" << std::endl;
      PrintDataTypeStrings();
      throw;
    }
  }
  else
  {
    ProcessAppType< PixelType, float >( command );
  }
}

template < typename InPixelType, typename OutPixelType >
void
VectorThirionFunction( const struct BRAINSDemonWarpAppParameters & command )
{
  constexpr int dims = 3;

  using ImageType = itk::Image< InPixelType, dims >;
  using TRealImage = itk::Image< float, dims >;
  using TVectorImage = itk::VectorImage< float, dims >;
  using OutputImageType = itk::Image< OutPixelType, dims >;
  using TDisplacementField = itk::Image< itk::Vector< float, dims >, dims >;
  //
  // If optional landmark files given, will use landmark registration to
  // generate
  // a deformation field to prime the thirion demons registration.

  typedef typename itk::VBRAINSDemonWarp< ImageType, TRealImage, OutputImageType > AppType;
  typename AppType::Pointer                                                        app = AppType::New();

  // Set up the demons filter
  using BaseRegistrationFilterType =
    typename itk::PDEDeformableRegistrationFilter< TRealImage, TRealImage, TDisplacementField >;
  // BaseRegistrationFilterType::Pointer filter =
  //   BaseRegistrationFilterType::New();
  typename BaseRegistrationFilterType::Pointer filter;
  if ( command.outputDebug )
  {
    std::cout << command.registrationFilterType << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }

  if ( command.registrationFilterType == "Demons" )
  {
    if ( command.vectorMovingVolume.size() == 1 )
    {
      using ActualRegistrationFilterType =
        typename itk::DemonsRegistrationFilter< TRealImage, TRealImage, TDisplacementField >;
      ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
      filter = actualfilter;
    }
    else
    {
      std::cout << "Thirion demons does not support multi-input images!" << std::endl;
      throw;
    }
  }

  else if ( command.registrationFilterType == "Diffeomorphic" )
  {
    //    std::cout << "Use Diffeomorphic Registration" << std::endl;
    if ( command.vectorMovingVolume.size() == 1 )
    {
      using ActualRegistrationFilterType =
        typename itk::DiffeomorphicDemonsRegistrationFilter< TRealImage, TRealImage, TDisplacementField >;
      using GradientType = typename ActualRegistrationFilterType::GradientType;
      typename ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
      // INFO:  Review this value setting with Insight Journal Diffeomorphic
      // implementation.
      actualfilter->SetMaximumUpdateStepLength( command.maxStepLength );
      actualfilter->SetUseGradientType( static_cast< GradientType >( command.gradientType ) );
      filter = actualfilter;
    }
    else
    {
      using ActualRegistrationFilterType =
        typename itk::VectorDiffeomorphicDemonsRegistrationFilter< TVectorImage, TVectorImage, TDisplacementField >;
      using GradientType = typename ActualRegistrationFilterType::GradientType;
      typename ActualRegistrationFilterType::Pointer VDDfilter = ActualRegistrationFilterType::New();
      // INFO:  Review this value setting with Insight Journal Diffeomorphic
      // implementation.
      VDDfilter->SetMaximumUpdateStepLength( command.maxStepLength );
      VDDfilter->SetUseGradientType( static_cast< GradientType >( command.gradientType ) );
      if ( command.smoothDisplacementFieldSigma > 0.1 )
      {
        if ( command.outputDebug )
        {
          std::cout << " Smoothing is on ....." << std::endl;
        }
        VDDfilter->SmoothDisplacementFieldOn();
        VDDfilter->SetStandardDeviations( command.smoothDisplacementFieldSigma );
      }
      else
      {
        VDDfilter->SmoothDisplacementFieldOff();
      }
      if ( command.smoothingUp > 0.1 )
      {
        if ( command.outputDebug )
        {
          std::cout << " Smoothing at update....." << std::endl;
        }
        VDDfilter->SmoothUpdateFieldOn();
        VDDfilter->SetUpdateFieldStandardDeviations( command.smoothingUp );
      }
      else
      {
        VDDfilter->SmoothUpdateFieldOff();
      }
      if ( command.outputDebug )
      {
        typename VCommandIterationUpdate< float, 3 >::Pointer observer = VCommandIterationUpdate< float, 3 >::New();
        VDDfilter->AddObserver( itk::IterationEvent(), observer );
      }

      app->SetVectorRegistrationFilter( VDDfilter );
    }
  }
  else if ( command.registrationFilterType == "FastSymmetricForces" )
  {
    // s <- s + u (ITK basic implementation)
    if ( command.vectorMovingVolume.size() == 1 )
    {
      using ActualRegistrationFilterType =
        typename itk::FastSymmetricForcesDemonsRegistrationFilter< TRealImage, TRealImage, TDisplacementField >;
      using GradientType = typename ActualRegistrationFilterType::GradientType;
      typename ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
      // INFO:  Review this value setting.
      actualfilter->SetMaximumUpdateStepLength( command.maxStepLength );
      actualfilter->SetUseGradientType( static_cast< GradientType >( command.gradientType ) );
      filter = actualfilter;
    }
    else
    {
      std::cout << "FastSymmetricForces demons does not support multi-input images!" << std::endl;
      throw;
    }
  }
  else
  {
    std::cerr << "Unknown Registration Filter type: " << command.registrationFilterType << std::endl;
    std::cerr.flush();
    throw;
  }

  // INFO:  Review this value setting with Insight Journal Diffeomorphic
  // implementation.
  if ( command.vectorMovingVolume.size() == 1 )
  {
    if ( command.smoothDisplacementFieldSigma > 0.1 )
    {
      if ( command.outputDebug )
      {
        std::cout << " Smoothing is on ....." << std::endl;
      }
      filter->SmoothDisplacementFieldOn();
      filter->SetStandardDeviations( command.smoothDisplacementFieldSigma );
    }
    else
    {
      filter->SmoothDisplacementFieldOff();
    }
    if ( command.smoothingUp > 0.1 )
    {
      if ( command.outputDebug )
      {
        std::cout << " Smoothing at update....." << std::endl;
      }
      filter->SmoothUpdateFieldOn();
      filter->SetUpdateFieldStandardDeviations( command.smoothingUp );
    }
    else
    {
      filter->SmoothUpdateFieldOff();
    }
    if ( command.outputDebug )
    {
      typename VCommandIterationUpdate< float, 3 >::Pointer observer = VCommandIterationUpdate< float, 3 >::New();
      filter->AddObserver( itk::IterationEvent(), observer );
    }

    app->SetRegistrationFilter( filter );
  }

  /* NOT YET IMPLEMENTED
  if ( command.fixedLandmarks != "none"
    && command.fixedLandmarks != ""
    && command.movingLandmarks != "none"
    && command.movingLandmarks != "" )
    {
    app->SetMovingLandmarkFilename(command.movingLandmarks);
    app->SetFixedLandmarkFilename(command.fixedLandmarks);
    }
  if ( command.initializeWithFourier != "" )
    {
    app->SetInitialCoefficientFilename( command.initializeWithFourier.c_str () );
    }
   app->SetForceCoronalZeroOrigin (command.forceCoronalZeroOrigin);
*/
  if ( command.initializeWithDisplacementField != "" )
  {
    app->SetInitialDisplacementFieldFilename( command.initializeWithDisplacementField.c_str() );
  }
  if ( command.initializeWithTransform != "" )
  {
    app->SetInitialTransformFilename( command.initializeWithTransform.c_str() );
  }

  std::vector< std::string > fixedImageNames = command.vectorFixedVolume;
  std::vector< std::string > movingImageNames = command.vectorMovingVolume;

  app->SetTheFixedImageFilename( fixedImageNames );
  app->SetTheMovingImageFilename( movingImageNames );
  app->SetWarpedImageName( command.outputVolume.c_str() );
  app->SetInterpolationMode( command.interpolationMode );
  app->SetMedianFilterSize( command.medianFilterSize );

  // Set the other optional arguments if specified by the user.
  if ( command.outputDisplacementFieldPrefix != "" )
  {
    app->SetDisplacementBaseName( command.outputDisplacementFieldPrefix.c_str() );
  }
  if ( command.outputDisplacementFieldVolume != "" )
  {
    app->SetDisplacementFieldOutputName( command.outputDisplacementFieldVolume.c_str() );
  }

  if ( command.outputCheckerboardVolume != "" )
  {
    app->SetCheckerBoardFilename( command.outputCheckerboardVolume.c_str() );
    unsigned int array[3] = { static_cast< unsigned int >( command.checkerboardPatternSubdivisions[0] ),
                              static_cast< unsigned int >( command.checkerboardPatternSubdivisions[1] ),
                              static_cast< unsigned int >( command.checkerboardPatternSubdivisions[2] ) };
    app->SetCheckerBoardPattern( array );
  }

  if ( command.outputNormalized )
  {
    std::string normalize = "ON"; // INFO:  SetOutNormalized should be a
                                  // boolean
                                  // not a string.
    app->SetOutNormalized( normalize.c_str() );
  }

  if ( command.outputDebug )
  {
    bool debug = true;
    app->SetOutDebug( debug ); // INFO:  SetOutDebug should be a boolean not a
                               // string.
  }

  app->SetTheMovingImageShrinkFactors( command.theMovingImageShrinkFactors );
  app->SetTheFixedImageShrinkFactors( command.theFixedImageShrinkFactors );

  app->SetUseHistogramMatching( command.histogramMatch );
  if ( app->GetUseHistogramMatching() )
  {
    if ( command.outputDebug )
    {
      std::cout << " Use Histogram Matching....." << std::endl;
    }
    app->SetNumberOfHistogramLevels( command.numberOfHistogramLevels );
    app->SetNumberOfMatchPoints( command.numberOfMatchPoints );
  }

  app->SetNumberOfLevels( command.numberOfLevels );
  app->SetNumberOfIterations( command.numberOfIterations );
  app->SetInterpolationMode( command.interpolationMode );

  app->SetWeightFactors( command.weightFactors );

  // If making BOBF option is specified Initialize its parameters
  if ( command.makeBOBF )
  {
    if ( ( command.fixedBinaryVolume == "" ) || ( command.movingBinaryVolume == "" ) )
    {
      std::cout
        << "Error: If BOBF option is set then the fixed mask name and moving mask file name should be specified. \n";
      throw;
    }

    app->SetFixedBinaryVolume( command.fixedBinaryVolume.c_str() );
    app->SetMovingBinaryVolume( command.movingBinaryVolume.c_str() );
    app->SetLower( command.lowerThresholdForBOBF );
    app->SetUpper( command.upperThresholdForBOBF );
    typename ImageType::SizeType radius;
    radius[0] = command.neighborhoodForBOBF[0]; // Radius along X
    radius[1] = command.neighborhoodForBOBF[1]; // Radius along Y
    radius[2] = command.neighborhoodForBOBF[2]; // Radius along Z
    app->SetRadius( radius );
    typename ImageType::IndexType seed;
    seed[0] = command.seedForBOBF[0]; // Seed in X dimension;
    seed[1] = command.seedForBOBF[1]; // Seed in Y dimension;
    seed[2] = command.seedForBOBF[2]; // Seed in Z dimension;
    app->SetSeed( seed );
  }
  if ( command.outputDebug )
  {
    std::cout << "Setting Default PixelValue: " << command.backgroundFillValue << "." << std::endl;
  }
  app->SetDefaultPixelValue( command.backgroundFillValue );
  if ( command.outputDebug )
  {
    std::cout << "Running Thirion Registration" << std::endl;
  }
  try
  {
    app->Execute();
  }
  catch ( itk::ExceptionObject & err )
  {
    std::cout << "Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
  }
  catch ( ... )
  {
    std::cout << "Caught a non-ITK exception " << __FILE__ << " " << __LINE__ << std::endl;
  }

  return;
}

// This function calls the Thirion registration filter setting all the
// parameters.
template < typename InPixelType, typename OutPixelType >
void
VectorProcessAppType( const struct BRAINSDemonWarpAppParameters & command )
{
  VectorThirionFunction< InPixelType, OutPixelType >( command );
}

// This function processes the output data type.
template < typename PixelType >
void
VectorProcssOutputType( struct BRAINSDemonWarpAppParameters & command )
{
  if ( command.outputPixelType != "" )
  {
    // process the string for the data type
    if ( CompareNoCase( command.outputPixelType, std::string( "uchar" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, unsigned char >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "short" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, short >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "ushort" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, unsigned short >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "int" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, int >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "float" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, float >( command );
    }
#ifdef _USE_UNCOMMON_TYPES // This is commented out because it causes too many
                           // segments in one object file for the intel
                           // compiler
    else if ( CompareNoCase( command.outputPixelType, std::string( "uint" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, unsigned int >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "double" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, double >( command );
    }
#endif
    else
    {
      std::cout << "Error. Invalid data type for -outtype!  Use one of these:" << std::endl;
      PrintDataTypeStrings();
      throw;
    }
  }
  else
  {
    VectorProcessAppType< PixelType, float >( command );
  }
}

// This function processes the output data type.
template < typename PixelType >
void
VectorProcessOutputType( struct BRAINSDemonWarpAppParameters & command )
{
  if ( command.outputPixelType != "" )
  {
    // process the string for the data type
    if ( CompareNoCase( command.outputPixelType, std::string( "uchar" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, unsigned char >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "short" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, short >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "ushort" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, unsigned short >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "int" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, int >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "float" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, float >( command );
    }
#ifdef _USE_UNCOMMON_TYPES // This is commented out because it causes too many
                           // segments in one object file for the intel
                           // compiler
    else if ( CompareNoCase( command.outputPixelType, std::string( "uint" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, unsigned int >( command );
    }
    else if ( CompareNoCase( command.outputPixelType, std::string( "double" ) ) == 0 )
    {
      VectorProcessAppType< PixelType, double >( command );
    }
#endif
    else
    {
      std::cout << "Error. Invalid data type for -outtype!  Use one of these:" << std::endl;
      PrintDataTypeStrings();
      throw;
    }
  }
  else
  {
    VectorProcessAppType< PixelType, float >( command );
  }
}

#endif // __BRAINSDemonWarpTemplates_h
