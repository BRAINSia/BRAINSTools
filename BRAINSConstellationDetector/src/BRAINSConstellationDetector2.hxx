/*=========================================================================
 Author: Wei Lu, Ali Ghayoor, Hans Johnson
 SINAPSE
 University of Iowa, 2010, 2012, 2013
 =========================================================================*/

#include "PrepareOutputImages.h"
#include "landmarksConstellationDetector.h"
#include "BRAINSConstellationDetector2.h"
#include "itkOrthogonalize3DRotationMatrix.h"
#include "itkLandmarkBasedTransformInitializer.h"

#include "itkResampleInPlaceImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkCastImageFilter.h"
#include <BRAINSFitHelper.h>


namespace itk
{


template <class TInputImage, class TOutputImage>
BRAINSConstellationDetector2<TInputImage, TOutputImage>
::BRAINSConstellationDetector2()
  {
  /** Essential Parameters */
  // Inputs
  this->m_InputTemplateModel = "";
  this->m_MspQualityLevel = 2;
  this->m_OtsuPercentileThreshold = 0.01;
  this->m_AcLowerBound = 1000.0;
  this->m_CutOutHeadInOutputVolume = false;
  this->m_RescaleIntensities = false;
  this->m_TrimRescaledIntensities = 4.4172;
  this->m_RescaleIntensitiesOutputRange.push_back(40);
  this->m_RescaleIntensitiesOutputRange.push_back(4000);
  this->m_BackgroundFillValueString = "0";
  this->m_InterpolationMode = "Linear";
  this->m_OriginalInputImage = NULL;

  // Outputs
  this->m_Transform = "";
  this->m_VersorTransform = NULL;
  this->m_OutputImage = NULL;
  this->m_OutputResampledImage = NULL;
  this->m_OutputUntransformedClippedVolume = NULL;

  /** Advanced parameters */
  /** Manual Override */
  // Inputs
  this->m_ForceACPoint.push_back(0);
  this->m_ForcePCPoint.push_back(0);
  this->m_ForceVN4Point.push_back(0);
  this->m_ForceRPPoint.push_back(0);

  /** Model Override */
  // Inputs
  this->m_RadiusMPJ = -1.;
  this->m_RadiusAC = -1.;
  this->m_RadiusPC = -1.;
  this->m_RadiusVN4 = -1.;

  /** Debug Options */
  // Inputs
  this->m_Debug = false;
  this->m_Verbose = false;
  this->m_WritedebuggingImagesLevel = 0;

  // Outputs
  this->m_WriteBranded2DImage = "";
  this->m_ResultsDir = "./";
  }

template <class TInputImage, class TOutputImage>
void
BRAINSConstellationDetector2<TInputImage, TOutputImage>
::GenerateData()
{
  // file pointer for opening the setup file
  // /////////////////////////////////////////////////////////////////////////////////////////////
  LMC::globalverboseFlag = this->m_Verbose;
  globalImagedebugLevel = this->m_WritedebuggingImagesLevel;

  // /////////////////////////////////////////////////////////////////////////////////////////////
  short BackgroundFillValue;
  if( this->m_BackgroundFillValueString == std::string("BIGNEG") )
    {
    BackgroundFillValue = -32768;
    }
  else
    {
    BackgroundFillValue = atoi( this->m_BackgroundFillValueString.c_str() );
    }
  // /////////////////////////////////////////////////////////////////////////////////////////////
  // read information from the setup file, and initialize some variables
  landmarksConstellationModelIO myModel;
  myModel.ReadModelFile( this->m_InputTemplateModel );

  if( LMC::globalverboseFlag )
    {
    std::cout << "Using Model File: " << this->m_InputTemplateModel << std::endl;
    myModel.PrintHeaderInfo();
    }

  // Override some landmarks by user
  vnl_vector<double> templateRadius;    // in units of mm
  templateRadius.set_size(4);
  templateRadius[0] = ( this->m_RadiusMPJ <= 0 ) ? myModel.GetRadius("RP") : this->m_RadiusMPJ;

  templateRadius[1] = ( this->m_RadiusAC <= 0 ) ? myModel.GetRadius("AC") : this->m_RadiusAC;
  templateRadius[2] = ( this->m_RadiusPC <= 0 ) ? myModel.GetRadius("PC") : this->m_RadiusPC;
  templateRadius[3] = ( this->m_RadiusVN4 <= 0 ) ? myModel.GetRadius("VN4") : this->m_RadiusVN4;
  myModel.SetRadius("RP", templateRadius[0]);
  myModel.SetRadius("AC", templateRadius[1]);
  myModel.SetRadius("PC", templateRadius[2]);
  myModel.SetRadius("VN4", templateRadius[3]);

  // Wei: We will get the input image from filter input rather than an external
  // file
  SImageType::Pointer volOrig;
    {
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(this->m_OriginalInputImage);
    duplicator->Update();
    volOrig = duplicator->GetModifiableOutput();
    }

  // RPPC is a vector on the MSP that points from the RP point to the PC.
  //   RP------->PC
  // //////////////////////////////////////////////////////////////////////////

  if( this->m_RescaleIntensities == true )
    {
    itk::StatisticsImageFilter<SImageType>::Pointer stats =
      itk::StatisticsImageFilter<SImageType>::New();
    stats->SetInput(volOrig);
    stats->Update();
    SImageType::PixelType minPixel( stats->GetMinimum() );
    SImageType::PixelType maxPixel( stats->GetMaximum() );

    if( this->m_TrimRescaledIntensities > 0.0 )
      {
      // REFACTOR: a histogram would be traditional here, but seems
      // over-the-top;
      // I did this because it seemed to me if I knew mean, sigma, max and min,
      // then I know Something about extreme outliers.
      // Look at the setLowHigh function in landmarksConstellationCommon.h as a
      // possible replacement

      const double meanOrig( stats->GetMean() );
      const double sigmaOrig( stats->GetSigma() );

      // REFACTOR:  In percentiles, 0.0005 two-tailed has worked in the past.
      // It only makes sense to trim the upper bound since the lower bound would
      // most likely
      // represent a large region of air around the head.  But this is not so
      // when using a mask.
      // For one-tailed, an error of 0.001 corresponds to 3.29052 standard
      // deviations of normal.
      // For one-tailed, an error of 0.0001 corresponds to 3.8906 standard
      // deviations of normal.
      // For one-tailed, an error of 0.00001 corresponds to 4.4172 standard
      // deviations of normal.
      // Naturally, the constant should default at the command line, ...

      const double variationBound( ( maxPixel - meanOrig ) / sigmaOrig );
      const double trimBound(variationBound - this->m_TrimRescaledIntensities);
      if( trimBound > 0.0 )
        {
        maxPixel = static_cast<SImageType::PixelType>( maxPixel - trimBound * sigmaOrig );
        }
      }

    itk::IntensityWindowingImageFilter<SImageType, SImageType>::Pointer remapIntensityFilter =
      itk::IntensityWindowingImageFilter<SImageType, SImageType>::New();
    remapIntensityFilter->SetInput(volOrig);
    remapIntensityFilter->SetOutputMaximum(this->m_RescaleIntensitiesOutputRange[1]);
    remapIntensityFilter->SetOutputMinimum(this->m_RescaleIntensitiesOutputRange[0]);
    remapIntensityFilter->SetWindowMinimum(minPixel);
    remapIntensityFilter->SetWindowMaximum(maxPixel);
    remapIntensityFilter->Update();

    this->m_ImageToBeResampled = remapIntensityFilter->GetOutput();
    }
  else
    {
    this->m_ImageToBeResampled = volOrig;
    }

  landmarksConstellationDetector myDetector;
    {
    // a little abuse of the duplicator here
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( this->GetInput() );
    duplicator->Update();
    // The detector will use the output image after the Hough eye detector
    myDetector.SetVolOrig( duplicator->GetModifiableOutput() );

    // The detector also needs the original input if it has to fix a bad estimation of the MSP
    duplicator->SetInputImage( this->m_ImageToBeResampled );
    duplicator->Update();
    myDetector.SetOriginalInput( duplicator->GetModifiableOutput() );
    }

  myDetector.SetInputTemplateModel( myModel );
  myDetector.SetLlsMatrices( this->m_LlsMatrices );
  myDetector.SetLlsMeans( this->m_LlsMeans );
  myDetector.SetSearchRadii( this->m_SearchRadii );
  myDetector.SetResultsDir( this->m_ResultsDir );
  myDetector.SetTemplateRadius( myModel.GetRadii() );
  myDetector.SetMSPQualityLevel( this->m_MspQualityLevel );
  myDetector.SetCenterOfHeadMass( this->m_CenterOfHeadMass );
  myDetector.SetHoughEyeFailure( this->m_HoughEyeFailure );

  LandmarksMapType LandmarksEMSP = this->GetLandmarksEMSP();
  if( !LandmarksEMSP.empty() )
    {
    myDetector.SetLandmarksEMSP( this->GetLandmarksEMSP() );
    }

  if( ( this->m_landmarksEMSP.find("LE") == this->m_landmarksEMSP.end() )
      || ( this->m_landmarksEMSP.find("RE") == this->m_landmarksEMSP.end() ) )
    {
    myDetector.SetHoughEyeTransform(this->m_HoughEyeTransform);
    myDetector.SetLEPoint(this->m_LEPoint);
    myDetector.SetREPoint(this->m_REPoint);
    }

    { /** Force setting the landmark points from the command line. */
    if( this->m_ForceACPoint.size() == 3 )
      {
      SImageType::PointType manualACPoint;
      for( int i = 0; i < 3; i++ )
        {
        manualACPoint[i] = this->m_ForceACPoint[i];
        }
      myDetector.SetOriginalSpaceNamedPoint(std::string("AC"), manualACPoint);
      }
    if( this->m_ForcePCPoint.size() == 3 )
      {
      SImageType::PointType manualPCPoint;
      for( int i = 0; i < 3; i++ )
        {
        manualPCPoint[i] = this->m_ForcePCPoint[i];
        }
      myDetector.SetOriginalSpaceNamedPoint(std::string("PC"), manualPCPoint);
      }
    if( this->m_ForceVN4Point.size() == 3 )
      {
      SImageType::PointType manualVN4Point;
      for( int i = 0; i < 3; i++ )
        {
        manualVN4Point[i] = this->m_ForceVN4Point[i];
        }
      myDetector.SetOriginalSpaceNamedPoint(std::string("VN4"), manualVN4Point);
      }
    if( this->m_ForceRPPoint.size() == 3 )
      {
      SImageType::PointType manualRPPoint;
      for( int i = 0; i < 3; i++ )
        {
        manualRPPoint[i] = this->m_ForceRPPoint[i];
        }
      myDetector.SetOriginalSpaceNamedPoint(std::string("RP"), manualRPPoint);
      }
    }

  myDetector.Compute();

    {
    RigidTransformType::Pointer ZeroCenteredTransform =
      myDetector.GetACPCAlignedZeroCenteredTransform();

    this->m_VersorTransform = VersorTransformType::New();
    this->m_VersorTransform->SetFixedParameters( ZeroCenteredTransform->GetFixedParameters() );
    itk::Versor<double>               versorRotation;
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized = itk::Orthogonalize3DRotationMatrix(
      ZeroCenteredTransform->GetMatrix() );
    versorRotation.Set( CleanedOrthogonalized );
    this->m_VersorTransform->SetRotation(versorRotation);
    this->m_VersorTransform->SetTranslation( ZeroCenteredTransform->GetTranslation() );
    }

  WriteTransformToDisk( this->m_VersorTransform.GetPointer(), std::string("/tmp/OrigToACPC_VersorTransformLINE_291.mat") );

  ////////////////////////////
  // START BRAINSFit alternative
    if( ! this->m_atlasVolume.empty() )
      {
      typedef itk::ImageFileReader<SImageType> AtlasReaderType;
      AtlasReaderType::Pointer atlasReader = AtlasReaderType::New();
      atlasReader->SetFileName( this->m_atlasVolume );
      try
        {
        atlasReader->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Error while reading atlasVolume file:\n "
          << err << std::endl;
        }

      std::cout << "read atlas" << std::endl;
      // TODO: prob needs a try-catch
      LandmarksMapType referenceAtlasLandmarks = ReadSlicer3toITKLmk( this->m_atlasLandmarks );
      std::cout << "read atlas landmarks " << std::endl;
      LandmarksMapType acpcLandmarks;
      itk::PrepareOutputLandmarks(
        this->m_VersorTransform.GetPointer(), //Input RO
        myDetector.GetNamedPoints(),
        acpcLandmarks
      );

      // Create a better version of this->m_VersorTransform using BRAINSFit.
      // take the the subjects landmarks in original space, and  landmarks from a reference Atlas, and compute an initial
      // affine transform
      // ( using logic from BRAINSLandmarkInitializer) and create initToAtlasAffineTransform.

      typedef std::map<std::string, float> WeightType;
      WeightType landmarkWeights;
      if( this->m_atlasLandmarkWeights != "" )
        {
        landmarkWeights = ReadLandmarkWeights( this->m_atlasLandmarkWeights.c_str() );
        }
      typedef itk::LandmarkBasedTransformInitializer<AffineTransformType, SImageType, SImageType> LandmarkBasedInitializerType;
      typedef typename LandmarkBasedInitializerType::LandmarkPointContainer LandmarkContainerType;
      LandmarkContainerType fixedLmks;
      LandmarkContainerType movingLmks;
      typedef typename  LandmarksMapType::const_iterator LandmarkConstIterator;
      typename LandmarkBasedInitializerType::LandmarkWeightType landmarkWgts;
      for( LandmarkConstIterator fixedIt = referenceAtlasLandmarks.begin(); fixedIt != referenceAtlasLandmarks.end();
        ++fixedIt )
        {
        LandmarkConstIterator movingIt = acpcLandmarks.find( fixedIt->first );
        if( movingIt != acpcLandmarks.end() )
          {
          fixedLmks.push_back( fixedIt->second);
          movingLmks.push_back( movingIt->second);
          if( !this->m_atlasLandmarkWeights.empty() )
            {
            if( landmarkWeights.find( fixedIt->first ) != landmarkWeights.end() )
              {
              landmarkWgts.push_back( landmarkWeights[fixedIt->first] );
              }
            else
              {
              std::cout << "Landmark for " << fixedIt->first << " does not exist. "
                << "Set the weight to 0.5 "
                << std::endl;
              landmarkWgts.push_back( 0.5F );
              }
            }
          }
        else
          {
          std::cout << "i shouldnt be here" << std::endl;
          exit(-1);
          //TODO:  Throw exception
          }
        }

      typedef itk::LandmarkBasedTransformInitializer<AffineTransformType, SImageType, SImageType> LandmarkBasedInitializerType;
      typename LandmarkBasedInitializerType::Pointer landmarkBasedInitializer = LandmarkBasedInitializerType::New();

      if( !this->m_atlasLandmarkWeights.empty() )
        {
        landmarkBasedInitializer->SetLandmarkWeight( landmarkWgts );
        }
      landmarkBasedInitializer->SetFixedLandmarks( fixedLmks );
      landmarkBasedInitializer->SetMovingLandmarks( movingLmks );

      typedef itk::AffineTransform<double, Dimension> AffineTransformType;
      typename AffineTransformType::Pointer initToAtlasAffineTransform = AffineTransformType::New();
      landmarkBasedInitializer->SetTransform( initToAtlasAffineTransform );
      landmarkBasedInitializer->InitializeTransform();

      //HACK:  THIS IS WRONG!  m_VersorTransform was used to set movingLmks! initToAtlasAffineTransform->Compose( this->m_VersorTransform, true );
      typedef itk::BRAINSFitHelper HelperType;
      HelperType::Pointer brainsFitHelper = HelperType::New();

      // Now Run BRAINSFitHelper class initialized with initToAtlasAffineTransform, original image, and atlas image
      // adapted from BRAINSABC/brainseg/AtlasRegistrationMethod.hxx - do I need to change any of these parameters?
      brainsFitHelper->SetNumberOfSamples(500000);
      brainsFitHelper->SetNumberOfHistogramBins(50);
      std::vector<int> numberOfIterations(1);
      numberOfIterations[0] = 1500;
      brainsFitHelper->SetNumberOfIterations(numberOfIterations);
      brainsFitHelper->SetTranslationScale(1000);
      brainsFitHelper->SetReproportionScale(1.0);
      brainsFitHelper->SetSkewScale(1.0);

      typedef itk::Image<float, 3>                            FloatImageType;
      typedef itk::CastImageFilter<SImageType, FloatImageType> CastFilterType;

        {
        typename CastFilterType::Pointer fixedCastFilter = CastFilterType::New();
        fixedCastFilter->SetInput( atlasReader->GetOutput() );
        fixedCastFilter->Update();
        brainsFitHelper->SetFixedVolume( fixedCastFilter->GetOutput() );

        typename CastFilterType::Pointer movingCastFilter = CastFilterType::New();
        movingCastFilter->SetInput( this->GetInput() );
        movingCastFilter->Update();
        brainsFitHelper->SetMovingVolume( movingCastFilter->GetOutput() );
        }

      std::vector<double> minimumStepSize(1);
      minimumStepSize[0] = 0.005;
      brainsFitHelper->SetMinimumStepLength(minimumStepSize);
      std::vector<std::string> transformType(1);
      transformType[0] = "Affine";
      brainsFitHelper->SetTransformType(transformType);

      brainsFitHelper->SetCurrentGenericTransform( initToAtlasAffineTransform.GetPointer() );
      brainsFitHelper->Update();

      this->m_VersorTransform =
        itk::ComputeRigidTransformFromGeneric( brainsFitHelper->GetCurrentGenericTransform().GetPointer() );
      if( this->m_VersorTransform.IsNull() )
        {
        // Fail if something weird happens.  TODO: This should throw an exception.
        std::cout << "this->m_VersorTransform is null. It means we're not registering to the atlas, after all."
          << std::endl;
        std::cout << "FAILIING" << std::endl;
        exit(-1);
        }

      //TODO: HACK:  Translate found ACPoint
#if 0
      // as a final step, translate the AC back to the origin.
        {
        LandmarkConstIterator                           acIter = acpcLandmarks.find( "AC" );
        const VersorRigid3DTransformType::OutputPointType acOrigPoint =
          "A transform of some sort here"->TransformPoint( acIter->second );

        VersorTransformType::Pointer finalTransform = VersorTransformType::New();
        finalTransform->SetFixedParameters( this->m_VersorTransform->GetFixedParameters() );
        finalTransform->SetParameters( this->m_VersorTransform->GetParameters() );
        finalTransform->GetInverse( invFinalTransform );

        // TODO:  CHECK if this can be less convoluted. Too many inverses used.  translate the forward by positive
        // rather than inverse by negative.
        //
        VersorRigid3DTransformType::OutputPointType acPoint = invFinalTransform->TransformPoint( acOrigPoint );
          {
          VersorRigid3DTransformType::OffsetType translation;
          translation[0] = -acPoint[0];
          translation[1] = -acPoint[1];
          translation[2] = -acPoint[2];
          invFinalTransform->Translate( translation, true );
          }
        invFinalTransform->GetInverse( finalTransform );

        // TODO: Remove VersorRigid3DTransformType::OutputPointType acFinalPoint =  invFinalTransform->TransformPoint (
        // acOrigPoint );
        }
#endif
      }
  ///END BRAINSFIT_ALTERNATIVE
  ////////////////////////////

  if( LMC::globalverboseFlag )
    {
    std::cout << "VersorRotation: " << this->m_VersorTransform->GetMatrix() << std::endl;
    std::cout << "itkVersorRigid3DTransform Parameters: " << this->m_VersorTransform->GetParameters() << std::endl;
    std::cout << "itkVersorRigid3DTransform FixedParameters: " << this->m_VersorTransform->GetFixedParameters()
      << std::endl;
    std::cout << "itkVersorRigid3DTransform GetCenter(): " << this->m_VersorTransform->GetCenter() << std::endl;
    std::cout << "itkVersorRigid3DTransform GetTranslation(): " << this->m_VersorTransform->GetTranslation()
                                                                   << std::endl;
    std::cout << "itkVersorRigid3DTransform GetMatrix(): " << this->m_VersorTransform->GetMatrix()
                                                              << std::endl;

    std::cout << "itkRigid3DTransform Parameters: " << this->m_VersorTransform->GetParameters() << std::endl;
    std::cout << "itkRigid3DTransform FixedParameters: " << this->m_VersorTransform->GetFixedParameters()
      << std::endl;
    std::cout << "itkRigid3DTransform GetCenter(): " << this->m_VersorTransform->GetCenter() << std::endl;
    std::cout << "itkRigid3DTransform GetTranslation(): " << this->m_VersorTransform->GetTranslation() << std::endl;
    std::cout << "itkRigid3DTransform GetMatrix(): " << this->m_VersorTransform->GetMatrix()
                                                        << std::endl;
    std::cout << "itkVersorRigid3DTransform: \n" <<  this->m_VersorTransform << std::endl;
    std::cout << "itkRigid3DTransform: \n" <<  this->m_VersorTransform << std::endl;
    }

  itk::PrepareOutputImages(this->m_OutputResampledImage,
    this->m_OutputImage,
    this->m_OutputUntransformedClippedVolume,
    myDetector.GetOriginalInput().GetPointer(), //Input RO
    this->m_VersorTransform.GetPointer(), //Input RO
    this->m_AcLowerBound, //Input RO
    BackgroundFillValue, //Input RO
    this->m_InterpolationMode, //Input RO
    this->m_CutOutHeadInOutputVolume, //Input RO
    this->m_OtsuPercentileThreshold //Input RO
  );

  itk::PrepareOutputLandmarks(
    this->m_VersorTransform.GetPointer(), //Input RO
    myDetector.GetNamedPoints(),
    this->m_AlignedPoints
  );

  if( globalImagedebugLevel > 3 )
    {
    SImageType::Pointer TaggedOriginalImage = myDetector.GetTaggedImage();
    itkUtil::WriteImage<SImageType>(TaggedOriginalImage, this->m_ResultsDir + "/TAGGED_POINTS.nii.gz");
      {
      SImageType::Pointer isoTaggedImage =
        TransformResample<SImageType, SImageType>(
          TaggedOriginalImage.GetPointer(), MakeIsoTropicReferenceImage().GetPointer(), BackgroundFillValue,
          GetInterpolatorFromString<SImageType>("Linear").GetPointer(), this->m_VersorTransform.GetPointer() );
      itkUtil::WriteImage<SImageType>(isoTaggedImage, this->m_ResultsDir + "/ISO_Lmk_MSP.nii.gz");
      }
      {
      SImageType::Pointer VersorisoTaggedImage =
        TransformResample<SImageType, SImageType>(
          TaggedOriginalImage.GetPointer(), MakeIsoTropicReferenceImage().GetPointer(), BackgroundFillValue,
          GetInterpolatorFromString<SImageType>("Linear").GetPointer(), this->m_VersorTransform.GetPointer() );
      itkUtil::WriteImage<SImageType>(VersorisoTaggedImage, this->m_ResultsDir + "/Versor_ISO_Lmk_MSP.nii.gz");
      }
      {
      RigidTransformType::Pointer OrigSpaceCenterOfGravityCentered = myDetector.GetTransformToMSP();
      SImageType::Pointer         RigidMSPImage =
        TransformResample<SImageType, SImageType>(
          TaggedOriginalImage.GetPointer(), MakeIsoTropicReferenceImage().GetPointer(), BackgroundFillValue,
          GetInterpolatorFromString<SImageType>("Linear").GetPointer(), OrigSpaceCenterOfGravityCentered.GetPointer() );
      itkUtil::WriteImage<SImageType>(RigidMSPImage, this->m_ResultsDir + "/RigidMSPImage_Lmk_MSP.nii.gz");
      }
    }
  if( this->m_WriteBranded2DImage.compare("") != 0 )
    {
    MakeBranded2DImage(this->m_OutputResampledImage.GetPointer(), myDetector,
      this->m_AlignedPoints["RP"],
      this->m_AlignedPoints["AC"],
      this->m_AlignedPoints["PC"],
      this->m_AlignedPoints["VN4"],
      this->m_AlignedPoints["CM"],
      this->m_WriteBranded2DImage);
    }
  this->GraftOutput(this->m_OutputImage);
}

template <class TInputImage, class TOutputImage>
void
BRAINSConstellationDetector2<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

}
