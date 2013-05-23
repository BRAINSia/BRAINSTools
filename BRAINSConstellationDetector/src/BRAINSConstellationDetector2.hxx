/*=========================================================================
 Author: Wei Lu
 at Psychiatry Imaging Lab,
 University of Iowa Health Care 2010

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 =========================================================================*/

#include "BRAINSConstellationDetector2.h"
#include "GenericTransformImage.h"
#include "itkOrthogonalize3DRotationMatrix.h"

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
  this->m_InvVersorTransform = NULL;
  this->m_OutputImage = NULL;
  this->m_OutputResampledImage = NULL;
  this->m_OutputUntransformedClippedVolume = NULL;
  this->m_ClippingFactorImage = NULL ;

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
  char modelfile[1024];
    {
    strcpy( modelfile, this->m_InputTemplateModel.c_str() );
    }

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
  myModel.ReadModelFile( std::string(modelfile) );

  if( LMC::globalverboseFlag )
    {
    std::cout << "Using Model File: " << modelfile << std::endl;
    myModel.PrintHeaderInfo();
    }

  // Override some landmarks by user
  vnl_vector<double> templateRadius;    // in units of mm
  templateRadius.set_size(4);
  templateRadius[0] = ( this->m_RadiusMPJ <= 0 ) ? myModel.GetRadius("RP") : this->m_RadiusMPJ; //
                                                                                                //
                                                                                                // =
                                                                                                //
                                                                                                // RP
                                                                                                //
                                                                                                // radius
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
  SImageType::Pointer image;
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

    image = remapIntensityFilter->GetOutput();
    }
  else
    {
    image = volOrig;
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
    duplicator->SetInputImage( image );
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

    this->m_InvVersorTransform = VersorTransformType::New();
    const SImageType::PointType centerPoint = this->m_VersorTransform->GetCenter();
    this->m_InvVersorTransform->SetCenter(centerPoint);
    this->m_InvVersorTransform->SetIdentity();
    this->m_VersorTransform->GetInverse(this->m_InvVersorTransform);

    // std::cout << "versor         transform parameters are " <<
    // ZeroCenteredTransform->GetParameters() << std::endl;
    // std::cout << "versor inverse transform parameters are " <<
    // this->m_InvVersorTransform->GetParameters() << std::endl;
    // std::cout << "               transform parameters are nice, huh?" <<
    // std::endl;
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

    if( globalImagedebugLevel > 3 )
      {
      SImageType::Pointer TaggedOriginalImage = myDetector.GetTaggedImage();
      itkUtil::WriteImage<SImageType>(TaggedOriginalImage, this->m_ResultsDir + "/TAGGED_POINTS.nii.gz");
        {
        SImageType::Pointer isoTaggedImage =
          TransformResample<SImageType, SImageType>(
            TaggedOriginalImage, MakeIsoTropicReferenceImage(), BackgroundFillValue,
            GetInterpolatorFromString<SImageType>("Linear"), this->m_VersorTransform.GetPointer() );
        itkUtil::WriteImage<SImageType>(isoTaggedImage, this->m_ResultsDir + "/ISO_Lmk_MSP.nii.gz");
        }
        {
        SImageType::Pointer VersorisoTaggedImage =
          TransformResample<SImageType, SImageType>(
            TaggedOriginalImage, MakeIsoTropicReferenceImage(), BackgroundFillValue,
            GetInterpolatorFromString<SImageType>("Linear"), this->m_VersorTransform.GetPointer() );
        itkUtil::WriteImage<SImageType>(VersorisoTaggedImage, this->m_ResultsDir + "/Versor_ISO_Lmk_MSP.nii.gz");
        }
        {
        RigidTransformType::Pointer OrigSpaceCenterOfGravityCentered = myDetector.GetTransformToMSP();
        SImageType::Pointer         RigidMSPImage =
          TransformResample<SImageType, SImageType>(
            TaggedOriginalImage, MakeIsoTropicReferenceImage(), BackgroundFillValue,
            GetInterpolatorFromString<SImageType>("Linear"), OrigSpaceCenterOfGravityCentered.GetPointer() );
        itkUtil::WriteImage<SImageType>(RigidMSPImage, this->m_ResultsDir + "/RigidMSPImage_Lmk_MSP.nii.gz");
        }
      }

    double PhysicalLowerBound = /* ACy when zero-centered is ... */ 0.0 - this->m_AcLowerBound;
    for( LandmarksMapType::const_iterator lit = myDetector.GetNamedPoints().begin();
         lit != myDetector.GetNamedPoints().end(); ++lit )
      {
      this->m_AlignedPoints[lit->first] =
        this->m_InvVersorTransform->TransformPoint
          ( myDetector.GetOriginalSpaceNamedPoint(lit->first) );
      }

      {
      image = myDetector.GetOriginalInput();   // image -> myDetector(modification May happen) -> image

      this->m_ImageToBeResampled = image ;
        {
        // const SImageType * constImage( this->m_OriginalInputImage.GetPointer() );
        const SImageType * constImage( image.GetPointer() );

        ResampleIPFilterPointer resampleIPFilter = ResampleIPFilterType::New();
        resampleIPFilter->SetInputImage( constImage );
        resampleIPFilter->SetRigidTransform( m_VersorTransform.GetPointer() );
        resampleIPFilter->Update();
        this->m_OutputImage = resampleIPFilter->GetOutput();
        this->GraftOutput(m_OutputImage);
        }

        {
        this->m_OutputResampledImage = TransformResample<SImageType, SImageType>( image,
                                                                                  MakeIsoTropicReferenceImage(),
                                                                                  BackgroundFillValue,
                                                                                  GetInterpolatorFromString<SImageType>(
                                                                                    this->m_InterpolationMode),
                                                                                  this->m_VersorTransform.GetPointer() );
        }

        {
        // ======================== Start
        // ======================== Start
        //  HACK -- chopping based on AcLowerBound
        //  This is ugly code that could be re-written much simpler.
        //
        typedef itk::ImageRegionIteratorWithIndex<SImageType> IteratorType;
        const double thousand = 1000.0;   // we need a DOUBLE constant, not a
        // FLOAT constant, for exact switch
        // comparisons.
        if( this->m_AcLowerBound < thousand )
          {
          // First Process the OutputResampledImage
          std::cout << "Chopping image below physical location: " << PhysicalLowerBound << "." << std::endl;
          ChopImageBelowLowerBound<SImageType>(this->m_OutputResampledImage, BackgroundFillValue, PhysicalLowerBound);
          ChopImageBelowLowerBound<SImageType>(this->m_OutputImage, BackgroundFillValue, PhysicalLowerBound);

          // Second Create a mask for inverse resampling to orignal space
          SImageType::Pointer ZeroOneImage = SImageType::New();
          ZeroOneImage->CopyInformation(this->m_OutputResampledImage);
          // don't forget to do SetRegions here!
          ZeroOneImage->SetRegions( this->m_OutputResampledImage->GetLargestPossibleRegion() );
          ZeroOneImage->Allocate();
          ZeroOneImage->FillBuffer(1);
          ChopImageBelowLowerBound<SImageType>(ZeroOneImage, BackgroundFillValue, PhysicalLowerBound);

          if( this->m_CutOutHeadInOutputVolume )  // Restrict mask to head
                                                  // tissue region if necessary
            {
            //  No double opportunity when generating both kinds of images.
            const unsigned int closingSize = 7;
            //            SImageType::Pointer HeadOutlineMaskImage =
            // FindLargestForgroundFilledMask<SImageType>(
            // this->m_OutputResampledImage, this->m_OtsuPercentileThreshold,
            // closingSize );
            typedef itk::LargestForegroundFilledMaskImageFilter<SImageType> LFFMaskFilterType;
            LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
            LFF->SetInput(this->m_OutputResampledImage);
            LFF->SetOtsuPercentileThreshold(this->m_OtsuPercentileThreshold);
            LFF->SetClosingSize(closingSize);
            LFF->Update();
            SImageType::Pointer HeadOutlineMaskImage =  LFF->GetOutput();

            IteratorType ItZeroOneImage( ZeroOneImage, ZeroOneImage->GetRequestedRegion() );
            ItZeroOneImage.GoToBegin();
            IteratorType ItOutputResampledImage( this->m_OutputResampledImage,
                                                 this->m_OutputResampledImage->GetRequestedRegion() );
            ItOutputResampledImage.GoToBegin();
            IteratorType ItHead( HeadOutlineMaskImage, HeadOutlineMaskImage->GetLargestPossibleRegion() );
            ItHead.GoToBegin();
            while( !ItHead.IsAtEnd() )
              {
              if( ItHead.Get() == 0 )
                {
                ItOutputResampledImage.Set(0);
                ItZeroOneImage.Set(0);
                }
              ++ItZeroOneImage;
              ++ItOutputResampledImage;
              ++ItHead;
              }
            }
          // Map the ZeroOne image through the inverse zero-centered transform
          // to
          // make the clipping factor image:
          typedef itk::NearestNeighborInterpolateImageFunction<SImageType, double> NearestNeighborInterpolatorType;
          // typename
          NearestNeighborInterpolatorType::Pointer interpolator = NearestNeighborInterpolatorType::New();
          typedef itk::ResampleImageFilter<SImageType, SImageType> ResampleFilterType;
          // typename
          ResampleFilterType::Pointer ResampleFilter = ResampleFilterType::New();
          ResampleFilter->SetInput(ZeroOneImage);
          ResampleFilter->SetInterpolator(interpolator);
          ResampleFilter->SetDefaultPixelValue(0);
          ResampleFilter->SetOutputParametersFromImage(image);

          ResampleFilter->SetTransform(this->m_InvVersorTransform);
          ResampleFilter->Update();
          // typename
          this->m_ClippingFactorImage = ResampleFilter->GetOutput();
          // itkUtil::WriteImage<SImageType>(clippingFactorImage,"ClippingImage.nii.gz");

          // Multiply the raw input image by the clipping factor image:
          typedef itk::MultiplyImageFilter<SImageType, SImageType> MultiplyFilterType;
          // typename
          MultiplyFilterType::Pointer MultiplyFilter = MultiplyFilterType::New();
          MultiplyFilter->SetInput1(image);
          MultiplyFilter->SetInput2(this->m_ClippingFactorImage);
          MultiplyFilter->Update();
          this->m_OutputUntransformedClippedVolume = MultiplyFilter->GetOutput();
          }
        // ======================== Stop
        // ======================== Stop
        }
      if( this->m_WriteBranded2DImage.compare("") != 0 )
        {
        /*
        for( LandmarksMapType::const_iterator lit = m_AlignedPoints.begin(); lit != m_AlignedPoints.end(); ++lit )
        {
            std::cout << lit->first << "=" << lit->second[0] << "," << lit->second[1] << "," << lit->second[2] << std::endl << std::endl;
        }

        std::cout << "filename : " << this->m_WriteBranded2DImage << std::endl;

        itkUtil::WriteImage<SImageType>(this->m_OutputResampledImage, "m_OutputResampledImage.nii.gz");

        double height = myDetector.GetModelHeight("AC");
        */

        MakeBranded2DImage(this->m_OutputResampledImage.GetPointer(), myDetector,
                           this->m_AlignedPoints["RP"],
                           this->m_AlignedPoints["AC"],
                           this->m_AlignedPoints["PC"],
                           this->m_AlignedPoints["VN4"],
                           this->m_AlignedPoints["CM"],
                           this->m_WriteBranded2DImage);
        }
      }     //  scope of this->m_OutputResampledImage
    }
}

template <class TInputImage, class TOutputImage>
void
BRAINSConstellationDetector2<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
}
