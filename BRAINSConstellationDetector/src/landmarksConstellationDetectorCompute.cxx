
#include "landmarksConstellationDetector.h"
// landmarkIO has to be included after landmarksConstellationDetector
#include "landmarkIO.h"
#include "itkOrthogonalize3DRotationMatrix.h"

#include "itkFindCenterOfBrainFilter.h"
#include "BRAINSHoughEyeDetector.h"

#include <BRAINSFitHelper.h>
#include "itkLandmarkBasedTransformInitializer.h"

/**
 * Diagnostic debugging function to print the differences between the landmarks of the ref to the items in cmp
 * @param ref
 * @param cmp
 * @param backward_compare use cmp as ref, and ref as cmp
 * @param filename __FILE__
 * @param lineno __LINE__
 * @return
 */
static bool
lmk_check_differences( const LandmarksMapType & ref, const LandmarksMapType & cmp, const bool backward_compare,
                       std::string filename, int lineno )
{
  std::cerr << "===" << lineno << " " << filename << std::endl;
  constexpr float landmark_sameness_tolerence = 0.001;
  bool            ref_lmks_same_in_cmp = true;
  for ( const auto lmk : ref )
  {
    const auto cmp_lmk_iter = cmp.find( lmk.first );
    if ( cmp_lmk_iter == cmp.end() )
    {
      std::cerr << "MISSING: " << lmk.first << std::endl;
      ref_lmks_same_in_cmp = false;
    }
    else if ( lmk.second.EuclideanDistanceTo( cmp_lmk_iter->second ) >= landmark_sameness_tolerence )
    {
      std::cerr << "ERROR:  " << lmk.first << " " << lmk.second << " != " << cmp_lmk_iter->second << std::endl;
      ref_lmks_same_in_cmp = false;
    }
    else
    {
      std::cout << "SAME: " << lmk.first << std::endl;
    }
  }
  if ( backward_compare )
  {
    std::cerr << "--- Backwards Compare" << std::endl;
    ref_lmks_same_in_cmp = lmk_check_differences( cmp, ref, false, filename, lineno );
  }
  const char * const status = ( ref_lmks_same_in_cmp ) ? "SUCCESS" : "FAILED";
  std::cerr << "============" << status << std::endl;
  return ref_lmks_same_in_cmp;
}

void
landmarksConstellationDetector::Compute( SImageType::Pointer orig_space_image )
{

  std::cout << "\nEstimating MSP..." << std::endl;


  if ( globalImagedebugLevel > 2 )
  {
    LandmarksMapType eyeFixed_lmks;
    eyeFixed_lmks["CM"] = this->m_eyeFixed_lmk_CenterOfHeadMass;
    const std::string roughlyAlignedCHMName( this->m_ResultsDir + "/eyeFixed_lmks.fcsv" );
    WriteITKtoSlicer3Lmk( roughlyAlignedCHMName, eyeFixed_lmks );

    const std::string roughlyAlignedVolumeName( this->m_ResultsDir + "/VolumeRoughAlignedWithHoughEye.nrrd" );
    itkUtil::WriteImage< SImageType >( this->m_eyeFixed_img, roughlyAlignedVolumeName );
  }

  // Compute the estimated MSP transform, and aligned image
  double c_c = 0;
  ComputeMSP( this->m_eyeFixed_img,
              this->m_test_orig2msp_img_tfm,
              this->m_msp_img,
              this->m_eyeFixed_lmk_CenterOfHeadMass,
              this->m_mspQualityLevel,
              c_c );
#if 0
  /*
   * If the MSP estimation is not good enough (i.e. c_c < -0.64), we try to compute the MSP again
   * using the previous m_test_orig2msp_img_tfm as an initial transform.
   * Threshold is set to 0.64 based on the statistical calculations on 23 successfully passed data.
   */
  if( c_c > -0.64 && !this->m_HoughEyeFailure )
    {
      unsigned int maxNumberOfIterations = 5;
      std::cout << "\n============================================================="
                << "\nBad Estimation for MSP Plane.\n"
                << "Repeat the Estimation Process up to " << maxNumberOfIterations
                << " More Times to Find a Better Estimation..." << std::endl;

      for (unsigned int i = 0; i<maxNumberOfIterations; i++)
        {
        std::cout << "\nTry " << i+1 << "..." << std::endl;

        // Rotate VolumeRoughAlignedWithHoughEye by finalTmsp again.
        SImageType::Pointer localRoughAlignedInput;
        DoResampleInPlace( this->m_eyeFixed_img.GetPointer(), // input
                           this->m_test_orig2msp_img_tfm.GetPointer(), // input
                           localRoughAlignedInput ); // output

        // Transform the eyeFixed_lmk_CenterOfHeadMass by finlaTmsp transform.
        VersorTransformType::Pointer localInvFinalTmsp = VersorTransformType::New();
        this->m_test_orig2msp_img_tfm->GetInverse( localInvFinalTmsp );
        SImageType::PointType localAlignedCOHM = localInvFinalTmsp->TransformPoint( this->m_eyeFixed_lmk_CenterOfHeadMass );

        if( globalImagedebugLevel > 2 )
          {
          LandmarksMapType locallyRotatedCHM;
          locallyRotatedCHM["CM"] = localAlignedCOHM;
          const std::string localRotatedCenterOfHeadMassNAME
                            ( this->m_ResultsDir + "/localRotatedCHM_" + local_to_string(i) + "_.fcsv" );
          WriteITKtoSlicer3Lmk( localRotatedCenterOfHeadMassNAME, locallyRotatedCHM );

          const std::string localRotatedRoughAligedVolumeNAME
                            ( this->m_ResultsDir + "/localRotatedRoughAlignedVolume_" + local_to_string(i) + "_.nrrd" );
          itkUtil::WriteImage<SImageType> ( localRoughAlignedInput, localRotatedRoughAligedVolumeNAME );
          }

        // Compute MSP again
        RigidTransformType::Pointer localTmsp;
        c_c = 0;
        ComputeMSP( localRoughAlignedInput, localTmsp,
                   this->m_msp_img, localAlignedCOHM, this->m_mspQualityLevel, c_c );

        // Compose the eyeFixed2msp_lmk_tfm transforms
        this->m_test_orig2msp_img_tfm->Compose(localTmsp);

        if ( c_c < -0.64 )
          {
          break;
          }
        }

      using StatisticsFilterType = itk::StatisticsImageFilter<SImageType>;
      StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
      statisticsFilter->SetInput(this->m_eyeFixed_img);
      statisticsFilter->Update();
      SImageType::PixelType minPixelValue = statisticsFilter->GetMinimum();

      this->m_msp_img = TransformResample<SImageType, SImageType>( this->m_eyeFixed_img.GetPointer(),
                                                                    MakeIsoTropicReferenceImage().GetPointer(),
                                                                    minPixelValue,
                                                                    GetInterpolatorFromString<SImageType>("Linear").GetPointer(),
                                                                    this->GetTransformToMSP().GetPointer() );
    }
#endif
  std::cout << "\n=============================================================" << std::endl;

  // Generate a warning if reflective correlation similarity measure is low.
  // It may be normal in some very diseased subjects, so don't throw an exception here.
  if ( c_c > -0.64 )
  {
    std::cout << "WARNING: Low reflective correlation between left/right hemispheres." << std::endl
              << "The estimated landmarks may not be reliable.\n"
              << std::endl;
  }

  // Throw an exception and stop BCD if RC metric is too low (less than 0.4) because results will not be reliable.
  if ( c_c > -0.40 )
  {
    itkGenericExceptionMacro( << "Too large MSP estimation error! reflective correlation metric is: " << c_c
                              << std::endl
                              << "Estimation of landmarks will not be reliable.\n"

                              << std::endl );
  }

  // In case hough eye detector failed
  if ( this->m_HoughEyeFailure || ( globalImagedebugLevel > 1 ) )
  {
    const std::string EMSP_Fiducial_file_name( "EMSP.fcsv" );
    std::stringstream failureMessageStream( "" );
    failureMessageStream << "EMSP aligned image and zero eye centers "
                         << "landmarks are written to " << std::endl
                         << this->m_ResultsDir << ". Use GUI corrector to "
                         << "correct the landmarks in." << EMSP_Fiducial_file_name
                         << "INITIAL LMKS: " << EMSP_Fiducial_file_name << "FOR IMAGE: " << this->m_msp_img
                         << "IN DIR: " << this->m_ResultsDir << std::endl;
    LandmarksMapType zeroEyeCenters;
    if ( this->m_HoughEyeFailure )
    {
      SImageType::PointType zeroPoint;
      zeroPoint.Fill( 0 );

      zeroEyeCenters["LE"] = zeroPoint;
      zeroEyeCenters["RE"] = zeroPoint;
    }
    WriteManualFixFiles( EMSP_Fiducial_file_name,
                         this->m_msp_img,
                         this->m_ResultsDir,
                         zeroEyeCenters,
                         failureMessageStream.str(),
                         this->m_HoughEyeFailure );
  }

  if ( globalImagedebugLevel > 2 )
  {
    const std::string MSP_ImagePlaneForOrig( this->m_ResultsDir + "/MSP_PLANE_For_Orig.nii.gz" );
    CreatedebugPlaneImage( this->m_eyeFixed_img, MSP_ImagePlaneForOrig );
    const std::string MSP_ImagePlane( this->m_ResultsDir + "/MSP_PLANE.nii.gz" );
    CreatedebugPlaneImage( this->m_msp_img, MSP_ImagePlane );
  }

  // HACK:  TODO:  DEBUG: Need to remove redundant need for volume and mask when
  // they can be the same image ( perhaps with a threshold );
  SImageType::Pointer mask_LR = this->m_msp_img;

  VersorTransformType::Pointer orig2eyeFixed_lmk_tfm;
  orig2eyeFixed_lmk_tfm = VersorTransformType::New();
  this->m_orig2eyeFixed_img_tfm->GetInverse( orig2eyeFixed_lmk_tfm );

  // TODO:  ERROR
  VersorTransformType::Pointer eyeFixed2msp_lmk_tfm = VersorTransformType::New();
  this->m_test_orig2msp_img_tfm->GetInverse( eyeFixed2msp_lmk_tfm );

  this->m_msp_lmk_CenterOfHeadMass = eyeFixed2msp_lmk_tfm->TransformPoint( this->m_eyeFixed_lmk_CenterOfHeadMass );
  this->m_msp_lmk_CenterOfHeadMass[0] = 0; // Search starts on the estimated MSP

  {
    SImageType::PointType msp_lmk_RP_Candidate;
    SImageType::PointType msp_lmk_PC_Candiate;
    SImageType::PointType msp_lmk_VN4_Candidate;
    SImageType::PointType msp_lmk_AC_Candidate;

    SImageType::PointType::VectorType RPtoCEC;
    SImageType::PointType::VectorType RPtoVN4;
    SImageType::PointType::VectorType RPtoAC;
    SImageType::PointType::VectorType RPtoPC;

    SImageType::PointType mspSpaceCEC;
    //{
    SImageType::PointType eyeFixed_LE = eyeFixed2msp_lmk_tfm->TransformPoint(
      orig2eyeFixed_lmk_tfm->TransformPoint( this->m_orig_lmks_constant.at( "LE" ) ) );
    SImageType::PointType eyeFixed_RE = eyeFixed2msp_lmk_tfm->TransformPoint(
      orig2eyeFixed_lmk_tfm->TransformPoint( this->m_orig_lmks_constant.at( "RE" ) ) );


    mspSpaceCEC.SetToMidPoint( eyeFixed_LE, eyeFixed_RE );
    mspSpaceCEC[0] = 0; // Search starts on the estimated MSP
    //}
    std::cout << "\nPerforming morpohmetric search + local search..." << std::endl;

    {
      /*
       * Search for MPJ ( RP )
       *
       * Use the knowledge of mean CMtoRP vector + local searching
       */
      std::cout << "Processing MPJ..." << std::endl;
      double searchRadiusLR = 4.; // in mm, and it is only for landmarks near
      // MSP
      double cc_RP_Max = 0;


      // TODO: Remove these constants.  just query the map.
      // save the result that whether we are going to process all the landmarks
      // in light of user-specified eye center info.
      const bool hasUserForcedRPPoint = mapHasKey( m_orig_lmks_constant, "RP" );
      const bool hasUserForcedACPoint = mapHasKey( m_orig_lmks_constant, "AC" );
      const bool hasUserForcedPCPoint = mapHasKey( m_orig_lmks_constant, "PC" );
      const bool hasUserForcedVN4Point = mapHasKey( m_orig_lmks_constant, "VN4" );

      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
      if ( hasUserForcedRPPoint )
      {
        std::cout << "Skip estimation of RP, directly forced by command line." << std::endl;
        msp_lmk_RP_Candidate = eyeFixed2msp_lmk_tfm->TransformPoint(
          orig2eyeFixed_lmk_tfm->TransformPoint( this->m_orig_lmks_constant.at( "RP" ) ) );
      }
      else
      {
        // The search radius of RP is set to 5 times larger than its template
        // radius in SI direction.
        // It's large because some scans have extra contents of neck or shoulder
        msp_lmk_RP_Candidate = FindCandidatePoints( this->m_msp_img,
                                                    mask_LR,
                                                    searchRadiusLR,
                                                    3. * this->m_TemplateRadius["RP"],
                                                    5. * this->m_TemplateRadius["RP"],
                                                    this->m_msp_lmk_CenterOfHeadMass.GetVectorFromOrigin() +
                                                      this->m_InputTemplateModel.GetCMtoRPMean(),
                                                    this->m_InputTemplateModel.GetTemplateMeans( "RP" ),
                                                    this->m_InputTemplateModel.m_VectorIndexLocations["RP"],
                                                    cc_RP_Max,
                                                    "RP" );
      }

      // Local search radius in LR direction is affected by the
      // estimated MSP error in LR direction
      const double err_MSP = std::abs( msp_lmk_RP_Candidate[0] - this->m_msp_lmk_CenterOfHeadMass[0] );
      std::cout << "The estimated MSP error in LR direction: " << err_MSP << " mm" << std::endl;

      if ( err_MSP < 1 )
      {
        searchRadiusLR = 1.;
        std::cout << "Local search radius in LR direction is set to 1 mm." << std::endl;
      }
      else if ( err_MSP < 2 )
      {
        searchRadiusLR = 2.;
        std::cout << "Local search radius in LR direction is set to 2 mm." << std::endl;
      }
      else if ( err_MSP > 6 )
      {
        itkGenericExceptionMacro( << "Bad MPJ estimation or too large MSP estimation error!"
                                  << "The estimation result is probably not reliable." );
      }
      else // not changed
      {
        std::cout << "Local search radius in LR direction is set to 4 mm." << std::endl;
      }

      /*
       * Search for 4VN
       *
       * we will use RPtoCEC ( RP-to-center of eye centers mean vector ) and RPtoVN4Mean
       * to determine the search center for VN4
       *
       * Assumption:
       * 1. Human brains share a similar angle from RPtoCEC to RPtoVN4
       * 2. Human brains share a similar RPtoVN4 norm
       * 3. The center of eye centers, MPJ, and AC are all very close to the MSP plane
       * 4. VN4 is always below CEC-MPJ line on MSP plane
       */
      std::cout << "Processing 4VN..." << std::endl;
      RPtoCEC = mspSpaceCEC - msp_lmk_RP_Candidate;

      // RPtoVN4 = this->m_InputTemplateModel.GetRPtoXMean( "VN4" );
      RPtoVN4 = FindVectorFromPointAndVectors(
        RPtoCEC, this->m_InputTemplateModel.GetRPtoCECMean(), this->m_InputTemplateModel.GetRPtoXMean( "VN4" ), -1 );

      double cc_VN4_Max = 0;

      if ( hasUserForcedVN4Point )
      {
        std::cout << "Skip estimation of VN4, directly forced by command line." << std::endl;
        msp_lmk_VN4_Candidate = eyeFixed2msp_lmk_tfm->TransformPoint(
          orig2eyeFixed_lmk_tfm->TransformPoint( this->m_orig_lmks_constant.at( "VN4" ) ) );
      }
      else
      {
        msp_lmk_VN4_Candidate = FindCandidatePoints( this->m_msp_img,
                                                     mask_LR,
                                                     searchRadiusLR,
                                                     1.6 * this->m_TemplateRadius["VN4"],
                                                     1.6 * this->m_TemplateRadius["VN4"],
                                                     msp_lmk_RP_Candidate.GetVectorFromOrigin() + RPtoVN4,
                                                     this->m_InputTemplateModel.GetTemplateMeans( "VN4" ),
                                                     this->m_InputTemplateModel.m_VectorIndexLocations["VN4"],
                                                     cc_VN4_Max,
                                                     "VN4" );
      }

      /*
       * Search for AC
       *
       * Assumption:
       * 1. Human brains share a similar angle from RPtoCEC to RPtoAC
       * 2. Human brains share a similar RPtoAC norm
       * 3. The center of eye centers, MPJ, and AC are all very close to the MSP plane
       * 4. AC is always above CEC-MPJ line on MSP plane
       */
      std::cout << "Processing AC..." << std::endl;

      // RPtoAC = this->m_InputTemplateModel.GetRPtoXMean( "AC" );
      RPtoAC = FindVectorFromPointAndVectors(
        RPtoCEC, this->m_InputTemplateModel.GetRPtoCECMean(), this->m_InputTemplateModel.GetRPtoXMean( "AC" ), 1 );
      double cc_AC_Max = 0;
      if ( hasUserForcedACPoint )
      {
        std::cout << "Skip estimation of AC , directly forced by command line." << std::endl;
        msp_lmk_AC_Candidate = eyeFixed2msp_lmk_tfm->TransformPoint(
          orig2eyeFixed_lmk_tfm->TransformPoint( this->m_orig_lmks_constant.at( "AC" ) ) );
      }
      else
      {
        msp_lmk_AC_Candidate = FindCandidatePoints( this->m_msp_img,
                                                    mask_LR,
                                                    searchRadiusLR,
                                                    1.6 * this->m_TemplateRadius["AC"],
                                                    1.6 * this->m_TemplateRadius["AC"],
                                                    msp_lmk_RP_Candidate.GetVectorFromOrigin() + RPtoAC,
                                                    this->m_InputTemplateModel.GetTemplateMeans( "AC" ),
                                                    this->m_InputTemplateModel.m_VectorIndexLocations["AC"],
                                                    cc_AC_Max,
                                                    "AC" );
      }

      /*
       * Search for PC
       *
       * Assumption:
       * 1. Human brains share a similar angle from RPtoCEC to RPtoPC
       * 2. Human brains share a similar RPtoPC norm
       * 3. The center of eye centers, MPJ, and PC are all very close to the MSP plane
       * 4. PC is always above CEC-MPJ line on MSP plane
       */
      std::cout << "Processing PC..." << std::endl;

      // RPtoPC = this->m_InputTemplateModel.GetRPtoXMean( "PC" );
      RPtoPC = FindVectorFromPointAndVectors(
        RPtoCEC, this->m_InputTemplateModel.GetRPtoCECMean(), this->m_InputTemplateModel.GetRPtoXMean( "PC" ), 1 );
      double cc_PC_Max = 0;
      if ( hasUserForcedPCPoint )
      {
        std::cout << "Skip estimation of PC, directly forced by command line." << std::endl;
        msp_lmk_PC_Candiate = eyeFixed2msp_lmk_tfm->TransformPoint(
          orig2eyeFixed_lmk_tfm->TransformPoint( this->m_orig_lmks_constant.at( "PC" ) ) );
      }
      else
      {
        msp_lmk_PC_Candiate = FindCandidatePoints( this->m_msp_img,
                                                   mask_LR,
                                                   searchRadiusLR,
                                                   4 * this->m_TemplateRadius["PC"],
                                                   4 * this->m_TemplateRadius["PC"],
                                                   msp_lmk_RP_Candidate.GetVectorFromOrigin() + RPtoPC,
                                                   this->m_InputTemplateModel.GetTemplateMeans( "PC" ),
                                                   this->m_InputTemplateModel.m_VectorIndexLocations["PC"],
                                                   cc_PC_Max,
                                                   "PC" );
      }
      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
      // A check point for base landmarks
      if ( LMC::globalverboseFlag )
      {
        std::cout << cc_RP_Max << " " << cc_PC_Max << " " << cc_VN4_Max << " " << cc_AC_Max << "\n"
                  << msp_lmk_RP_Candidate << " " << msp_lmk_PC_Candiate << " " << msp_lmk_VN4_Candidate << " "
                  << msp_lmk_AC_Candidate << std::endl;
      }

      if ( globalImagedebugLevel > 1 )
      {
        std::string BrandedImageAName( this->m_ResultsDir + "/BrandedImage.png" );

        MakeBrandeddebugImage( this->m_msp_img.GetPointer(),
                               this->m_InputTemplateModel,
                               this->m_msp_lmk_CenterOfHeadMass + this->m_InputTemplateModel.GetCMtoRPMean(),
                               msp_lmk_RP_Candidate + RPtoAC,
                               msp_lmk_RP_Candidate + RPtoPC,
                               msp_lmk_RP_Candidate + RPtoVN4,
                               BrandedImageAName,
                               msp_lmk_RP_Candidate,
                               msp_lmk_AC_Candidate,
                               msp_lmk_PC_Candiate,
                               msp_lmk_VN4_Candidate );
        lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
        if ( globalImagedebugLevel > 3 )
        {
          std::string LabelImageAName( this->m_ResultsDir + "/MSP_Mask.nii.gz" );
          MakeLabelImage( this->m_msp_img,
                          msp_lmk_RP_Candidate,
                          msp_lmk_AC_Candidate,
                          msp_lmk_PC_Candiate,
                          msp_lmk_VN4_Candidate,
                          LabelImageAName );
        }
      }
      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
      // ============================================================================================
      //  Update landmark points


      // After finding RP, AC, PC, we can compute the versor transform for
      // registration
      // Also in this stage, we store some results for later use
      // Save named points in original space

      if ( !hasUserForcedRPPoint )
      {
        const SImagePointType tempLoc = this->m_test_orig2msp_img_tfm->TransformPoint( msp_lmk_RP_Candidate );

        this->m_orig_lmks_updated["RP"] = this->m_orig2eyeFixed_img_tfm->TransformPoint( tempLoc );
      }

      if ( !hasUserForcedVN4Point )
      {
        const SImagePointType tempLoc = this->m_test_orig2msp_img_tfm->TransformPoint( msp_lmk_VN4_Candidate );

        this->m_orig_lmks_updated["VN4"] = this->m_orig2eyeFixed_img_tfm->TransformPoint( tempLoc );
      }

      if ( !hasUserForcedACPoint )
      {
        const SImagePointType tempLoc = this->m_test_orig2msp_img_tfm->TransformPoint( msp_lmk_AC_Candidate );

        this->m_orig_lmks_updated["AC"] = this->m_orig2eyeFixed_img_tfm->TransformPoint( tempLoc );
        std::cout << "msp_lmk_AC:         " << msp_lmk_AC_Candidate << std::endl;
        std::cout << "tempLoc eyeFixed?:  " << tempLoc << std::endl;
        std::cout << "orig AC:            " << this->m_orig_lmks_updated["AC"] << std::endl;
      }

      if ( !hasUserForcedPCPoint )
      {
        const SImagePointType tempLoc = this->m_test_orig2msp_img_tfm->TransformPoint( msp_lmk_PC_Candiate );

        this->m_orig_lmks_updated["PC"] = this->m_orig2eyeFixed_img_tfm->TransformPoint( tempLoc );
      }
      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
      {
        // XXXXXXX
        // HACK TODO: "CM" seems to be wrong here

        this->m_orig_lmks_updated["CM"] =
          this->m_test_orig2msp_img_tfm->TransformPoint( this->m_msp_lmk_CenterOfHeadMass );
        this->m_orig_lmks_updated["CM"] =
          this->m_orig2eyeFixed_img_tfm->TransformPoint( this->m_orig_lmks_updated.at( "CM" ) );
      }
      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
      this->m_orig_lmks_updated["LE"] = this->m_orig_lmks_constant.at( "LE" );
      this->m_orig_lmks_updated["RE"] = this->m_orig_lmks_constant.at( "RE" );

      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );

      // Write some debug images
      {
        if ( globalImagedebugLevel > 3 )
        {
          std::string ResampledMaskmageName( this->m_ResultsDir + "/Resampled_Mask.nii.gz" );
          MakeLabelImage( this->m_eyeFixed_img,
                          msp_lmk_RP_Candidate,
                          msp_lmk_AC_Candidate,
                          msp_lmk_PC_Candiate,
                          msp_lmk_VN4_Candidate,
                          ResampledMaskmageName );

          std::string OrigMaskImageName( this->m_ResultsDir + "/Orig_Mask.nii.gz" );
          MakeLabelImage( this->m_eyeFixed_img,
                          this->m_orig_lmks_updated.at( "RP" ),
                          this->m_orig_lmks_updated.at( "AC" ),
                          this->m_orig_lmks_updated.at( "PC" ),
                          this->m_orig_lmks_updated.at( "VN4" ),
                          OrigMaskImageName );
        }
      }
      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
      const std::vector< std::string > base_lmk_names{ "RP", "AC", "PC", "VN4", "LE", "RE" };

      // Compute the AC-PC aligned transform
      // Note for the sake of EPCA, we need the transform at this stage
      // so that the method is robust against rotation
      this->m_orig2msp_img_tfm = this->Compute_orig2msp_img_tfm(
        this->m_orig_lmks_updated["RP"], this->m_orig_lmks_updated["AC"], this->m_orig_lmks_updated["PC"] );
      // NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)
      VersorTransformType::Pointer orig2msp_lmk_tfm =
        GetLandmarkTransformFromImageTransform( this->m_orig2msp_img_tfm.GetPointer() );


      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );

      // Save some named points in EMSP space mainly for debug use
      LandmarksMapType local_msp_lmks_algo_found; // named points in EMSP space
      {
        local_msp_lmks_algo_found["RP"] = msp_lmk_RP_Candidate;
        local_msp_lmks_algo_found["AC"] = msp_lmk_AC_Candidate;
        local_msp_lmks_algo_found["PC"] = msp_lmk_PC_Candiate;
        local_msp_lmks_algo_found["VN4"] = msp_lmk_VN4_Candidate;
        local_msp_lmks_algo_found["CM"] = this->m_msp_lmk_CenterOfHeadMass;
#if 1                                                  // OLD TEST PASSS
        local_msp_lmks_algo_found["LE"] = eyeFixed_LE; // TODO, We may have a better estimate of this now.
        local_msp_lmks_algo_found["RE"] = eyeFixed_RE; // TODO, We may have a better estimate of this now.
#else
        local_msp_lmks_algo_found["LE"] = orig2msp_lmk_tfm->TransformPoint( this->m_orig_lmk_LE );
        local_msp_lmks_algo_found["RE"] = orig2msp_lmk_tfm->TransformPoint( this->m_orig_lmk_RE );
#endif
      }

      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
      /*
       * For the rest of the landmarks
       *
       * The search center of the rest of the landmarks will be estimated
       * by linear model estimation with EPCA
       * Note The current version use landmarks in acpc-aligned space.
       */
      {
        // Get a copy of landmarks on ACPC plane for eliminating accumulative
        // errors of local search process

        // Save some named points in AC-PC aligned space
        LandmarksMapType msp_lmks_from_orig_space; // named points in ACPC-aligned space
        {
          for ( const auto & lmk_name : base_lmk_names )
          {
            msp_lmks_from_orig_space[lmk_name] =
              orig2msp_lmk_tfm->TransformPoint( this->m_orig_lmks_updated.at( lmk_name ) );
          }
        }
        lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
        // Build up an evolutionary processing list
        // order: RP, AC, PC, VN4, LE, RE, ...
        // Note this order should comply with the order we defined in LLS model
        std::vector< std::string > processingList = base_lmk_names;
        unsigned int               numBaseLandmarks = processingList.size();
        // Initialize the iteratively updated lmks (addinga new lmk each time through
        // the loop, from the initial 6 main lmks.
        LandmarksMapType iteratively_updated_msp_lmks = msp_lmks_from_orig_space; // Initialize with lmks from
                                                                                  // orig_space
        unsigned int dim = 3;
        for ( unsigned int LlsMat_idx = 1; LlsMat_idx <= m_LlsMatrices.size(); ++LlsMat_idx )
        {
          // The processing order is indicated by the length of EPCA coefficient
          auto Lls_info_pair = this->m_LlsMatrices.begin();
          while ( Lls_info_pair->second.columns() != ( numBaseLandmarks + LlsMat_idx - 2 ) * dim )
          {
            if ( Lls_info_pair != this->m_LlsMatrices.end() )
            {
              ++Lls_info_pair; // transversal
            }
            else
            {
              std::cerr << "Error: wrong number of parameters for linear model!" << std::endl;
              std::cerr << "Some EPCA landmarks may not be detected." << std::endl;
              return;
            }
          }
          lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
          const std::string LlsMatrix_name = Lls_info_pair->first;
          std::cout << "Processing iterative update: " << LlsMatrix_name << "..." << std::endl;
          {
            // in every iteration, the last name in the processing list
            // indicates the landmark to be estimated
            processingList.push_back( LlsMatrix_name );

            // Find search center by linear model estimation with
            // dimension reduction.
            // The result will be stored into m_msp_lmks[Lls_info_pair->first]
            LinearEstimation( iteratively_updated_msp_lmks, processingList, numBaseLandmarks );

            // check whether it is midline landmark, set search range
            // and modify search center accordingly
            double localSearchRadiusLR = this->m_SearchRadii[LlsMatrix_name];
            {
              std::vector< std::string >::const_iterator midlineIt = this->m_MidlinePointsList.begin();
              for ( ; midlineIt != this->m_MidlinePointsList.end(); ++midlineIt )
              {
                if ( ( *midlineIt ).compare( LlsMatrix_name ) == 0 )
                {
                  localSearchRadiusLR = searchRadiusLR;                 // a variable changed with err_MSP
                  iteratively_updated_msp_lmks[LlsMatrix_name][0] = 0.; // search starts near EMSP
                  break;
                }
              }
            }
            lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
            // TODO: Simplify below to msp_lmks_from_orig_space[Lls_info_pair->first] =
            // iteratively_updated_msp_lmks.at(Lls_info_pair->first)
            msp_lmks_from_orig_space[Lls_info_pair->first][0] = iteratively_updated_msp_lmks[Lls_info_pair->first][0];
            msp_lmks_from_orig_space[Lls_info_pair->first][1] = iteratively_updated_msp_lmks[Lls_info_pair->first][1];
            msp_lmks_from_orig_space[Lls_info_pair->first][2] = iteratively_updated_msp_lmks[Lls_info_pair->first][2];

            // Obtain the position of the current landmark in other spaces
            this->m_orig_lmks_updated[Lls_info_pair->first] =
              this->m_orig2msp_img_tfm->TransformPoint( msp_lmks_from_orig_space.at( Lls_info_pair->first ) );
            {
              local_msp_lmks_algo_found[Lls_info_pair->first] =
                eyeFixed2msp_lmk_tfm->TransformPoint( orig2eyeFixed_lmk_tfm->TransformPoint(
                  this->m_orig_lmks_updated.at( Lls_info_pair->first ) ) ); // TODO: Verify this
            }

            lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
            // Enable local search
            if ( 1 )
            {
              // local search
              double cc_Max = 0;
              local_msp_lmks_algo_found[Lls_info_pair->first] =
                FindCandidatePoints( this->m_msp_img,
                                     mask_LR,
                                     localSearchRadiusLR,
                                     this->m_SearchRadii[Lls_info_pair->first],
                                     this->m_SearchRadii[Lls_info_pair->first],
                                     local_msp_lmks_algo_found[Lls_info_pair->first].GetVectorFromOrigin(),
                                     this->m_InputTemplateModel.GetTemplateMeans( Lls_info_pair->first ),
                                     this->m_InputTemplateModel.m_VectorIndexLocations[Lls_info_pair->first],
                                     cc_Max,
                                     Lls_info_pair->first );

              // Update landmarks in input and ACPC-aligned space
              this->m_orig_lmks_updated[Lls_info_pair->first] =
                this->m_test_orig2msp_img_tfm->TransformPoint( local_msp_lmks_algo_found.at( Lls_info_pair->first ) );
              {
                this->m_orig_lmks_updated[Lls_info_pair->first] =
                  this->m_orig2eyeFixed_img_tfm->TransformPoint( this->m_orig_lmks_updated.at( Lls_info_pair->first ) );
              }
              msp_lmks_from_orig_space[Lls_info_pair->first] =
                orig2msp_lmk_tfm->TransformPoint( this->m_orig_lmks_updated.at( Lls_info_pair->first ) );
            }
            lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
          }
        } // End of arbitrary landmarks detection for the rest of "new" ones
      }   // End of arbitrary landmarks detection by linear model estimation

      this->ComputeFinalRefinedACPCAlignedTransform( orig_space_image, this->m_orig_lmks_updated );
      if ( globalImagedebugLevel > 1 ) // This may be very useful for GUI
      // corrector
      {
        WriteITKtoSlicer3Lmk( this->m_ResultsDir + "/EMSP.fcsv", local_msp_lmks_algo_found );
      }
      lmk_check_differences( m_orig_lmks_constant, m_orig_lmks_updated, false, __FILE__, __LINE__ );
    } // End of local searching kernel
  }   // End of local searching
}
