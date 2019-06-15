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
/*
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "landmarksConstellationDetector.h"
// landmarkIO has to be included after landmarksConstellationDetector
#include "landmarkIO.h"
#include "itkOrthogonalize3DRotationMatrix.h"

#include "itkFindCenterOfBrainFilter.h"
#include "BRAINSHoughEyeDetector.h"

#include <BRAINSFitHelper.h>
#include "itkLandmarkBasedTransformInitializer.h"

std::string
local_to_string( unsigned int i )
{
  std::stringstream localStream;
  localStream << i;
  return localStream.str();
}

// NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)

VersorTransformType::Pointer
landmarksConstellationDetector::GetLandmarkTransformFromImageTransform(
  VersorTransformType::ConstPointer orig2msp_img_tfm )
{
  VersorTransformType::Pointer orig2msp_lmk_tfm = VersorTransformType::New();
  SImageType::PointType        centerPoint = orig2msp_img_tfm->GetCenter();
  orig2msp_lmk_tfm->SetCenter( centerPoint );
  orig2msp_lmk_tfm->SetIdentity();
  orig2msp_img_tfm->GetInverse( orig2msp_lmk_tfm );
  return orig2msp_lmk_tfm;
}

VersorTransformType::Pointer
landmarksConstellationDetector::Compute_orig2msp_img_tfm( const SImagePointType & RP, const SImagePointType & AC,
                                                          const SImagePointType & PC )
{
  SImageType::PointType ZeroCenter;
  ZeroCenter.Fill( 0.0 );

  RigidTransformType::Pointer orig2msp_lmk_tfm_estimate = computeTmspFromPoints( RP, AC, PC, ZeroCenter );

  VersorTransformType::Pointer orig2msp_lmk_tfm_cleaned = VersorTransformType::New();
  orig2msp_lmk_tfm_cleaned->SetFixedParameters( orig2msp_lmk_tfm_estimate->GetFixedParameters() );

  itk::Versor< double >               versorRotation;
  const itk::Matrix< double, 3, 3 > & CleanedOrthogonalized =
    itk::Orthogonalize3DRotationMatrix( orig2msp_lmk_tfm_estimate->GetMatrix() );
  versorRotation.Set( CleanedOrthogonalized );

  orig2msp_lmk_tfm_cleaned->SetRotation( versorRotation );
  orig2msp_lmk_tfm_cleaned->SetTranslation( orig2msp_lmk_tfm_estimate->GetTranslation() );
  return orig2msp_lmk_tfm_cleaned;
}

void
landmarksConstellationDetector::ComputeFinalRefinedACPCAlignedTransform( SImageType::Pointer      original_space_image,
                                                                         const LandmarksMapType & updated_orig_lmks )
{
  ////////////////////////////
  // START BRAINSFit alternative
  if ( !this->m_atlasVolume.empty() )
  {
    using AtlasReaderType = itk::ImageFileReader< SImageType >;
    AtlasReaderType::Pointer atlasReader = AtlasReaderType::New();
    atlasReader->SetFileName( this->m_atlasVolume );
    std::cout << "read atlas: " << this->m_atlasVolume << std::endl;
    try
    {
      atlasReader->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
      std::cerr << "Error while reading atlasVolume file:\n " << err << std::endl;
    }

    // TODO: prob needs a try-catch
    std::cout << "read atlas landmarks:  " << this->m_atlasLandmarks << std::endl;
    LandmarksMapType referenceAtlasLandmarks = ReadSlicer3toITKLmk( this->m_atlasLandmarks );

    // Create a better version of this->m_orig2msp_img_tfm using BRAINSFit.
    // take the the subjects landmarks in original space, and  landmarks from a reference Atlas, and compute an initial
    // affine transform
    // ( using logic from BRAINSLandmarkInitializer) and create initToAtlasAffineTransform.

    using WeightType = std::map< std::string, float >;
    WeightType landmarkWeights;
    if ( this->m_atlasLandmarkWeights != "" )
    {
      landmarkWeights = ReadLandmarkWeights( this->m_atlasLandmarkWeights.c_str() );
      std::cout << "read atlas landmarksweights:  " << this->m_atlasLandmarkWeights << std::endl;
    }
    // TEST turning this back on.
    using LmkInitTransformType = itk::AffineTransform< double, Dimension >;
    using LandmarkBasedInitializerType =
      itk::LandmarkBasedTransformInitializer< LmkInitTransformType, SImageType, SImageType >;
    using LandmarkContainerType = LandmarkBasedInitializerType::LandmarkPointContainer;
    LandmarkContainerType atlasLmks;
    LandmarkContainerType movingLmks;
    using LandmarkConstIterator = LandmarksMapType::const_iterator;
    LandmarkBasedInitializerType::LandmarkWeightType landmarkWgts;
    for ( LandmarkConstIterator fixedIt = referenceAtlasLandmarks.begin(); fixedIt != referenceAtlasLandmarks.end();
          ++fixedIt )
    {
      LandmarkConstIterator movingIt = updated_orig_lmks.find( fixedIt->first );
      if ( movingIt != updated_orig_lmks.cend() )
      {
        atlasLmks.push_back( fixedIt->second );
        movingLmks.push_back( movingIt->second );
        if ( !this->m_atlasLandmarkWeights.empty() )
        {
          if ( landmarkWeights.find( fixedIt->first ) != landmarkWeights.end() )
          {
            landmarkWgts.push_back( landmarkWeights[fixedIt->first] );
          }
          else
          {
            std::cout << "Landmark for " << fixedIt->first << " does not exist. "
                      << "Set the weight to 0.5 " << std::endl;
            landmarkWgts.push_back( 0.5F );
          }
        }
      }
      else
      {
        itkGenericExceptionMacro( << "Could not find " << fixedIt->first << " in originalSpaceLandmarksPreBRAINSFit "
                                  << std::endl
                                  << "MIS MATCHED MOVING AND FIXED LANDMARKS!" << std::endl );
      }
    }

    LandmarkBasedInitializerType::Pointer landmarkBasedInitializer = LandmarkBasedInitializerType::New();

    if ( !this->m_atlasLandmarkWeights.empty() )
    {
      landmarkBasedInitializer->SetLandmarkWeight( landmarkWgts );
    }
    landmarkBasedInitializer->SetFixedLandmarks( atlasLmks );
    landmarkBasedInitializer->SetMovingLandmarks( movingLmks );

    LmkInitTransformType::Pointer initToAtlasAffineTransform = LmkInitTransformType::New();
    landmarkBasedInitializer->SetTransform( initToAtlasAffineTransform );
    landmarkBasedInitializer->InitializeTransform();

    using HelperType = itk::BRAINSFitHelper;
    HelperType::Pointer brainsFitHelper = HelperType::New();

    // Now Run BRAINSFitHelper class initialized with initToAtlasAffineTransform, original image, and atlas image
    // adapted from BRAINSABC/brainseg/AtlasRegistrationMethod.hxx - do I need to change any of these parameters?
    brainsFitHelper->SetSamplingPercentage( 0.05 ); // Use 5% of voxels for samples
    brainsFitHelper->SetNumberOfHistogramBins( 50 );
    const std::vector< int > numberOfIterations( 1, 1500 );
    brainsFitHelper->SetNumberOfIterations( numberOfIterations );
    brainsFitHelper->SetTranslationScale( 1000 );
    brainsFitHelper->SetReproportionScale( 1.0 );
    brainsFitHelper->SetSkewScale( 1.0 );

    using FloatImageType = itk::Image< float, 3 >;
    using CastFilterType = itk::CastImageFilter< SImageType, FloatImageType >;

    {
      CastFilterType::Pointer fixedCastFilter = CastFilterType::New();
      fixedCastFilter->SetInput( atlasReader->GetOutput() );
      fixedCastFilter->Update();
      brainsFitHelper->SetFixedVolume( fixedCastFilter->GetOutput() );

      CastFilterType::Pointer movingCastFilter = CastFilterType::New();
      movingCastFilter->SetInput( original_space_image );
      movingCastFilter->Update();
      brainsFitHelper->SetMovingVolume( movingCastFilter->GetOutput() );
      itkUtil::WriteImage< FloatImageType >( movingCastFilter->GetOutput(), "./DEBUGFloatCastMovingImage.nii.gz" );
    }

    const std::vector< double > minimumStepSize( 1, 0.0005 );
    brainsFitHelper->SetMinimumStepLength( minimumStepSize );
    std::vector< std::string > transformType( 1 );
    transformType[0] = "Affine";
    brainsFitHelper->SetTransformType( transformType );

    using CompositeTransformType = itk::CompositeTransform< double, 3 >;
    CompositeTransformType::Pointer initToAtlasAffineCompositeTransform =
      dynamic_cast< CompositeTransformType * >( initToAtlasAffineTransform.GetPointer() );
    if ( initToAtlasAffineCompositeTransform.IsNull() )
    {
      initToAtlasAffineCompositeTransform = CompositeTransformType::New();
      initToAtlasAffineCompositeTransform->AddTransform( initToAtlasAffineTransform );
    }
    brainsFitHelper->SetCurrentGenericTransform( initToAtlasAffineCompositeTransform );

    // Provide better masking of images
    {
      // if( fixedBinaryVolume != "" || movingBinaryVolume != "" )
      //  {
      //  std::cout
      //    << "ERROR:  Can not specify mask file names when ROIAUTO is used for the maskProcessingMode"
      //    << std::endl;
      //  return EXIT_FAILURE;
      //  }
      static constexpr unsigned int ROIAutoClosingSize = 4;
      static constexpr unsigned int ROIAutoDilateSize = 6;
      {
        using ROIAutoType = itk::BRAINSROIAutoImageFilter< SImageType, itk::Image< unsigned char, 3 > >;
        ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
        ROIFilter->SetInput( atlasReader->GetOutput() );
        ROIFilter->SetClosingSize( ROIAutoClosingSize );
        ROIFilter->SetDilateSize( ROIAutoDilateSize );
        ROIFilter->Update();
        ImageMaskPointer fixedMask = ROIFilter->GetSpatialObjectROI();
        brainsFitHelper->SetFixedBinaryVolume( fixedMask );
      }
      {
        using ROIAutoType = itk::BRAINSROIAutoImageFilter< SImageType, itk::Image< unsigned char, 3 > >;
        ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
        ROIFilter->SetInput( original_space_image );
        ROIFilter->SetClosingSize( ROIAutoClosingSize );
        ROIFilter->SetDilateSize( ROIAutoDilateSize );
        ROIFilter->Update();
        ImageMaskPointer movingMask = ROIFilter->GetSpatialObjectROI();
        brainsFitHelper->SetMovingBinaryVolume( movingMask );
      }
    }
    brainsFitHelper->SetDebugLevel( 10 );
    brainsFitHelper->Update();

    this->m_orig2msp_img_tfm = itk::ComputeRigidTransformFromGeneric(
      brainsFitHelper->GetCurrentGenericTransform()->GetNthTransform( 0 ).GetPointer() );
    if ( this->m_orig2msp_img_tfm.IsNull() )
    {
      // Fail if something weird happens.
      itkGenericExceptionMacro( << "this->m_orig2msp_img_tfm is null. "
                                << "It means we're not registering to the atlas, after all." << std::endl );
    }

    {
      // NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)
      using VersorRigid3DTransformType = itk::VersorRigid3DTransform< double >;
      VersorTransformType::Pointer LandmarkOrigToACPCTransform =
        GetLandmarkTransformFromImageTransform( this->m_orig2msp_img_tfm.GetPointer() );
      // TODO: Change name of LandmarkOrigToACPCTransform
      const VersorRigid3DTransformType::OutputPointType acPointInACPCSpace =
        LandmarkOrigToACPCTransform->TransformPoint( GetNamedPointFromLandmarkList( updated_orig_lmks, "AC" ) );
      // std::cout << "HACK: PRE-FIXING" << acPointInACPCSpace << std::endl;
      {
        VersorRigid3DTransformType::OffsetType translation;
        translation[0] =
          +acPointInACPCSpace[0]; // NOTE: Positive translation in ImageTransform is Negative in LandmarkTransform
        translation[1] = +acPointInACPCSpace[1];
        translation[2] = +acPointInACPCSpace[2];
        // First shift the transform
        this->m_orig2msp_img_tfm->Translate( translation, false );
      }

      // LandmarkOrigToACPCTransform = GetLandmarkTransformFromImageTransform( this->m_orig2msp_img_tfm.GetPointer()  );
      // VersorRigid3DTransformType::OutputPointType newacPointInACPCSpace =
      // LandmarkOrigToACPCTransform->TransformPoint(GetNamedPointFromLandmarkList(updated_orig_lmks,"AC")); std::cout
      // << "HACK: POST-FIXING" << newacPointInACPCSpace << std::endl;
      // TODO:  This still does not put it to (0,0,0) and it should.
    }
  }
  /// END BRAINSFIT_ALTERNATIVE
  ////////////////////////////
}

VersorTransformType::Pointer
landmarksConstellationDetector::GetImageOrigToACPCVersorTransform( void ) const
{
  return m_orig2msp_img_tfm;
}

SImageType::PointType
landmarksConstellationDetector::FindCandidatePoints(
  SImageType::Pointer volumeMSP, SImageType::Pointer mask_LR, const double LR_restrictions,
  const double PA_restrictions, const double SI_restrictions,
  // TODO: restrictions should really be ellipsoidal values
  const SImageType::PointType::VectorType &                      CenterOfSearchArea,
  const std::vector< std::vector< float > > &                    TemplateMean,
  const landmarksConstellationModelIO::IndexLocationVectorType & model, double & cc_Max, const std::string & mapID )
{
  cc_Max = -123456789.0;

  LinearInterpolatorType::Pointer imInterp = LinearInterpolatorType::New();
  imInterp->SetInputImage( volumeMSP );
  LinearInterpolatorType::Pointer maskInterp = LinearInterpolatorType::New();
  maskInterp->SetInputImage( mask_LR );

  // Final location is initialized with the center of search area
  SImageType::PointType GuessPoint;
  GuessPoint[0] = CenterOfSearchArea[0];
  GuessPoint[1] = CenterOfSearchArea[1];
  GuessPoint[2] = CenterOfSearchArea[2];

  // Boundary check
  {
    SImageType::PointType boundaryL( GuessPoint );
    boundaryL[0] += LR_restrictions;
    SImageType::PointType boundaryR( GuessPoint );
    boundaryR[0] -= LR_restrictions;
    SImageType::PointType boundaryP( GuessPoint );
    boundaryP[1] += PA_restrictions;
    SImageType::PointType boundaryA( GuessPoint );
    boundaryA[1] -= PA_restrictions;
    SImageType::PointType boundaryS( GuessPoint );
    boundaryS[2] += SI_restrictions;
    SImageType::PointType boundaryI( GuessPoint );
    boundaryI[2] -= SI_restrictions;
    if ( ( !maskInterp->IsInsideBuffer( boundaryL ) ) || ( !maskInterp->IsInsideBuffer( boundaryR ) ) ||
         ( !maskInterp->IsInsideBuffer( boundaryP ) ) || ( !maskInterp->IsInsideBuffer( boundaryA ) ) ||
         ( !maskInterp->IsInsideBuffer( boundaryS ) ) || ( !maskInterp->IsInsideBuffer( boundaryI ) ) )
    {
      std::cout << "WARNING: search region outside of the image region." << std::endl;
      std::cout << "The detection has probably large error!" << std::endl;
      return GuessPoint;
    }
  }

  // height and radius of the moving template
  const double height = this->m_InputTemplateModel.GetHeight( mapID );
  double       radii = this->m_InputTemplateModel.GetRadius( mapID );

  // HACK:
  // When, rmpj, rac, rpc and rvn4 are used, they help to increase the bounding box
  // restrictions; however, we do not want that landmark template image be affected
  // by that.
  if ( mapID == "RP" || mapID == "AC" || mapID == "PC" || mapID == "VN4" )
  {
    if ( radii > 5 )
    {
      radii = 5;
    }
  }

  // Bounding box around the center point.
  // To compute the correlation for the border points, the bounding box needs
  // to be expanded by the template size.
  //
  const double LeftToRight_BEGIN = CenterOfSearchArea[0] - LR_restrictions - height / 2;
  const double AnteriorToPostierior_BEGIN = CenterOfSearchArea[1] - PA_restrictions - radii;
  const double InferiorToSuperior_BEGIN = CenterOfSearchArea[2] - SI_restrictions - radii;

  const double LeftToRight_END = CenterOfSearchArea[0] + LR_restrictions + height / 2;
  const double AnteriorToPostierior_END = CenterOfSearchArea[1] + PA_restrictions + radii;
  const double InferiorToSuperior_END = CenterOfSearchArea[2] + SI_restrictions + radii;

  // Now bounding area will be converted to an image
  //
  SImageType::Pointer roiImage = SImageType::New();
  // origin
  SImageType::PointType roiOrigin;
  roiOrigin[0] = LeftToRight_BEGIN;
  roiOrigin[1] = AnteriorToPostierior_BEGIN;
  roiOrigin[2] = InferiorToSuperior_BEGIN;
  roiImage->SetOrigin( roiOrigin );
  // size
  SImageType::SizeType roiSize;
  roiSize[0] = ( LeftToRight_END - LeftToRight_BEGIN ) / ( volumeMSP->GetSpacing()[0] ) + 1;
  roiSize[1] = ( AnteriorToPostierior_END - AnteriorToPostierior_BEGIN ) / ( volumeMSP->GetSpacing()[1] ) + 1;
  roiSize[2] = ( InferiorToSuperior_END - InferiorToSuperior_BEGIN ) / ( volumeMSP->GetSpacing()[2] ) + 1;
  // start index
  SImageType::IndexType roiStart;
  roiStart.Fill( 0 );
  // region
  SImageType::RegionType roiRegion( roiStart, roiSize );
  roiImage->SetRegions( roiRegion );
  // spacing
  roiImage->SetSpacing( volumeMSP->GetSpacing() );
  // direction
  roiImage->SetDirection( volumeMSP->GetDirection() );
  roiImage->Allocate();
  roiImage->FillBuffer( 0 );

  // Since the actual bounding box is a rounded area, a Roi mask is also needed.
  //
  SImageType::Pointer roiMask = SImageType::New();
  roiMask->CopyInformation( roiImage );
  roiMask->SetRegions( roiImage->GetLargestPossibleRegion() );
  roiMask->Allocate();
  roiMask->FillBuffer( 0 );

  // roiImage is filled with values from volumeMSP
  //
  SImageType::PointType currentPointLocation;
  currentPointLocation[0] = CenterOfSearchArea[0];
  currentPointLocation[1] = CenterOfSearchArea[1];
  currentPointLocation[2] = CenterOfSearchArea[2];

  constexpr double deltaLR = 1; // in mm
  constexpr double deltaAP = 1; // in mm
  constexpr double deltaIS = 1; // in mm

  for ( double LeftToRight = LeftToRight_BEGIN; LeftToRight < LeftToRight_END; LeftToRight += deltaLR )
  {
    currentPointLocation[0] = LeftToRight;
    for ( double AnteriorToPostierior = AnteriorToPostierior_BEGIN; AnteriorToPostierior < AnteriorToPostierior_END;
          AnteriorToPostierior += deltaAP )
    {
      currentPointLocation[1] = AnteriorToPostierior;
      for ( double InferiorToSuperior = InferiorToSuperior_BEGIN; InferiorToSuperior < InferiorToSuperior_END;
            InferiorToSuperior += deltaIS )
      {
        currentPointLocation[2] = InferiorToSuperior;

        // Is current point within the input mask
        if ( maskInterp->Evaluate( currentPointLocation ) > 0.5 )
        {
          // Is current point inside the boundary box
          const SImageType::PointType::VectorType temp =
            currentPointLocation.GetVectorFromOrigin() - CenterOfSearchArea;
          const double inclusionDistance = temp.GetNorm();
          if ( ( inclusionDistance < ( SI_restrictions + radii ) ) &&
               ( std::abs( temp[1] ) < ( PA_restrictions + radii ) ) )
          {
            SImageType::IndexType index3D;
            roiImage->TransformPhysicalPointToIndex( currentPointLocation, index3D );
            roiImage->SetPixel( index3D, imInterp->Evaluate( currentPointLocation ) );
            roiMask->SetPixel( index3D, 1 );
          }
        }
      }
    }
  }

  ////////
  // Now we need to normalize only bounding region inside the roiImage
  ///////
  using BinaryImageToLabelMapFilterType = itk::BinaryImageToLabelMapFilter< SImageType >;
  BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
  binaryImageToLabelMapFilter->SetInput( roiMask );
  binaryImageToLabelMapFilter->Update();

  using LabelMapToLabelImageFilterType =
    itk::LabelMapToLabelImageFilter< BinaryImageToLabelMapFilterType::OutputImageType, SImageType >;
  LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
  labelMapToLabelImageFilter->SetInput( binaryImageToLabelMapFilter->GetOutput() );
  labelMapToLabelImageFilter->Update();

  using LabelStatisticsImageFilterType = itk::LabelStatisticsImageFilter< SImageType, SImageType >;
  LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetLabelInput( labelMapToLabelImageFilter->GetOutput() );
  labelStatisticsImageFilter->SetInput( roiImage );
  labelStatisticsImageFilter->Update();

  if ( labelStatisticsImageFilter->GetNumberOfLabels() != 1 )
  {
    itkGenericExceptionMacro( << "The bounding box mask should be connected." );
  }

  using LabelPixelType = LabelStatisticsImageFilterType::LabelPixelType;
  LabelPixelType      labelValue = labelStatisticsImageFilter->GetValidLabelValues()[0];
  const double        ROImean = labelStatisticsImageFilter->GetMean( labelValue );
  const double        ROIvar = labelStatisticsImageFilter->GetVariance( labelValue );
  const unsigned long ROIcount = labelStatisticsImageFilter->GetCount( labelValue );

  // The area inside the bounding box is normalized using the mean and variance statistics
  using SubtractImageFilterType = itk::SubtractImageFilter< SImageType, SImageType, FImageType3D >;
  SubtractImageFilterType::Pointer subtractConstantFromImageFilter = SubtractImageFilterType::New();
  subtractConstantFromImageFilter->SetInput( roiImage );
  subtractConstantFromImageFilter->SetConstant2( ROImean );
  subtractConstantFromImageFilter->Update();

  if ( std::sqrt( ROIcount * ROIvar ) < std::numeric_limits< double >::epsilon() )
  {
    itkGenericExceptionMacro( << "Zero norm for bounding area." );
  }
  const double normInv = 1 / ( std::sqrt( ROIcount * ROIvar ) );

  using MultiplyImageFilterType = itk::MultiplyImageFilter< FImageType3D, FImageType3D, FImageType3D >;
  MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput( subtractConstantFromImageFilter->GetOutput() );
  multiplyImageFilter->SetConstant( normInv );

  FImageType3D::Pointer normalizedRoiImage = multiplyImageFilter->GetOutput();
  /////////////// End of normalization of roiImage //////////////

  if ( globalImagedebugLevel > 8 )
  {
    std::string normalized_roiImage_name( this->m_ResultsDir + "/NormalizedRoiImage_" +
                                          itksys::SystemTools::GetFilenameName( mapID ) + ".nii.gz" );
    itkUtil::WriteImage< FImageType3D >( normalizedRoiImage, normalized_roiImage_name );

    std::string roiMask_name( this->m_ResultsDir + "/roiMask_" + itksys::SystemTools::GetFilenameName( mapID ) +
                              ".nii.gz" );
    itkUtil::WriteImage< SImageType >( roiMask, roiMask_name );
  }

  // Now each landmark template should be converted to a moving template image
  //
  FImageType3D::Pointer lmkTemplateImage = FImageType3D::New();
  lmkTemplateImage->SetOrigin( roiImage->GetOrigin() );
  lmkTemplateImage->SetSpacing( roiImage->GetSpacing() );
  lmkTemplateImage->SetDirection( roiImage->GetDirection() );
  FImageType3D::SizeType mi_size;
  mi_size[0] = 2 * height + 1;
  mi_size[1] = 2 * radii + 1;
  mi_size[2] = 2 * radii + 1;
  FImageType3D::IndexType mi_start;
  mi_start.Fill( 0 );
  FImageType3D::RegionType mi_region( mi_start, mi_size );
  lmkTemplateImage->SetRegions( mi_region );
  lmkTemplateImage->Allocate();

  // Since each landmark template is a cylinder, a template mask is needed.
  //
  SImageType::Pointer templateMask = SImageType::New();
  templateMask->CopyInformation( lmkTemplateImage );
  templateMask->SetRegions( lmkTemplateImage->GetLargestPossibleRegion() );
  templateMask->Allocate();

  // Fill the template moving image based on the vector index locations
  //  and template mean values for different angle rotations
  //
  double cc_rotation_max = 0.0;
  for ( unsigned int curr_rotationAngle = 0; curr_rotationAngle < TemplateMean.size(); curr_rotationAngle++ )
  {
    lmkTemplateImage->FillBuffer( 0 );
    templateMask->FillBuffer( 0 );
    // iterate over mean values for the current rotation angle
    std::vector< float >::const_iterator mean_iter = TemplateMean[curr_rotationAngle].begin();
    // Fill the lmk template image using the mean values
    for ( landmarksConstellationModelIO::IndexLocationVectorType::const_iterator it = model.begin(); it != model.end();
          ++it, ++mean_iter )
    {
      FImageType3D::IndexType pixelIndex;
      pixelIndex[0] = ( *it )[0] + height;
      pixelIndex[1] = ( *it )[1] + radii;
      pixelIndex[2] = ( *it )[2] + radii;
      lmkTemplateImage->SetPixel( pixelIndex, *mean_iter );
      templateMask->SetPixel( pixelIndex, 1 );
    }
    if ( globalImagedebugLevel > 8 )
    {
      std::string tmpImageName( this->m_ResultsDir + "/lmkTemplateImage_" +
                                itksys::SystemTools::GetFilenameName( mapID ) + "_" +
                                local_to_string( curr_rotationAngle ) + ".nii.gz" );
      itkUtil::WriteImage< FImageType3D >( lmkTemplateImage, tmpImageName );

      std::string tmpMaskName( this->m_ResultsDir + "/templateMask_" + itksys::SystemTools::GetFilenameName( mapID ) +
                               "_" + local_to_string( curr_rotationAngle ) + ".nii.gz" );
      itkUtil::WriteImage< SImageType >( templateMask, tmpMaskName );
    }

    // Finally NCC is calculated in frequency domain
    //
    using CorrelationFilterType =
      itk::MaskedFFTNormalizedCorrelationImageFilter< FImageType3D, FImageType3D, SImageType >;
    CorrelationFilterType::Pointer correlationFilter = CorrelationFilterType::New();
    correlationFilter->SetFixedImage( normalizedRoiImage );
    correlationFilter->SetFixedImageMask( roiMask );
    correlationFilter->SetMovingImage( lmkTemplateImage );
    correlationFilter->SetMovingImageMask( templateMask );
    correlationFilter->SetRequiredFractionOfOverlappingPixels( 1 );
    correlationFilter->Update();
    if ( globalImagedebugLevel > 8 )
    {
      std::string ncc_output_name( this->m_ResultsDir + "/NCCOutput_" + itksys::SystemTools::GetFilenameName( mapID ) +
                                   "_" + local_to_string( curr_rotationAngle ) + ".nii.gz" );
      itkUtil::WriteImage< FImageType3D >( correlationFilter->GetOutput(), ncc_output_name );
    }

    // Maximum NCC for current rotation angle
    using MinimumMaximumImageCalculatorType = itk::MinimumMaximumImageCalculator< FImageType3D >;
    MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter =
      MinimumMaximumImageCalculatorType::New();
    minimumMaximumImageCalculatorFilter->SetImage( correlationFilter->GetOutput() );
    minimumMaximumImageCalculatorFilter->Compute();
    double cc = minimumMaximumImageCalculatorFilter->GetMaximum();
    if ( cc > cc_rotation_max )
    {
      cc_rotation_max = cc;
      // Where maximum happens
      FImageType3D::IndexType maximumCorrelationPatchCenter = minimumMaximumImageCalculatorFilter->GetIndexOfMaximum();
      correlationFilter->GetOutput()->TransformIndexToPhysicalPoint( maximumCorrelationPatchCenter, GuessPoint );
    }
  }
  cc_Max = cc_rotation_max;

  if ( LMC::globalverboseFlag )
  {
    std::cout << "cc max: " << cc_Max << std::endl;
    std::cout << "guessed point in physical space: " << GuessPoint << std::endl;
  }
  return GuessPoint;
}

void
landmarksConstellationDetector::EulerToVersorRigid( VersorTransformType::Pointer &         result,
                                                    const RigidTransformType::ConstPointer eulerRigid )
{
  if ( result.IsNotNull() && eulerRigid.IsNotNull() )
  {
    result->SetFixedParameters( eulerRigid->GetFixedParameters() );
    itk::Versor< double >               versorRotation;
    const itk::Matrix< double, 3, 3 > & CleanedOrthogonalized =
      itk::Orthogonalize3DRotationMatrix( eulerRigid->GetMatrix() );
    versorRotation.Set( CleanedOrthogonalized );
    result->SetRotation( versorRotation );
    result->SetTranslation( eulerRigid->GetTranslation() );
  }
  else
  {
    itkGenericExceptionMacro( << "Error missing Pointer data, assigning "
                              << "Euler3DTransformPointer to VersorRigid3DTransformPointer." << std::endl );
  }
}


void
landmarksConstellationDetector::DoResampleInPlace( const SImageType::ConstPointer         inputImg,
                                                   const RigidTransformType::ConstPointer rigidTx,
                                                   SImageType::Pointer &                  inPlaceResampledImg )
{
  VersorTransformType::Pointer versorRigidTx = VersorTransformType::New();
  EulerToVersorRigid( versorRigidTx, rigidTx.GetPointer() );

  using ResampleIPFilterType = itk::ResampleInPlaceImageFilter< SImageType, SImageType >;
  using ResampleIPFilterPointer = ResampleIPFilterType::Pointer;
  ResampleIPFilterPointer InPlaceResampler = ResampleIPFilterType::New();
  InPlaceResampler->SetInputImage( inputImg );
  InPlaceResampler->SetRigidTransform( versorRigidTx.GetPointer() );
  InPlaceResampler->Update();

  inPlaceResampledImg = InPlaceResampler->GetOutput();
}

void
landmarksConstellationDetector::LinearEstimation( LandmarksMapType &                 msp_lmks_linearly_estimated,
                                                  const std::vector< std::string > & processingList,
                                                  unsigned                           numBasePoints )
{
  unsigned int dim = msp_lmks_linearly_estimated[processingList[0]].GetPointDimension();

  if ( processingList.size() <= numBasePoints )
  {
    std::cout << "No EPCA landmarks to be estimated." << std::endl;
    return;
  }
  std::string           newPointName = processingList[processingList.size() - 1];
  SImageType::PointType newPoint;
  newPoint.Fill( 0 );

  // Construct Xi_t
  VectorType Xi_t;
  Xi_t.set_size( dim * ( processingList.size() - 2 ) );
  for ( unsigned int k = 1; k <= processingList.size() - 2; ++k )
  {
    for ( unsigned int d = 0; d < dim; ++d )
    {
      Xi_t( ( k - 1 ) * dim + d ) = msp_lmks_linearly_estimated[processingList[k]][d] -
                                    msp_lmks_linearly_estimated[processingList[0]][d] -
                                    this->m_LlsMeans[newPointName][d];
    }
  }
  SImageType::PointType::VectorType tmp;
  tmp.SetVnlVector( this->m_LlsMatrices[newPointName] * Xi_t );
  newPoint = msp_lmks_linearly_estimated[processingList[0]] + tmp;
  msp_lmks_linearly_estimated[newPointName] = newPoint;

  // debug
  // std::cout << "Mi' = " << this->m_LlsMatrices[newPointName] << std::endl;
  // std::cout << "Xi_t = " << Xi_t << std::endl;
  // std::cout << "MPJ = " << msp_lmks_linearly_estimated[processingList[0]] << std::endl;
  // std::cout << newPointName << " = " << newPoint << std::endl;
}

SImageType::PointType::VectorType
landmarksConstellationDetector::FindVectorFromPointAndVectors( SImageType::PointType::VectorType BA,
                                                               SImageType::PointType::VectorType BAMean,
                                                               SImageType::PointType::VectorType BCMean, int sign )
{
  SImageType::PointType::VectorType BC;

  double cosTheta; // cosine of the angle from BA to BC

  cosTheta = BAMean * BCMean / BAMean.GetNorm() / BCMean.GetNorm();

  // Start searching on MSP
  BC[0] = 0;
  double a = BA * BA;
  double b = -2. * BA.GetNorm() * BCMean.GetNorm() * BA[2] * cosTheta;
  double c = BCMean * BCMean * ( BA * BA * cosTheta * cosTheta - BA[1] * BA[1] );
  double delta = b * b - 4 * a * c;
  if ( delta < 0 )
  {
    itkGenericExceptionMacro( << "Failed to solve a 2rd-order equation!" );
  }
  else if ( sign == 1 || sign == -1 )
  {
    BC[2] = -( b - sign * sqrt( delta ) ) / 2. / a;
  }
  else
  {
    itkGenericExceptionMacro( << "Bad parameter! sign = 1 or sign = -1 please" );
  }
  BC[1] = ( BA.GetNorm() * BCMean.GetNorm() * cosTheta - BA[2] * BC[2] ) / BA[1];
  return BC;
}


void
WriteManualFixFiles( const std::string & EMSP_Fiducial_file_name, SImageType * const mspVolume,
                     const std::string & resultDir, const LandmarksMapType & errorLmks,
                     const std::string & failureMessage, const bool throwException )
{ // ADD MetaData for EMSP_FCSV_FILENAME
  itk::MetaDataDictionary & dict = mspVolume->GetMetaDataDictionary();
  const char * const        metaDataEMSP_FCSVName = "EMSP_FCSV_FILENAME";
  itk::EncapsulateMetaData< std::string >( dict, metaDataEMSP_FCSVName, EMSP_Fiducial_file_name.c_str() );

  // write EMSP aligned image
  itkUtil::WriteImage< SImageType >( mspVolume, resultDir + "/EMSP.nrrd" );

  if ( errorLmks.size() > 0 )
  {
    WriteITKtoSlicer3Lmk( resultDir + "/" + EMSP_Fiducial_file_name, errorLmks );
  }
  if ( throwException )
  {
    itkGenericExceptionMacro( << failureMessage );
  }
}
