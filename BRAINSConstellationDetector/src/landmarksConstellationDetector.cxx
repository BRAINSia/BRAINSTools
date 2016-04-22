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

#include <sstream>

std::string local_to_string(unsigned int i)
{
  std::stringstream localStream;
  localStream << i;
  return localStream.str();
}

//NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)
static
VersorTransformType::Pointer
GetLandmarkTransformFromImageTransform( VersorTransformType::ConstPointer ImageTransform )
{
  VersorTransformType::Pointer LandmarkOrigToACPCTransform = VersorTransformType::New();
  SImageType::PointType        centerPoint = ImageTransform->GetCenter();
  LandmarkOrigToACPCTransform->SetCenter( centerPoint );
  LandmarkOrigToACPCTransform->SetIdentity();
  ImageTransform->GetInverse( LandmarkOrigToACPCTransform );
  return LandmarkOrigToACPCTransform;
}

VersorTransformType::Pointer
landmarksConstellationDetector::ComputeACPCAlignedZeroCenteredTransform(void)
{
  SImageType::PointType ZeroCenter;
  ZeroCenter.Fill(0.0);

  RigidTransformType::Pointer landmarkDefinedACPCAlignedToZeroTransform =
    computeTmspFromPoints(GetNamedPointFromLandmarkList(this->GetOriginalSpaceNamedPoints(),"RP"),
                          GetNamedPointFromLandmarkList(this->GetOriginalSpaceNamedPoints(),"AC"),
                          GetNamedPointFromLandmarkList(this->GetOriginalSpaceNamedPoints(),"PC"),
                          ZeroCenter);

  VersorTransformType::Pointer ACPCAlignedZeroCenteredTransform = VersorTransformType::New();
  ACPCAlignedZeroCenteredTransform->SetFixedParameters( landmarkDefinedACPCAlignedToZeroTransform->GetFixedParameters() );

  itk::Versor<double>               versorRotation;
  const itk::Matrix<double, 3, 3> & CleanedOrthogonalized =
    itk::Orthogonalize3DRotationMatrix( landmarkDefinedACPCAlignedToZeroTransform->GetMatrix() );
  versorRotation.Set( CleanedOrthogonalized );

  ACPCAlignedZeroCenteredTransform->SetRotation( versorRotation );
  ACPCAlignedZeroCenteredTransform->SetTranslation( landmarkDefinedACPCAlignedToZeroTransform->GetTranslation() );
  return ACPCAlignedZeroCenteredTransform;
}

void
landmarksConstellationDetector::ComputeFinalRefinedACPCAlignedTransform(void)
{
  ////////////////////////////
  // START BRAINSFit alternative
  if( ! this->m_atlasVolume.empty() )
    {
    typedef itk::ImageFileReader<SImageType> AtlasReaderType;
    AtlasReaderType::Pointer atlasReader = AtlasReaderType::New();
    atlasReader->SetFileName( this->m_atlasVolume );
    std::cout << "read atlas: " << this->m_atlasVolume << std::endl;
    try
      {
      atlasReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Error while reading atlasVolume file:\n "
        << err << std::endl;
      }

  // TODO: prob needs a try-catch
  std::cout << "read atlas landmarks:  " << this->m_atlasLandmarks << std::endl;
  LandmarksMapType referenceAtlasLandmarks = ReadSlicer3toITKLmk( this->m_atlasLandmarks );

  // Create a better version of this->m_ImageOrigToACPCVersorTransform using BRAINSFit.
  // take the the subjects landmarks in original space, and  landmarks from a reference Atlas, and compute an initial
  // affine transform
  // ( using logic from BRAINSLandmarkInitializer) and create initToAtlasAffineTransform.

  typedef std::map<std::string, float> WeightType;
  WeightType landmarkWeights;
  if( this->m_atlasLandmarkWeights != "" )
    {
    landmarkWeights = ReadLandmarkWeights( this->m_atlasLandmarkWeights.c_str() );
    std::cout << "read atlas landmarksweights:  " << this->m_atlasLandmarkWeights << std::endl;
    }
  // TEST turning this back on.
  typedef itk::AffineTransform<double, Dimension> LmkInitTransformType;
  typedef itk::LandmarkBasedTransformInitializer<LmkInitTransformType, SImageType, SImageType> LandmarkBasedInitializerType;
  typedef LandmarkBasedInitializerType::LandmarkPointContainer LandmarkContainerType;
  LandmarkContainerType atlasLmks;
  LandmarkContainerType movingLmks;
  typedef  LandmarksMapType::const_iterator LandmarkConstIterator;
  LandmarkBasedInitializerType::LandmarkWeightType landmarkWgts;
  for( LandmarkConstIterator fixedIt = referenceAtlasLandmarks.begin(); fixedIt != referenceAtlasLandmarks.end();
    ++fixedIt )
    {
    LandmarkConstIterator movingIt = this->m_OriginalSpaceNamedPoints.find( fixedIt->first );
    if( movingIt != this->m_OriginalSpaceNamedPoints.end() )
      {
      atlasLmks.push_back( fixedIt->second);
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
      itkGenericExceptionMacro(<< "Could not find " << fixedIt->first
                               << " in originalSpaceLandmarksPreBRAINSFit " << std::endl
                               << "MIS MATCHED MOVING AND FIXED LANDMARKS!" << std::endl);
      }
    }

  LandmarkBasedInitializerType::Pointer landmarkBasedInitializer = LandmarkBasedInitializerType::New();

  if( !this->m_atlasLandmarkWeights.empty() )
    {
    landmarkBasedInitializer->SetLandmarkWeight( landmarkWgts );
    }
  landmarkBasedInitializer->SetFixedLandmarks( atlasLmks );
  landmarkBasedInitializer->SetMovingLandmarks( movingLmks );

  LmkInitTransformType::Pointer initToAtlasAffineTransform = LmkInitTransformType::New();
  landmarkBasedInitializer->SetTransform( initToAtlasAffineTransform );
  landmarkBasedInitializer->InitializeTransform();

  typedef itk::BRAINSFitHelper HelperType;
  HelperType::Pointer brainsFitHelper = HelperType::New();

  // Now Run BRAINSFitHelper class initialized with initToAtlasAffineTransform, original image, and atlas image
  // adapted from BRAINSABC/brainseg/AtlasRegistrationMethod.hxx - do I need to change any of these parameters?
  brainsFitHelper->SetSamplingPercentage(0.05); // Use 5% of voxels for samples
  brainsFitHelper->SetNumberOfHistogramBins(50);
  const std::vector<int> numberOfIterations(1,1500);
  brainsFitHelper->SetNumberOfIterations(numberOfIterations);
  brainsFitHelper->SetTranslationScale(1000);
  brainsFitHelper->SetReproportionScale(1.0);
  brainsFitHelper->SetSkewScale(1.0);

  typedef itk::Image<float, 3>                            FloatImageType;
  typedef itk::CastImageFilter<SImageType, FloatImageType> CastFilterType;

    {
    CastFilterType::Pointer fixedCastFilter = CastFilterType::New();
    fixedCastFilter->SetInput( atlasReader->GetOutput() );
    fixedCastFilter->Update();
    brainsFitHelper->SetFixedVolume( fixedCastFilter->GetOutput() );

    CastFilterType::Pointer movingCastFilter = CastFilterType::New();
    movingCastFilter->SetInput( this->m_OriginalInputImage );
    movingCastFilter->Update();
    brainsFitHelper->SetMovingVolume( movingCastFilter->GetOutput() );
    itkUtil::WriteImage<FloatImageType>(movingCastFilter->GetOutput() , "./DEBUGFloatCastMovingImage.nii.gz");
    }

  const std::vector<double> minimumStepSize(1,0.0005);
  brainsFitHelper->SetMinimumStepLength(minimumStepSize);
  std::vector<std::string> transformType(1);
  transformType[0] = "Affine";
  brainsFitHelper->SetTransformType(transformType);

  typedef itk::CompositeTransform<double, 3> CompositeTransformType;
  CompositeTransformType::Pointer initToAtlasAffineCompositeTransform =
    dynamic_cast<CompositeTransformType *>( initToAtlasAffineTransform.GetPointer() );
  if( initToAtlasAffineCompositeTransform.IsNull() )
    {
    initToAtlasAffineCompositeTransform = CompositeTransformType::New();
    initToAtlasAffineCompositeTransform->AddTransform( initToAtlasAffineTransform );
    }
  brainsFitHelper->SetCurrentGenericTransform( initToAtlasAffineCompositeTransform );

  // Provide better masking of images
    {
    //if( fixedBinaryVolume != "" || movingBinaryVolume != "" )
    //  {
    //  std::cout
    //    << "ERROR:  Can not specify mask file names when ROIAUTO is used for the maskProcessingMode"
    //    << std::endl;
    //  return EXIT_FAILURE;
    //  }
    static const unsigned int ROIAutoClosingSize=4;
    static const unsigned int ROIAutoDilateSize=6;
      {
      typedef itk::BRAINSROIAutoImageFilter<SImageType, itk::Image<unsigned char, 3> > ROIAutoType;
      ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
      ROIFilter->SetInput( atlasReader->GetOutput() );
      ROIFilter->SetClosingSize(ROIAutoClosingSize);
      ROIFilter->SetDilateSize(ROIAutoDilateSize);
      ROIFilter->Update();
      ImageMaskPointer fixedMask = ROIFilter->GetSpatialObjectROI();
      brainsFitHelper->SetFixedBinaryVolume( fixedMask );
      }
      {
      typedef itk::BRAINSROIAutoImageFilter<SImageType, itk::Image<unsigned char, 3> > ROIAutoType;
      ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
      ROIFilter->SetInput( this->m_OriginalInputImage );
      ROIFilter->SetClosingSize(ROIAutoClosingSize);
      ROIFilter->SetDilateSize(ROIAutoDilateSize);
      ROIFilter->Update();
      ImageMaskPointer movingMask = ROIFilter->GetSpatialObjectROI();
      brainsFitHelper->SetMovingBinaryVolume( movingMask );
      }
    }
  brainsFitHelper->SetDebugLevel(10);
  brainsFitHelper->Update();

  this->m_ImageOrigToACPCVersorTransform = itk::ComputeRigidTransformFromGeneric( brainsFitHelper->GetCurrentGenericTransform()->GetNthTransform(0).GetPointer() );
  if( this->m_ImageOrigToACPCVersorTransform.IsNull() )
    {
    // Fail if something weird happens.
    itkGenericExceptionMacro(<< "this->m_ImageOrigToACPCVersorTransform is null. "
                             << "It means we're not registering to the atlas, after all." << std::endl);
    }

    {
    //NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)
    typedef itk::VersorRigid3DTransform<double>              VersorRigid3DTransformType;
    VersorTransformType::Pointer LandmarkOrigToACPCTransform = GetLandmarkTransformFromImageTransform( this->m_ImageOrigToACPCVersorTransform.GetPointer()  );

    const VersorRigid3DTransformType::OutputPointType acPointInACPCSpace = LandmarkOrigToACPCTransform->TransformPoint(GetNamedPointFromLandmarkList(this->m_OriginalSpaceNamedPoints,"AC"));
    //std::cout << "HACK: PRE-FIXING" << acPointInACPCSpace << std::endl;
      {
      VersorRigid3DTransformType::OffsetType translation;
      translation[0] = +acPointInACPCSpace[0]; //NOTE: Positive translation in ImageTransform is Negative in LandmarkTransform
      translation[1] = +acPointInACPCSpace[1];
      translation[2] = +acPointInACPCSpace[2];
      //First shift the transform
      this->m_ImageOrigToACPCVersorTransform->Translate( translation, false );
      }

    //LandmarkOrigToACPCTransform = GetLandmarkTransformFromImageTransform( this->m_ImageOrigToACPCVersorTransform.GetPointer()  );
    //VersorRigid3DTransformType::OutputPointType newacPointInACPCSpace = LandmarkOrigToACPCTransform->TransformPoint(GetNamedPointFromLandmarkList(this->m_OriginalSpaceNamedPoints,"AC"));
    //std::cout << "HACK: POST-FIXING" << newacPointInACPCSpace << std::endl;
    //TODO:  This still does not put it to (0,0,0) and it should.
    }
    }
  ///END BRAINSFIT_ALTERNATIVE
  ////////////////////////////
}

VersorTransformType::Pointer
landmarksConstellationDetector::GetImageOrigToACPCVersorTransform(void) const
{
  return m_ImageOrigToACPCVersorTransform;
}

SImageType::PointType
landmarksConstellationDetector::FindCandidatePoints
  ( SImageType::Pointer volumeMSP,
  SImageType::Pointer mask_LR,
  const double LR_restrictions,
  const double PA_restrictions,
  const double SI_restrictions,
  // TODO: restrictions should really be ellipsoidal values
  const SImageType::PointType::VectorType & CenterOfSearchArea,
  const std::vector<std::vector<float> > & TemplateMean,
  const landmarksConstellationModelIO::IndexLocationVectorType & model,
  double & cc_Max, const std::string & mapID )
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
  if( ( !maskInterp->IsInsideBuffer( boundaryL ) ) ||
     ( !maskInterp->IsInsideBuffer( boundaryR ) ) ||
     ( !maskInterp->IsInsideBuffer( boundaryP ) ) ||
     ( !maskInterp->IsInsideBuffer( boundaryA ) ) ||
     ( !maskInterp->IsInsideBuffer( boundaryS ) ) ||
     ( !maskInterp->IsInsideBuffer( boundaryI ) ) )
    {
    std::cout << "WARNING: search region outside of the image region." << std::endl;
    std::cout << "The detection has probably large error!" << std::endl;
    return GuessPoint;
    }
  }

  // height and radius of the moving template
  const double height = this->m_InputTemplateModel.GetHeight(mapID);
  double radii = this->m_InputTemplateModel.GetRadius(mapID);

  // HACK:
  // When, rmpj, rac, rpc and rvn4 are used, they help to increase the bounding box
  // restrictions; however, we do not want that landmark template image be affected
  // by that.
  if( mapID == "RP" || mapID == "AC" || mapID == "PC" || mapID == "VN4")
    {
    if( radii > 5 )
      {
      radii = 5;
      }
    }

  // Bounding box around the center point.
  // To compute the correlation for the border points, the bounding box needs
  // to be expanded by the template size.
  //
  const double LeftToRight_BEGIN = CenterOfSearchArea[0] - LR_restrictions - height/2;
  const double AnteriorToPostierior_BEGIN = CenterOfSearchArea[1] - PA_restrictions - radii;
  const double InferiorToSuperior_BEGIN = CenterOfSearchArea[2] - SI_restrictions - radii;

  const double LeftToRight_END = CenterOfSearchArea[0] + LR_restrictions + height/2;
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
  roiSize[0] = (LeftToRight_END - LeftToRight_BEGIN)/(volumeMSP->GetSpacing()[0]) + 1;
  roiSize[1] = (AnteriorToPostierior_END - AnteriorToPostierior_BEGIN)/(volumeMSP->GetSpacing()[1]) + 1;
  roiSize[2] = (InferiorToSuperior_END - InferiorToSuperior_BEGIN)/(volumeMSP->GetSpacing()[2]) + 1;
  // start index
  SImageType::IndexType roiStart;
  roiStart.Fill(0);
  // region
  SImageType::RegionType roiRegion(roiStart,roiSize);
  roiImage->SetRegions(roiRegion);
  // spacing
  roiImage->SetSpacing( volumeMSP->GetSpacing() );
  // direction
  roiImage->SetDirection( volumeMSP->GetDirection() );
  roiImage->Allocate();
  roiImage->FillBuffer(0);

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

  const double deltaLR = 1; // in mm
  const double deltaAP = 1; // in mm
  const double deltaIS = 1; // in mm

  for( double LeftToRight = LeftToRight_BEGIN;
      LeftToRight < LeftToRight_END; LeftToRight += deltaLR )
    {
    currentPointLocation[0] = LeftToRight;
    for( double AnteriorToPostierior = AnteriorToPostierior_BEGIN;
        AnteriorToPostierior < AnteriorToPostierior_END; AnteriorToPostierior += deltaAP )
      {
      currentPointLocation[1] = AnteriorToPostierior;
      for( double InferiorToSuperior = InferiorToSuperior_BEGIN;
          InferiorToSuperior < InferiorToSuperior_END; InferiorToSuperior += deltaIS )
        {
        currentPointLocation[2] = InferiorToSuperior;

          // Is current point within the input mask
        if( maskInterp->Evaluate( currentPointLocation ) > 0.5 )
          {
            // Is current point inside the boundary box
          const SImageType::PointType::VectorType temp =
                                          currentPointLocation.GetVectorFromOrigin() - CenterOfSearchArea;
          const double inclusionDistance = temp.GetNorm();
          if( ( inclusionDistance < (SI_restrictions+radii) ) && ( std::abs( temp[1] ) < (PA_restrictions+radii) ) )
            {
            SImageType::IndexType index3D;
            roiImage->TransformPhysicalPointToIndex( currentPointLocation, index3D );
            roiImage->SetPixel(index3D,imInterp->Evaluate(currentPointLocation));
            roiMask->SetPixel(index3D,1);
            }
          }
        }
      }
    }

  ////////
  // Now we need to normalize only bounding region inside the roiImage
  ///////
  typedef itk::BinaryImageToLabelMapFilter<SImageType> BinaryImageToLabelMapFilterType;
  BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
  binaryImageToLabelMapFilter->SetInput( roiMask );
  binaryImageToLabelMapFilter->Update();

  typedef itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, SImageType> LabelMapToLabelImageFilterType;
  LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
  labelMapToLabelImageFilter->SetInput( binaryImageToLabelMapFilter->GetOutput() );
  labelMapToLabelImageFilter->Update();

  typedef itk::LabelStatisticsImageFilter< SImageType, SImageType > LabelStatisticsImageFilterType;
  LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetLabelInput( labelMapToLabelImageFilter->GetOutput() );
  labelStatisticsImageFilter->SetInput( roiImage );
  labelStatisticsImageFilter->Update();

  if( labelStatisticsImageFilter->GetNumberOfLabels() != 1 )
    {
    itkGenericExceptionMacro(<< "The bounding box mask should be connected.");
    }

  typedef LabelStatisticsImageFilterType::LabelPixelType                LabelPixelType;
  LabelPixelType labelValue = labelStatisticsImageFilter->GetValidLabelValues()[0];
  const double ROImean = labelStatisticsImageFilter->GetMean( labelValue );
  const double ROIvar = labelStatisticsImageFilter->GetVariance( labelValue );
  const unsigned long ROIcount = labelStatisticsImageFilter->GetCount( labelValue );

  // The area inside the bounding box is normalized using the mean and variance statistics
  typedef itk::SubtractImageFilter <SImageType, SImageType, FImageType3D> SubtractImageFilterType;
  SubtractImageFilterType::Pointer subtractConstantFromImageFilter = SubtractImageFilterType::New();
  subtractConstantFromImageFilter->SetInput( roiImage );
  subtractConstantFromImageFilter->SetConstant2( ROImean );
  subtractConstantFromImageFilter->Update();

  if( std::sqrt( ROIcount * ROIvar ) < std::numeric_limits<double>::epsilon() )
    {
    itkGenericExceptionMacro(<< "Zero norm for bounding area.");
    }
  const double normInv = 1 / (std::sqrt( ROIcount * ROIvar ));

  typedef itk::MultiplyImageFilter<FImageType3D, FImageType3D, FImageType3D> MultiplyImageFilterType;
  MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput( subtractConstantFromImageFilter->GetOutput() );
  multiplyImageFilter->SetConstant( normInv );

  FImageType3D::Pointer normalizedRoiImage = multiplyImageFilter->GetOutput();
  /////////////// End of normalization of roiImage //////////////

  if( globalImagedebugLevel > 8 )
    {
    std::string normalized_roiImage_name( this->m_ResultsDir + "/NormalizedRoiImage_"
                                        + itksys::SystemTools::GetFilenameName( mapID )
                                        + ".nii.gz" );
    itkUtil::WriteImage<FImageType3D>( normalizedRoiImage, normalized_roiImage_name );

    std::string roiMask_name( this->m_ResultsDir + "/roiMask_"
                          + itksys::SystemTools::GetFilenameName( mapID )
                          + ".nii.gz" );
    itkUtil::WriteImage<SImageType>( roiMask, roiMask_name );
    }

  // Now each landmark template should be converted to a moving template image
  //
  FImageType3D::Pointer lmkTemplateImage = FImageType3D::New();
  lmkTemplateImage->SetOrigin( roiImage->GetOrigin() );
  lmkTemplateImage->SetSpacing( roiImage->GetSpacing() );
  lmkTemplateImage->SetDirection( roiImage->GetDirection() );
  FImageType3D::SizeType mi_size;
  mi_size[0] = 2*height+1;
  mi_size[1] = 2*radii+1;
  mi_size[2] = 2*radii+1;
  FImageType3D::IndexType mi_start;
  mi_start.Fill(0);
  FImageType3D::RegionType mi_region(mi_start,mi_size);
  lmkTemplateImage->SetRegions(mi_region);
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
  for( unsigned int curr_rotationAngle = 0;
      curr_rotationAngle < TemplateMean.size(); curr_rotationAngle++ )
    {
    lmkTemplateImage->FillBuffer(0);
    templateMask->FillBuffer( 0 );
    // iterate over mean values for the current rotation angle
    std::vector<float>::const_iterator mean_iter = TemplateMean[curr_rotationAngle].begin();
    // Fill the lmk template image using the mean values
    for( landmarksConstellationModelIO::IndexLocationVectorType::const_iterator it = model.begin();
        it != model.end();
        ++it, ++mean_iter )
      {
      FImageType3D::IndexType pixelIndex;
      pixelIndex[0] = (*it)[0]+height;
      pixelIndex[1] = (*it)[1]+radii;
      pixelIndex[2] = (*it)[2]+radii;
      lmkTemplateImage->SetPixel( pixelIndex, *mean_iter );
      templateMask->SetPixel( pixelIndex, 1 );
      }
    if( globalImagedebugLevel > 8 )
      {
      std::string tmpImageName( this->m_ResultsDir + "/lmkTemplateImage_"
                          + itksys::SystemTools::GetFilenameName( mapID ) + "_"
                          + local_to_string(curr_rotationAngle) + ".nii.gz" );
      itkUtil::WriteImage<FImageType3D>( lmkTemplateImage, tmpImageName );

      std::string tmpMaskName( this->m_ResultsDir + "/templateMask_"
                          + itksys::SystemTools::GetFilenameName( mapID ) + "_"
                          + local_to_string(curr_rotationAngle) + ".nii.gz" );
      itkUtil::WriteImage<SImageType>( templateMask, tmpMaskName );
      }

    // Finally NCC is calculated in frequency domain
    //
    typedef itk::MaskedFFTNormalizedCorrelationImageFilter<FImageType3D, FImageType3D, SImageType> CorrelationFilterType;
    CorrelationFilterType::Pointer correlationFilter = CorrelationFilterType::New();
    correlationFilter->SetFixedImage( normalizedRoiImage );
    correlationFilter->SetFixedImageMask( roiMask );
    correlationFilter->SetMovingImage( lmkTemplateImage );
    correlationFilter->SetMovingImageMask( templateMask );
    correlationFilter->SetRequiredFractionOfOverlappingPixels( 1 );
    correlationFilter->Update();
    if( globalImagedebugLevel > 8 )
      {
      std::string ncc_output_name( this->m_ResultsDir + "/NCCOutput_"
                          + itksys::SystemTools::GetFilenameName( mapID ) + "_"
                          + local_to_string(curr_rotationAngle) + ".nii.gz" );
      itkUtil::WriteImage<FImageType3D>( correlationFilter->GetOutput(), ncc_output_name );
      }

    // Maximum NCC for current rotation angle
    typedef itk::MinimumMaximumImageCalculator<FImageType3D> MinimumMaximumImageCalculatorType;
    MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter =
    MinimumMaximumImageCalculatorType::New ();
    minimumMaximumImageCalculatorFilter->SetImage( correlationFilter->GetOutput() );
    minimumMaximumImageCalculatorFilter->Compute();
    double cc = minimumMaximumImageCalculatorFilter->GetMaximum();
    if( cc > cc_rotation_max )
      {
      cc_rotation_max = cc;
      // Where maximum happens
      FImageType3D::IndexType maximumCorrelationPatchCenter =
                                      minimumMaximumImageCalculatorFilter->GetIndexOfMaximum();
      correlationFilter->GetOutput()->TransformIndexToPhysicalPoint( maximumCorrelationPatchCenter, GuessPoint );
      }
    }
  cc_Max = cc_rotation_max;

  if( LMC::globalverboseFlag )
    {
    std::cout << "cc max: " << cc_Max << std::endl;
    std::cout << "guessed point in physical space: " << GuessPoint << std::endl;
    }
  return GuessPoint;
}

void
landmarksConstellationDetector::EulerToVersorRigid( VersorTransformType::Pointer & result,
                                                    const RigidTransformType::ConstPointer eulerRigid )
{
  if( result.IsNotNull() && eulerRigid.IsNotNull() )
    {
    result->SetFixedParameters( eulerRigid->GetFixedParameters() );
    itk::Versor<double>               versorRotation;
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized = itk::Orthogonalize3DRotationMatrix( eulerRigid->GetMatrix() );
    versorRotation.Set( CleanedOrthogonalized );
    result->SetRotation(versorRotation);
    result->SetTranslation( eulerRigid->GetTranslation() );
    }
  else
    {
    itkGenericExceptionMacro(<< "Error missing Pointer data, assigning "
                             << "Euler3DTransformPointer to VersorRigid3DTransformPointer." << std::endl);
    }
}

void landmarksConstellationDetector::DoResampleInPlace( const SImageType::ConstPointer inputImg,
                                                        const RigidTransformType::ConstPointer rigidTx,
                                                        SImageType::Pointer & inPlaceResampledImg )
{
  VersorTransformType::Pointer versorRigidTx = VersorTransformType::New();
  EulerToVersorRigid( versorRigidTx, rigidTx.GetPointer() );

  typedef itk::ResampleInPlaceImageFilter<SImageType, SImageType> ResampleIPFilterType;
  typedef ResampleIPFilterType::Pointer                           ResampleIPFilterPointer;
  ResampleIPFilterPointer InPlaceResampler = ResampleIPFilterType::New();
  InPlaceResampler->SetInputImage( inputImg );
  InPlaceResampler->SetRigidTransform( versorRigidTx.GetPointer() );
  InPlaceResampler->Update();

  inPlaceResampledImg = InPlaceResampler->GetOutput();
}

void landmarksConstellationDetector::Compute( void )
{
  std::cout << "\nEstimating MSP..." << std::endl;

  // save the result that whether we are going to process all the landmarks
  // in light of user-specified eye center info.
  bool hasUserSpecEyeCenterInfo = true;
  bool hasUserForcedRPPoint = true;
  bool hasUserForcedACPoint = true;
  bool hasUserForcedPCPoint = true;
  bool hasUserForcedVN4Point = true;

  if( ( this->m_NamedPointEMSP.find( "LE" ) == this->m_NamedPointEMSP.end() )
      || ( this->m_NamedPointEMSP.find( "RE" ) == this->m_NamedPointEMSP.end() ) )
    {
    hasUserSpecEyeCenterInfo = false;
    }

  if( this->m_OriginalSpaceNamedPoints.find("RP") == this->m_OriginalSpaceNamedPoints.end() )
    {
    hasUserForcedRPPoint = false;
    }

  if( this->m_OriginalSpaceNamedPoints.find("AC") == this->m_OriginalSpaceNamedPoints.end() )
    {
    hasUserForcedACPoint = false;
    }

  if( this->m_OriginalSpaceNamedPoints.find("PC") == this->m_OriginalSpaceNamedPoints.end() )
    {
    hasUserForcedPCPoint = false;
    }

  if( this->m_OriginalSpaceNamedPoints.find("VN4") == this->m_OriginalSpaceNamedPoints.end() )
    {
    hasUserForcedVN4Point = false;
    }

  if( globalImagedebugLevel > 2 )
    {
    LandmarksMapType roughlyAlignedCHM;
    roughlyAlignedCHM["CM"] = this->m_CenterOfHeadMass;
    const std::string roughlyAlignedCHMName( this->m_ResultsDir + "/roughlyAlignedCHM.fcsv" );
    WriteITKtoSlicer3Lmk( roughlyAlignedCHMName, roughlyAlignedCHM );

    const std::string roughlyAlignedVolumeName( this->m_ResultsDir + "/VolumeRoughAlignedWithHoughEye.nrrd" );
    itkUtil::WriteImage<SImageType>( this->m_VolumeRoughAlignedWithHoughEye, roughlyAlignedVolumeName );
    }

  // Compute the estimated MSP transform, and aligned image
  double c_c = 0;
  ComputeMSP( this->m_VolumeRoughAlignedWithHoughEye, this->m_finalTmsp,
              this->m_VolumeMSP, this->m_CenterOfHeadMass, this->m_mspQualityLevel, c_c );
#if 0
  /*
   * If the MSP estimation is not good enough (i.e. c_c < -0.64), we try to compute the MSP again
   * using the previous m_finalTmsp as an initial transform.
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
        DoResampleInPlace( this->m_VolumeRoughAlignedWithHoughEye.GetPointer(), // input
                           this->m_finalTmsp.GetPointer(), // input
                           localRoughAlignedInput ); // output

        // Transform the centerOfHeadMass by finlaTmsp transform.
        VersorTransformType::Pointer localInvFinalTmsp = VersorTransformType::New();
        this->m_finalTmsp->GetInverse( localInvFinalTmsp );
        SImageType::PointType localAlignedCOHM = localInvFinalTmsp->TransformPoint( this->m_CenterOfHeadMass );

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
                   this->m_VolumeMSP, localAlignedCOHM, this->m_mspQualityLevel, c_c );

        // Compose the Tmsp transforms
        this->m_finalTmsp->Compose(localTmsp);

        if ( c_c < -0.64 )
          {
          break;
          }
        }

      typedef itk::StatisticsImageFilter<SImageType> StatisticsFilterType;
      StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
      statisticsFilter->SetInput(this->m_VolumeRoughAlignedWithHoughEye);
      statisticsFilter->Update();
      SImageType::PixelType minPixelValue = statisticsFilter->GetMinimum();

      this->m_VolumeMSP = TransformResample<SImageType, SImageType>( this->m_VolumeRoughAlignedWithHoughEye.GetPointer(),
                                                                    MakeIsoTropicReferenceImage().GetPointer(),
                                                                    minPixelValue,
                                                                    GetInterpolatorFromString<SImageType>("Linear").GetPointer(),
                                                                    this->GetTransformToMSP().GetPointer() );
    }
#endif
  std::cout << "\n=============================================================" << std::endl;

  // Generate a warning if reflective correlation similarity measure is low.
  // It may be normal in some very diseased subjects, so don't throw an exception here.
  if( c_c > -0.64 )
    {
    std::cout << "WARNING: Low reflective correlation between left/right hemispheres." << std::endl
    << "The estimated landmarks may not be reliable.\n" << std::endl;
    }

  // Throw an exception and stop BCD if RC metric is too low (less than 0.4) because results will not be reliable.
  if( c_c > -0.4 )
    {
    itkGenericExceptionMacro(<< "Too large MSP estimation error! reflective correlation metric is: "
                             << c_c << std::endl
                             << "Estimation of landmarks will not be reliable.\n" << std::endl);
    }

  // In case hough eye detector failed
  if( this->m_HoughEyeFailure || ( globalImagedebugLevel > 1 ) )
    {
    const std::string EMSP_Fiducial_file_name("EMSP.fcsv");
    //ADD MetaData for EMSP_FCSV_FILENAME
    itk::MetaDataDictionary &dict = this->m_VolumeMSP->GetMetaDataDictionary();
    const char * const metaDataEMSP_FCSVName = "EMSP_FCSV_FILENAME";
    itk::EncapsulateMetaData<std::string>(dict,metaDataEMSP_FCSVName,EMSP_Fiducial_file_name.c_str());

    // write EMSP aligned image
    itkUtil::WriteImage<SImageType> ( this->m_VolumeMSP, this->m_ResultsDir + "/EMSP.nrrd" );

    if( this->m_HoughEyeFailure )
      {
      SImageType::PointType zeroPoint;
      zeroPoint.Fill( 0 );
      LandmarksMapType zeroEyeCenters;
      zeroEyeCenters["LE"] = zeroPoint;
      zeroEyeCenters["RE"] = zeroPoint;
      WriteITKtoSlicer3Lmk
        ( this->m_ResultsDir + "/"+EMSP_Fiducial_file_name, zeroEyeCenters );
      itkGenericExceptionMacro(<< "EMSP aligned image and zero eye centers "
                               << "landmarks are written to " << std::endl
                               << this->m_ResultsDir << ". Use GUI corrector to "
                               << "correct the landmarks in." << EMSP_Fiducial_file_name);
      }
    }

  if( globalImagedebugLevel > 2 )
    {
    const std::string MSP_ImagePlaneForOrig( this->m_ResultsDir + "/MSP_PLANE_For_Orig.nii.gz" );
    CreatedebugPlaneImage( this->m_VolumeRoughAlignedWithHoughEye, MSP_ImagePlaneForOrig );
    const std::string MSP_ImagePlane( this->m_ResultsDir + "/MSP_PLANE.nii.gz" );
    CreatedebugPlaneImage( this->m_VolumeMSP, MSP_ImagePlane );
    }

  // HACK:  TODO:  DEBUG: Need to remove redundant need for volume and mask when
  // they can be the same image ( perhaps with a threshold );
  SImageType::Pointer mask_LR = this->m_VolumeMSP;

  VersorTransformType::Pointer InvHoughEyeTransform;
  if( !hasUserSpecEyeCenterInfo )
    {
    InvHoughEyeTransform = VersorTransformType::New();
    this->m_HoughEyeTransform->GetInverse( InvHoughEyeTransform );
    }

  VersorTransformType::Pointer InvFinalTmsp = VersorTransformType::New();
  this->m_finalTmsp->GetInverse( InvFinalTmsp );

  this->m_CenterOfHeadMassEMSP =
    InvFinalTmsp->TransformPoint( this->m_CenterOfHeadMass );
  this->m_CenterOfHeadMassEMSP[0] = 0; // Search starts on the estimated MSP

    {
    SImageType::PointType CandidateRPPoint;
    SImageType::PointType CandidatePCPoint;
    SImageType::PointType CandidateVN4Point;
    SImageType::PointType CandidateACPoint;

    SImageType::PointType::VectorType RPtoCEC;
    SImageType::PointType::VectorType RPtoVN4;
    SImageType::PointType::VectorType RPtoAC;
    SImageType::PointType::VectorType RPtoPC;

    SImageType::PointType mspSpaceCEC;

    if( !hasUserSpecEyeCenterInfo )
      {
      m_NamedPointEMSP["LE"] =
        InvHoughEyeTransform->TransformPoint( this->m_LEPoint );
      m_NamedPointEMSP["LE"] =
        InvFinalTmsp->TransformPoint( this->m_NamedPointEMSP["LE"] );
      m_NamedPointEMSP["RE"] =
        InvHoughEyeTransform->TransformPoint( this->m_REPoint );
      m_NamedPointEMSP["RE"] =
        InvFinalTmsp->TransformPoint( this->m_NamedPointEMSP["RE"] );
      }

    mspSpaceCEC.SetToMidPoint( this->m_NamedPointEMSP["LE"],
                               this->m_NamedPointEMSP["RE"] );
    mspSpaceCEC[0] = 0; // Search starts on the estimated MSP

    std::cout << "\nPerforming morphmetric search + local search..." << std::endl;
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
      if( this->m_NamedPointEMSP.find( "RP" ) != this->m_NamedPointEMSP.end() )
        {
        CandidateRPPoint = this->m_NamedPointEMSP["RP"];
        std::cout << "Skip estimation, directly load from file." << std::endl;
        }
      else if( hasUserForcedRPPoint )
        {
        std::cout << "Skip estimation, directly forced by command line." << std::endl;
        // Points chosen by Slicer should be converted to proper ITK representation
        this->m_OriginalSpaceNamedPoints["RP"][0] = -(this->m_OriginalSpaceNamedPoints["RP"][0]);
        this->m_OriginalSpaceNamedPoints["RP"][1] = -(this->m_OriginalSpaceNamedPoints["RP"][1]);
        CandidateRPPoint = InvFinalTmsp->TransformPoint( InvHoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["RP"] ) );
        }
      else
        {
        // The search radius of RP is set to 5 times larger than its template
        // radius in SI direction.
        // It's large because some scans have extra contents of neck or shoulder
        CandidateRPPoint =
          FindCandidatePoints( this->m_VolumeMSP, mask_LR, searchRadiusLR,
                               3. * this->m_TemplateRadius["RP"],
                               5. * this->m_TemplateRadius["RP"],
                               this->m_CenterOfHeadMassEMSP.GetVectorFromOrigin()
                               + this->m_InputTemplateModel.GetCMtoRPMean(),
                               this->m_InputTemplateModel.GetTemplateMeans( "RP" ),
                               this->m_InputTemplateModel.m_VectorIndexLocations["RP"],
                               cc_RP_Max, "RP" );
        }

      // Local search radius in LR direction is affected by the
      // estimated MSP error in LR direction
      const double err_MSP = std::abs( CandidateRPPoint[0]
                                      - this->m_CenterOfHeadMassEMSP[0] );
      std::cout << "The estimated MSP error in LR direction: "
                << err_MSP << " mm" << std::endl;

      if( err_MSP < 1 )
        {
        searchRadiusLR = 1.;
        std::cout << "Local search radius in LR direction is set to 1 mm." << std::endl;
        }
      else if( err_MSP < 2 )
        {
        searchRadiusLR = 2.;
        std::cout << "Local search radius in LR direction is set to 2 mm." << std::endl;
        }
      else if( err_MSP > 6 )
        {
        itkGenericExceptionMacro(<< "Bad MPJ estimation or too large MSP estimation error!"
                                 << "The estimation result is probably not reliable.");
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
      RPtoCEC = mspSpaceCEC - CandidateRPPoint;

      // RPtoVN4 = this->m_InputTemplateModel.GetRPtoXMean( "VN4" );
      RPtoVN4 = FindVectorFromPointAndVectors
          ( RPtoCEC, this->m_InputTemplateModel.GetRPtoCECMean(),
          this->m_InputTemplateModel.GetRPtoXMean( "VN4" ), -1 );

      double cc_VN4_Max = 0;

      if( this->m_NamedPointEMSP.find( "VN4" ) != this->m_NamedPointEMSP.end() )
        {
        CandidateVN4Point = this->m_NamedPointEMSP["VN4"];
        std::cout << "Skip estimation, directly load from file." << std::endl;
        }
      else if( hasUserForcedVN4Point )
        {
        std::cout << "Skip estimation, directly forced by command line." << std::endl;
        // Points chosen by Slicer should be converted to proper ITK representation
        this->m_OriginalSpaceNamedPoints["VN4"][0] = -(this->m_OriginalSpaceNamedPoints["VN4"][0]);
        this->m_OriginalSpaceNamedPoints["VN4"][1] = -(this->m_OriginalSpaceNamedPoints["VN4"][1]);
        CandidateVN4Point = InvFinalTmsp->TransformPoint( InvHoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["VN4"] ) );
        }
      else
        {
        CandidateVN4Point =
          FindCandidatePoints( this->m_VolumeMSP, mask_LR, searchRadiusLR,
                               1.6 * this->m_TemplateRadius["VN4"],
                               1.6 * this->m_TemplateRadius["VN4"],
                               CandidateRPPoint.GetVectorFromOrigin() + RPtoVN4,
                               this->m_InputTemplateModel.GetTemplateMeans( "VN4" ),
                               this->m_InputTemplateModel.m_VectorIndexLocations["VN4"],
                               cc_VN4_Max, "VN4" );
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
      RPtoAC = FindVectorFromPointAndVectors
          ( RPtoCEC, this->m_InputTemplateModel.GetRPtoCECMean(),
          this->m_InputTemplateModel.GetRPtoXMean( "AC" ), 1 );
      double cc_AC_Max = 0;
      if( this->m_NamedPointEMSP.find( "AC" ) != this->m_NamedPointEMSP.end() )
        {
        CandidateACPoint = this->m_NamedPointEMSP["AC"];
        std::cout << "Skip estimation, directly load from file." << std::endl;
        }
      else if( hasUserForcedACPoint )
        {
        std::cout << "Skip estimation, directly forced by command line." << std::endl;
        // Points chosen by Slicer should be converted to proper ITK representation
        this->m_OriginalSpaceNamedPoints["AC"][0] = -(this->m_OriginalSpaceNamedPoints["AC"][0]);
        this->m_OriginalSpaceNamedPoints["AC"][1] = -(this->m_OriginalSpaceNamedPoints["AC"][1]);
        CandidateACPoint = InvFinalTmsp->TransformPoint( InvHoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["AC"] ) );
        }
      else
        {
        CandidateACPoint =
          FindCandidatePoints( this->m_VolumeMSP, mask_LR, searchRadiusLR,
                               1.6 * this->m_TemplateRadius["AC"],
                               1.6 * this->m_TemplateRadius["AC"],
                               CandidateRPPoint.GetVectorFromOrigin() + RPtoAC,
                               this->m_InputTemplateModel.GetTemplateMeans( "AC" ),
                               this->m_InputTemplateModel.m_VectorIndexLocations["AC"],
                               cc_AC_Max, "AC" );
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
      RPtoPC = FindVectorFromPointAndVectors
          ( RPtoCEC, this->m_InputTemplateModel.GetRPtoCECMean(),
          this->m_InputTemplateModel.GetRPtoXMean( "PC" ), 1 );
      double cc_PC_Max = 0;
      if( this->m_NamedPointEMSP.find( "PC" ) != this->m_NamedPointEMSP.end() )
        {
        CandidatePCPoint = this->m_NamedPointEMSP["PC"];
        std::cout << "Skip estimation, directly load from file." << std::endl;
        }
      else if( hasUserForcedPCPoint )
        {
        std::cout << "Skip estimation, directly forced by command line." << std::endl;
        // Points chosen by Slicer should be converted to proper ITK representation
        this->m_OriginalSpaceNamedPoints["PC"][0] = -(this->m_OriginalSpaceNamedPoints["PC"][0]);
        this->m_OriginalSpaceNamedPoints["PC"][1] = -(this->m_OriginalSpaceNamedPoints["PC"][1]);
        CandidatePCPoint = InvFinalTmsp->TransformPoint( InvHoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["PC"] ) );
        }
      else
        {
        CandidatePCPoint =
          FindCandidatePoints( this->m_VolumeMSP, mask_LR, searchRadiusLR,
                               4 * this->m_TemplateRadius["PC"],
                               4 * this->m_TemplateRadius["PC"],
                               CandidateRPPoint.GetVectorFromOrigin() + RPtoPC,
                               this->m_InputTemplateModel.GetTemplateMeans( "PC" ),
                               this->m_InputTemplateModel.m_VectorIndexLocations["PC"],
                               cc_PC_Max, "PC" );
        }

      // A check point for base landmarks
      if( LMC::globalverboseFlag )
        {
        std::cout << cc_RP_Max << " " << cc_PC_Max << " " << cc_VN4_Max
                  << " " << cc_AC_Max << "\n"
                  << CandidateRPPoint << " " << CandidatePCPoint << " " << CandidateVN4Point << " "
                  << CandidateACPoint << std::endl;
        }

      if( globalImagedebugLevel > 1 )
        {
        std::string BrandedImageAName( this->m_ResultsDir + "/BrandedImage.png" );

        MakeBrandeddebugImage( this->m_VolumeMSP.GetPointer(),
                               this->m_InputTemplateModel,
                               this->m_CenterOfHeadMassEMSP + this->m_InputTemplateModel.GetCMtoRPMean(),
                               CandidateRPPoint + RPtoAC,
                               CandidateRPPoint + RPtoPC,
                               CandidateRPPoint + RPtoVN4,
                               BrandedImageAName,
                               CandidateRPPoint,
                               CandidateACPoint,
                               CandidatePCPoint,
                               CandidateVN4Point
                               );

        if( globalImagedebugLevel > 3 )
          {
          std::string LabelImageAName( this->m_ResultsDir + "/MSP_Mask.nii.gz" );
          MakeLabelImage( this->m_VolumeMSP, CandidateRPPoint, CandidateACPoint,
                          CandidatePCPoint, CandidateVN4Point, LabelImageAName );
          }
        }

      // After finding RP, AC, PC, we can compute the versor transform for
      // registration
      // Also in this stage, we store some results for later use
      // Save named points in original space
      if( !hasUserForcedRPPoint )
        {
        this->m_OriginalSpaceNamedPoints["RP"] =
          this->m_finalTmsp->TransformPoint( CandidateRPPoint );
        }

      if( !hasUserForcedVN4Point )
        {
        this->m_OriginalSpaceNamedPoints["VN4"] =
          this->m_finalTmsp->TransformPoint( CandidateVN4Point );
        }

      if( !hasUserForcedACPoint )
        {
        this->m_OriginalSpaceNamedPoints["AC"] =
          this->m_finalTmsp->TransformPoint( CandidateACPoint );
        }

      if( !hasUserForcedPCPoint )
        {
        this->m_OriginalSpaceNamedPoints["PC"] =
          this->m_finalTmsp->TransformPoint( CandidatePCPoint );
        }

      this->m_OriginalSpaceNamedPoints["CM"] =
        this->m_finalTmsp->TransformPoint( this->m_CenterOfHeadMassEMSP );

      if( !hasUserSpecEyeCenterInfo )
        {
        if( !hasUserForcedRPPoint )
          {
          this->m_OriginalSpaceNamedPoints["RP"] =
            this->m_HoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["RP"] );
          }
        if( !hasUserForcedVN4Point )
          {
          this->m_OriginalSpaceNamedPoints["VN4"] =
          this->m_HoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["VN4"] );
          }
        if( !hasUserForcedACPoint )
          {
          this->m_OriginalSpaceNamedPoints["AC"] =
            this->m_HoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["AC"] );
          }
        if( !hasUserForcedPCPoint )
          {
          this->m_OriginalSpaceNamedPoints["PC"] =
            this->m_HoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["PC"] );
          }
        this->m_OriginalSpaceNamedPoints["CM"] =
          this->m_HoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["CM"] );
        this->m_OriginalSpaceNamedPoints["LE"] = this->m_LEPoint;
        this->m_OriginalSpaceNamedPoints["RE"] = this->m_REPoint;
        }
      else
        {
        this->m_OriginalSpaceNamedPoints["LE"] =
          this->m_finalTmsp->TransformPoint( this->m_NamedPointEMSP["LE"] );
        this->m_OriginalSpaceNamedPoints["RE"] =
          this->m_finalTmsp->TransformPoint( this->m_NamedPointEMSP["RE"] );
        }

      // Write some debug images
        {
        if( globalImagedebugLevel > 3 )
          {
          std::string ResampledMaskmageName
            ( this->m_ResultsDir + "/Resampled_Mask.nii.gz" );
          MakeLabelImage( this->m_VolumeRoughAlignedWithHoughEye, CandidateRPPoint,
                          CandidateACPoint, CandidatePCPoint,
                          CandidateVN4Point, ResampledMaskmageName );

          std::string OrigMaskImageName
            ( this->m_ResultsDir + "/Orig_Mask.nii.gz" );
          MakeLabelImage( this->m_VolumeRoughAlignedWithHoughEye, this->m_OriginalSpaceNamedPoints["RP"],
                          this->m_OriginalSpaceNamedPoints["AC"], this->m_OriginalSpaceNamedPoints["PC"],
                          this->m_OriginalSpaceNamedPoints["VN4"], OrigMaskImageName );
          }
        }

      // Save some named points in EMSP space mainly for debug use
        {
        this->m_NamedPointEMSP["AC"] = CandidateACPoint;
        this->m_NamedPointEMSP["PC"] = CandidatePCPoint;
        this->m_NamedPointEMSP["RP"] = CandidateRPPoint;
        this->m_NamedPointEMSP["VN4"] = CandidateVN4Point;
        this->m_NamedPointEMSP["CM"] = this->m_CenterOfHeadMassEMSP;
        // Eye centers in EMSP have been saved in a earlier time
        }

      // Compute the AC-PC aligned transform
      // Note for the sake of EPCA, we need the transform at this stage
      // so that the method is robust against rotation
      this->m_ImageOrigToACPCVersorTransform = this->ComputeACPCAlignedZeroCenteredTransform();
      //NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)
      VersorTransformType::Pointer LandmarkOrigToACPCTransform = GetLandmarkTransformFromImageTransform( this->m_ImageOrigToACPCVersorTransform.GetPointer()  );

      // Save some named points in AC-PC aligned space
      LandmarksMapType NamedPointACPC;    // named points in ACPC-aligned space
        {
        NamedPointACPC["AC"] =
          LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["AC"] );
        NamedPointACPC["PC"] =
          LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["PC"] );
        NamedPointACPC["RP"] =
          LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["RP"] );
        NamedPointACPC["VN4"] =
          LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["VN4"] );
        NamedPointACPC["CM"] =
          LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["CM"] );
        NamedPointACPC["LE"] =
          LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["LE"] );
        NamedPointACPC["RE"] =
          LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints["RE"] );
        }

      // Get a copy of landmarks on ACPC plane for eliminating accumulative
      // errors of local search process
      LandmarksMapType NamedPointACPCRaw; // reserve this map to avoid the accumulation of local search errors
        {
        NamedPointACPCRaw["AC"] = NamedPointACPC["AC"];
        NamedPointACPCRaw["PC"] = NamedPointACPC["PC"];
        NamedPointACPCRaw["RP"] = NamedPointACPC["RP"];
        NamedPointACPCRaw["VN4"] = NamedPointACPC["VN4"];
        NamedPointACPCRaw["LE"] = NamedPointACPC["LE"];
        NamedPointACPCRaw["RE"] = NamedPointACPC["RE"];
        }

      /*
     * For the rest of the landmarks
     *
     * The search center of the rest of the landmarks will be estimated
     * by linear model estimation with EPCA
     * Note The current version use landmarks in acpc-aligned space.
     */
        {
        // Build up an evolutionary processing list
        // order: RP, AC, PC, VN4, LE, RE, ...
        // Note this order should comply with the order we defined in LLS model
        std::vector<std::string> processingList;
        processingList.push_back( "RP" );
        processingList.push_back( "AC" );
        processingList.push_back( "PC" );
        processingList.push_back( "VN4" );
        processingList.push_back( "LE" );
        processingList.push_back( "RE" );
        unsigned int numBaseLandmarks = 6;
        unsigned int dim = 3;
        for( unsigned int ii = 1; ii <= m_LlsMatrices.size(); ++ii )
          {
          // The processing order is indicated by the length of EPCA coefficient
          std::map<std::string, MatrixType>::iterator
            iit = this->m_LlsMatrices.begin();
          while( iit->second.columns() != ( numBaseLandmarks + ii - 2 ) * dim )
            {
            if( iit != this->m_LlsMatrices.end() )
              {
              ++iit; // transversal
              }
            else
              {
              std::cerr << "Error: wrong number of parameters for linear model!" << std::endl;
              std::cerr << "Some EPCA landmarks may not be detected." << std::endl;
              return;
              }
            }

          std::cout << "Processing " << iit->first << "..." << std::endl;
          if( this->m_NamedPointEMSP.find( iit->first ) != this->m_NamedPointEMSP.end() )
            {
            std::cout << "Skip estimation, directly load from file." << std::endl;
            this->m_OriginalSpaceNamedPoints[iit->first] =
              this->m_finalTmsp->TransformPoint( this->m_NamedPointEMSP[iit->first] );
            NamedPointACPC[iit->first] =
              LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints[iit->first] );
            NamedPointACPCRaw[iit->first] = NamedPointACPC[iit->first];
            }
          else
            {
            // in every iteration, the last name in the processing list
            // indicates the landmark to be estimated
            processingList.push_back( iit->first );

            // Find search center by linear model estimation with
            // dimension reduction.
            // The result will be stored into this->m_NamedPointEMSP[iit->first]
            LinearEstimation( NamedPointACPCRaw, processingList, numBaseLandmarks );

            // check whether it is midline landmark, set search range
            // and modify search center accordingly
            double localSearchRadiusLR = this->m_SearchRadii[iit->first];
              {
              std::vector<std::string>::const_iterator midlineIt =
                this->m_MidlinePointsList.begin();
              for( ; midlineIt != this->m_MidlinePointsList.end(); ++midlineIt )
                {
                if( ( *midlineIt ).compare( iit->first ) == 0 )
                  {
                  localSearchRadiusLR = searchRadiusLR;          // a variable
                                                                 // changed with
                                                                 // err_MSP
                  NamedPointACPCRaw[iit->first][0] = 0.; // search
                                                                 // starts
                                                                 // near EMSP
                  break;
                  }
                }
              }

            NamedPointACPC[iit->first][0] = NamedPointACPCRaw[iit->first][0];
            NamedPointACPC[iit->first][1] = NamedPointACPCRaw[iit->first][1];
            NamedPointACPC[iit->first][2] = NamedPointACPCRaw[iit->first][2];

            // Obtain the position of the current landmark in other spaces
            this->m_OriginalSpaceNamedPoints[iit->first] =
              this->m_ImageOrigToACPCVersorTransform->TransformPoint( NamedPointACPC[iit->first] );
            if( !hasUserSpecEyeCenterInfo )
              {
              this->m_NamedPointEMSP[iit->first] =
                InvHoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints[iit->first] );
              this->m_NamedPointEMSP[iit->first] =
                InvFinalTmsp->TransformPoint( this->m_NamedPointEMSP[iit->first] );
              }
            else
              {
              this->m_NamedPointEMSP[iit->first] =
                InvFinalTmsp->TransformPoint( this->m_OriginalSpaceNamedPoints[iit->first] );
              }

            // Enable local search
            if( 1 )
              {
              // local search
              double cc_Max = 0;
              this->m_NamedPointEMSP[iit->first] =
                FindCandidatePoints( this->m_VolumeMSP, mask_LR, localSearchRadiusLR,
                                     this->m_SearchRadii[iit->first],
                                     this->m_SearchRadii[iit->first],
                                     this->m_NamedPointEMSP[iit->first].GetVectorFromOrigin(),
                                     this->m_InputTemplateModel.GetTemplateMeans( iit->first ),
                                     this->m_InputTemplateModel.m_VectorIndexLocations[iit->first],
                                     cc_Max, iit->first );

              // Update landmarks in input and ACPC-aligned space
              this->m_OriginalSpaceNamedPoints[iit->first] = this->m_finalTmsp->TransformPoint( this->m_NamedPointEMSP[iit->first] );
              if( !hasUserSpecEyeCenterInfo )
                {
                this->m_OriginalSpaceNamedPoints[iit->first] =
                  this->m_HoughEyeTransform->TransformPoint( this->m_OriginalSpaceNamedPoints[iit->first] );
                }
              NamedPointACPC[iit->first] =
                LandmarkOrigToACPCTransform->TransformPoint( this->m_OriginalSpaceNamedPoints[iit->first] );
              }
            }
          } // End of arbitrary landmarks detection for the rest of "new" ones
        }   // End of arbitrary landmarks detection by linear model estimation

      this->ComputeFinalRefinedACPCAlignedTransform();
      if( globalImagedebugLevel > 1 )  // This may be very useful for GUI
                                       // corrector
        {
        WriteITKtoSlicer3Lmk
          ( this->m_ResultsDir + "/EMSP.fcsv",
          this->m_NamedPointEMSP );
        }
      } // End of local searching kernel
    }   // End of local searching
}

void landmarksConstellationDetector::LinearEstimation
  ( LandmarksMapType & namedPoints,
  const std::vector<std::string> & processingList,
  unsigned numBasePoints )
{
  unsigned int dim = namedPoints[processingList[0]].GetPointDimension();

  if( processingList.size() <= numBasePoints )
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
  for( unsigned int k = 1; k <= processingList.size() - 2; ++k )
    {
    for( unsigned int d = 0; d < dim; ++d )
      {
      Xi_t( ( k - 1 ) * dim + d ) = namedPoints[processingList[k]][d]
        - namedPoints[processingList[0]][d] - this->m_LlsMeans[newPointName][d];
      }
    }
  SImageType::PointType::VectorType tmp;
  tmp.SetVnlVector( this->m_LlsMatrices[newPointName] * Xi_t );
  newPoint = namedPoints[processingList[0]] + tmp;
  namedPoints[newPointName] = newPoint;

  // debug
  // std::cout << "Mi' = " << this->m_LlsMatrices[newPointName] << std::endl;
  // std::cout << "Xi_t = " << Xi_t << std::endl;
  // std::cout << "MPJ = " << namedPoints[processingList[0]] << std::endl;
  // std::cout << newPointName << " = " << newPoint << std::endl;
}

SImageType::PointType::VectorType
landmarksConstellationDetector::FindVectorFromPointAndVectors
  ( SImageType::PointType::VectorType BA,
  SImageType::PointType::VectorType BAMean,
  SImageType::PointType::VectorType BCMean,
  int sign )
{
  SImageType::PointType::VectorType BC;

  double cosTheta; // cosine of the angle from BA to BC

  cosTheta = BAMean * BCMean
    / BAMean.GetNorm() / BCMean.GetNorm();

  // Start searching on MSP
  BC[0] = 0;
  double a = BA * BA;
  double b = -2. * BA.GetNorm() * BCMean.GetNorm() * BA[2] * cosTheta;
  double c = BCMean * BCMean
    * ( BA * BA * cosTheta * cosTheta - BA[1] * BA[1] );
  double delta = b * b - 4 * a * c;
  if( delta < 0 )
    {
    itkGenericExceptionMacro(<< "Failed to solve a 2rd-order equation!");
    }
  else if( sign == 1 || sign == -1 )
    {
    BC[2] = -( b - sign * sqrt( delta ) ) / 2. / a;
    }
  else
    {
    itkGenericExceptionMacro(<< "Bad parameter! sign = 1 or sign = -1 please");
    }
  BC[1] = ( BA.GetNorm() * BCMean.GetNorm()
            * cosTheta - BA[2] * BC[2] ) / BA[1];
  return BC;
}
