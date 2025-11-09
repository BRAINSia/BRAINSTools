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

#include "itkLandmarkBasedTransformInitializer.h"
#include <BRAINSFitHelper.h>


std::string
local_to_string(unsigned int i)
{
  std::stringstream localStream;
  localStream << i;
  return localStream.str();
}

// NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)

VersorRigidTransformType::Pointer
landmarksConstellationDetector::GetLandmarkTransformFromImageTransform(
  const VersorRigidTransformType::ConstPointer & orig2msp_img_tfm)
{
  VersorRigidTransformType::Pointer orig2msp_lmk_tfm = VersorRigidTransformType::New();
  const SImageType::PointType       centerPoint = orig2msp_img_tfm->GetCenter();
  orig2msp_lmk_tfm->SetCenter(centerPoint);
  orig2msp_lmk_tfm->SetIdentity();
  orig2msp_img_tfm->GetInverse(orig2msp_lmk_tfm);
  return orig2msp_lmk_tfm;
}

VersorRigidTransformType::Pointer
landmarksConstellationDetector::Compute_orig2msp_img_tfm(const SImagePointType & RP,
                                                         const SImagePointType & AC,
                                                         const SImagePointType & PC)
{
  SImageType::PointType ZeroCenter;
  ZeroCenter.Fill(0.0);

  const Euler3DTransformType::Pointer orig2msp_lmk_tfm_estimate = computeTmspFromPoints(RP, AC, PC, ZeroCenter);

  VersorRigidTransformType::Pointer orig2msp_lmk_tfm_cleaned = VersorRigidTransformType::New();
  orig2msp_lmk_tfm_cleaned->SetFixedParameters(orig2msp_lmk_tfm_estimate->GetFixedParameters());

  itk::Versor<double>               versorRotation;
  const itk::Matrix<double, 3, 3> & CleanedOrthogonalized =
    itk::Orthogonalize3DRotationMatrix(orig2msp_lmk_tfm_estimate->GetMatrix());
  versorRotation.Set(CleanedOrthogonalized);

  orig2msp_lmk_tfm_cleaned->SetRotation(versorRotation);
  orig2msp_lmk_tfm_cleaned->SetTranslation(orig2msp_lmk_tfm_estimate->GetTranslation());
  return orig2msp_lmk_tfm_cleaned;
}

void
landmarksConstellationDetector::ComputeFinalRefinedACPCAlignedTransform(
  const SImageType::Pointer & original_space_image,
  const LandmarksMapType &    updated_orig_lmks)
{
  ////////////////////////////
  // START BRAINSFit alternative
  if (!this->m_atlasVolume.empty())
  {
    using AtlasReaderType = itk::ImageFileReader<SImageType>;
    const AtlasReaderType::Pointer atlasReader = AtlasReaderType::New();
    atlasReader->SetFileName(this->m_atlasVolume);
    std::cout << "read atlas: " << this->m_atlasVolume << std::endl;
    try
    {
      atlasReader->Update();
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cerr << "Error while reading atlasVolume file:\n " << err << std::endl;
    }

    std::cout << "read atlas landmarks:  " << this->m_atlasLandmarks << std::endl;
    LandmarksMapType referenceAtlasLandmarks = ReadSlicer3toITKLmk(this->m_atlasLandmarks);

    // Create a better version of this->m_orig2msp_img_tfm using BRAINSFit.
    // take the the subjects landmarks in original space, and  landmarks from a reference Atlas, and compute an initial
    // affine transform
    // ( using logic from BRAINSLandmarkInitializer) and create initToAtlasAffineTransform.

    LandmarksWeightMapType landmarkWeights;
    if (!this->m_atlasLandmarkWeights.empty())
    {
      landmarkWeights = ReadLandmarkWeights(this->m_atlasLandmarkWeights.c_str());
      std::cout << "read atlas landmarksweights:  " << this->m_atlasLandmarkWeights << std::endl;
    }
    // TEST turning this back on.
    using LmkInitTransformType = itk::AffineTransform<double, Dimension>;
    using LandmarkBasedInitializerType =
      itk::LandmarkBasedTransformInitializer<LmkInitTransformType, SImageType, SImageType>;
    using LandmarkContainerType = LandmarkBasedInitializerType::LandmarkPointContainer;
    LandmarkContainerType                            atlasLmks;
    LandmarkContainerType                            movingLmks;
    LandmarkBasedInitializerType::LandmarkWeightType landmarkWgts;
    for (auto fixedIt = referenceAtlasLandmarks.begin(); fixedIt != referenceAtlasLandmarks.end(); ++fixedIt)
    {
      auto movingIt = updated_orig_lmks.find(fixedIt->first);
      if (movingIt != updated_orig_lmks.cend())
      {
        atlasLmks.push_back(fixedIt->second);
        movingLmks.push_back(movingIt->second);
        if (!this->m_atlasLandmarkWeights.empty())
        {
          if (landmarkWeights.find(fixedIt->first) != landmarkWeights.end())
          {
            landmarkWgts.push_back(landmarkWeights[fixedIt->first]);
          }
          else
          {
            std::cout << "Landmark for " << fixedIt->first << " does not exist. "
                      << "Set the weight to 0.5 " << std::endl;
            landmarkWgts.push_back(0.5F);
          }
        }
      }
      else
      {
        itkGenericExceptionMacro(<< "Could not find " << fixedIt->first << " in originalSpaceLandmarksPreBRAINSFit "
                                 << std::endl
                                 << "MIS MATCHED MOVING AND FIXED LANDMARKS!" << std::endl);
      }
    }

    const LandmarkBasedInitializerType::Pointer landmarkBasedInitializer = LandmarkBasedInitializerType::New();

    if (!this->m_atlasLandmarkWeights.empty())
    {
      landmarkBasedInitializer->SetLandmarkWeight(landmarkWgts);
    }
    landmarkBasedInitializer->SetFixedLandmarks(atlasLmks);
    landmarkBasedInitializer->SetMovingLandmarks(movingLmks);

    const LmkInitTransformType::Pointer initToAtlasAffineTransform = LmkInitTransformType::New();
    landmarkBasedInitializer->SetTransform(initToAtlasAffineTransform);
    landmarkBasedInitializer->InitializeTransform();

    using HelperType = itk::BRAINSFitHelper;
    const HelperType::Pointer brainsFitHelper = HelperType::New();

    // Now Run BRAINSFitHelper class initialized with initToAtlasAffineTransform, original image, and atlas image
    // adapted from BRAINSABC/brainseg/AtlasRegistrationMethod.hxx - do I need to change any of these parameters?
    brainsFitHelper->SetSamplingPercentage(0.05); // Use 5% of voxels for samples
    brainsFitHelper->SetNumberOfHistogramBins(50);
    const std::vector<int> numberOfIterations(1, 1500);
    brainsFitHelper->SetNumberOfIterations(numberOfIterations);
    brainsFitHelper->SetTranslationScale(1000);
    brainsFitHelper->SetReproportionScale(1.0);
    brainsFitHelper->SetSkewScale(1.0);

    using FloatImageType = itk::Image<float, 3>;
    using CastFilterType = itk::CastImageFilter<SImageType, FloatImageType>;

    {
      const CastFilterType::Pointer fixedCastFilter = CastFilterType::New();
      fixedCastFilter->SetInput(atlasReader->GetOutput());
      fixedCastFilter->Update();
      brainsFitHelper->SetFixedVolume(fixedCastFilter->GetOutput());

      const CastFilterType::Pointer movingCastFilter = CastFilterType::New();
      movingCastFilter->SetInput(original_space_image);
      movingCastFilter->Update();
      brainsFitHelper->SetMovingVolume(movingCastFilter->GetOutput());
      itkUtil::WriteImage<FloatImageType>(movingCastFilter->GetOutput(), "./DEBUGFloatCastMovingImage.nii.gz");
    }

    const std::vector<double> minimumStepSize(1, 0.0005);
    brainsFitHelper->SetMinimumStepLength(minimumStepSize);
    std::vector<std::string> transformType(1);
    transformType[0] = "Affine";
    brainsFitHelper->SetTransformType(transformType);

    using CompositeTransformType = itk::CompositeTransform<double, 3>;
    CompositeTransformType::Pointer initToAtlasAffineCompositeTransform =
      dynamic_cast<CompositeTransformType *>(initToAtlasAffineTransform.GetPointer());
    if (initToAtlasAffineCompositeTransform.IsNull())
    {
      initToAtlasAffineCompositeTransform = CompositeTransformType::New();
      initToAtlasAffineCompositeTransform->AddTransform(initToAtlasAffineTransform);
    }
    brainsFitHelper->SetCurrentGenericTransform(initToAtlasAffineCompositeTransform);

    // Provide better masking of images
    {
      static constexpr unsigned int ROIAutoClosingSize = 4;
      static constexpr unsigned int ROIAutoDilateSize = 6;
      {
        using ROIAutoType = itk::BRAINSROIAutoImageFilter<SImageType, itk::Image<unsigned char, 3>>;
        const ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
        ROIFilter->SetInput(atlasReader->GetOutput());
        ROIFilter->SetClosingSize(ROIAutoClosingSize);
        ROIFilter->SetDilateSize(ROIAutoDilateSize);
        ROIFilter->Update();
        const ImageMaskPointer fixedMask = ROIFilter->GetSpatialObjectROI();
        brainsFitHelper->SetFixedBinaryVolume(fixedMask);
      }
      {
        using ROIAutoType = itk::BRAINSROIAutoImageFilter<SImageType, itk::Image<unsigned char, 3>>;
        const ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
        ROIFilter->SetInput(original_space_image);
        ROIFilter->SetClosingSize(ROIAutoClosingSize);
        ROIFilter->SetDilateSize(ROIAutoDilateSize);
        ROIFilter->Update();
        const ImageMaskPointer movingMask = ROIFilter->GetSpatialObjectROI();
        brainsFitHelper->SetMovingBinaryVolume(movingMask);
      }
    }
    brainsFitHelper->SetDebugLevel(10);
    brainsFitHelper->Update();

    this->m_orig2msp_img_tfm = itk::ComputeRigidTransformFromGeneric(
      brainsFitHelper->GetCurrentGenericTransform()->GetNthTransform(0).GetPointer());
    if (this->m_orig2msp_img_tfm.IsNull())
    {
      // Fail if something weird happens.
      itkGenericExceptionMacro(<< "this->m_orig2msp_img_tfm is null. "
                               << "It means we're not registering to the atlas, after all." << std::endl);
    }

    {
      // NOTE: LandmarkTransforms are inverse of ImageTransforms, (You pull images, you push landmarks)
      using VersorRigid3DTransformType = itk::VersorRigid3DTransform<double>;
      const VersorRigidTransformType::Pointer orig2msp_lmk_tfm =
        GetLandmarkTransformFromImageTransform(this->m_orig2msp_img_tfm.GetPointer());
      const VersorRigid3DTransformType::OutputPointType acPointInACPCSpace =
        orig2msp_lmk_tfm->TransformPoint(GetNamedPointFromLandmarkList(updated_orig_lmks, "AC"));
      {
        VersorRigid3DTransformType::OffsetType translation;
        translation[0] =
          +acPointInACPCSpace[0]; // NOTE: Positive translation in ImageTransform is Negative in LandmarkTransform
        translation[1] = +acPointInACPCSpace[1];
        translation[2] = +acPointInACPCSpace[2];
        // First shift the transform
        this->m_orig2msp_img_tfm->Translate(translation, false);
      }
      // INFO:  This still does not put it to (0,0,0) and it should.
    }
  }
  /// END BRAINSFIT_ALTERNATIVE
  ////////////////////////////
}

VersorRigidTransformType::Pointer
landmarksConstellationDetector::GetImageOrigToACPCVersorTransform() const
{
  return m_orig2msp_img_tfm;
}

SImageType::PointType
landmarksConstellationDetector::FindCandidatePoints(
  const SImageType::Pointer & volumeMSP,
  const SImageType::Pointer & mask_LR,
  const double                LR_restrictions,
  const double                PA_restrictions,
  const double                SI_restrictions,
  // INFO: restrictions should really be ellipsoidal values
  const SImageType::PointType::VectorType &                      CenterOfSearchArea,
  const std::vector<std::vector<float>> &                        TemplateMean,
  const landmarksConstellationModelIO::IndexLocationVectorType & model,
  double &                                                       cc_Max,
  const std::string &                                            mapID)
{
  cc_Max = -123456789.0;

  const LinearInterpolatorType::Pointer imInterp = LinearInterpolatorType::New();
  imInterp->SetInputImage(volumeMSP);
  const LinearInterpolatorType::Pointer maskInterp = LinearInterpolatorType::New();
  maskInterp->SetInputImage(mask_LR);

  // Final location is initialized with the center of search area
  const SImageType::PointType InitialGuessPoint{ CenterOfSearchArea };

  // Boundary check
  {
    SImageType::PointType boundaryL(InitialGuessPoint);
    boundaryL[0] += LR_restrictions;
    SImageType::PointType boundaryR(InitialGuessPoint);
    boundaryR[0] -= LR_restrictions;
    SImageType::PointType boundaryP(InitialGuessPoint);
    boundaryP[1] += PA_restrictions;
    SImageType::PointType boundaryA(InitialGuessPoint);
    boundaryA[1] -= PA_restrictions;
    SImageType::PointType boundaryS(InitialGuessPoint);
    boundaryS[2] += SI_restrictions;
    SImageType::PointType boundaryI(InitialGuessPoint);
    boundaryI[2] -= SI_restrictions;
    if ((!maskInterp->IsInsideBuffer(boundaryL)) || (!maskInterp->IsInsideBuffer(boundaryR)) ||
        (!maskInterp->IsInsideBuffer(boundaryP)) || (!maskInterp->IsInsideBuffer(boundaryA)) ||
        (!maskInterp->IsInsideBuffer(boundaryS)) || (!maskInterp->IsInsideBuffer(boundaryI)))
    {
      std::cout << "WARNING: search region outside of the image region for " << mapID << "." << std::endl;
      std::cout << "The detection has probably large error!"
                << "\nUsing point " << InitialGuessPoint << " for " << mapID << std::endl;
      return InitialGuessPoint;
    }
  }

  // height and radius of the moving template
  const double height = this->m_InputTemplateModel.GetHeight(mapID);
  // HACK:
  // When, rmpj, rac, rpc and rvn4 are used, they help to increase the bounding box
  // restrictions; however, we do not want that landmark template image be affected
  // by that.
  const double radii = ((mapID == "RP" || mapID == "AC" || mapID == "PC" || mapID == "VN4") &&
                        this->m_InputTemplateModel.GetRadius(mapID) > 5.0)
                         ? 5.0
                         : this->m_InputTemplateModel.GetRadius(mapID);

  constexpr double deltaLR = 1; // in mm
  constexpr double deltaPA = 1; // in mm
  constexpr double deltaSI = 1; // in mm

  // Bounding box around the center point.
  // To compute the correlation for the border points, the bounding box needs
  // to be expanded by the template size.
  //
  SImageType::PointType LPS_BEGIN;
  LPS_BEGIN[0] = CenterOfSearchArea[0] - LR_restrictions - height / 2;
  LPS_BEGIN[1] = CenterOfSearchArea[1] - PA_restrictions - radii;
  LPS_BEGIN[2] = CenterOfSearchArea[2] - SI_restrictions - radii;

  SImageType::PointType LPS_END;
  LPS_END[0] = CenterOfSearchArea[0] + LR_restrictions + height / 2;
  LPS_END[1] = CenterOfSearchArea[1] + PA_restrictions + radii;
  LPS_END[2] = CenterOfSearchArea[2] + SI_restrictions + radii;

  // Now bounding area will be converted to an image
  //
  const SImageType::Pointer roiImage = SImageType::New();
  {
    // origin
    SImageType::PointType roiOrigin;
    roiOrigin[0] = LPS_BEGIN[0];
    roiOrigin[1] = LPS_BEGIN[1];
    roiOrigin[2] = LPS_BEGIN[2];
    roiImage->SetOrigin(roiOrigin);
  }
  {
    // size
    SImageType::SizeType roiSize;
    roiSize[0] = static_cast<SImageType::SizeType::SizeValueType>(1.0 + LPS_END[0] - LPS_BEGIN[0]);
    roiSize[1] = static_cast<SImageType::SizeType::SizeValueType>(1.0 + LPS_END[1] - LPS_BEGIN[1]);
    roiSize[2] = static_cast<SImageType::SizeType::SizeValueType>(1.0 + LPS_END[2] - LPS_BEGIN[2]);
    // start index
    SImageType::IndexType roiStart;
    roiStart.Fill(0);
    // region
    const SImageType::RegionType roiRegion(roiStart, roiSize);
    roiImage->SetRegions(roiRegion);
  }
  { // Default to spacing of 1 and size identity direction cosine
    SImageType::SpacingType spacing;
    spacing[0] = deltaLR;
    spacing[1] = deltaPA;
    spacing[2] = deltaSI;
    // roiImage->SetSpacing(volumeMSP->GetSpacing());  // spacing
    // roiImage->SetDirection(volumeMSP->GetDirection());  // direction
  }
  roiImage->Allocate(true); // true implies roiImage->FillBuffer(0);

  // Since the actual bounding box is a rounded area, a Roi mask is also needed.
  //
  const SImageType::Pointer roiMask = SImageType::New();
  roiMask->CopyInformation(roiImage);
  roiMask->SetRegions(roiImage->GetLargestPossibleRegion());
  roiMask->Allocate(true); // true implies roiMask->FillBuffer(0);

  // roiImage is filled with values from volumeMSP
  //
  {
    SImageType::PointType currentPointLocation;
    for (currentPointLocation[0] = LPS_BEGIN[0]; currentPointLocation[0] < LPS_END[0];
         currentPointLocation[0] += deltaLR)
    {
      for (currentPointLocation[1] = LPS_BEGIN[1]; currentPointLocation[1] < LPS_END[1];
           currentPointLocation[1] += deltaPA)
      {
        for (currentPointLocation[2] = LPS_BEGIN[2]; currentPointLocation[2] < LPS_END[2];
             currentPointLocation[2] += deltaSI)
        {
          // Is current point within the input mask
          if (maskInterp->Evaluate(currentPointLocation) > 0.5)
          {
            // Is current point inside the boundary box
            const SImageType::PointType::VectorType temp =
              currentPointLocation.GetVectorFromOrigin() - CenterOfSearchArea;
            const double inclusionDistance = temp.GetNorm();
            if ((inclusionDistance < (SI_restrictions + radii)) && (std::abs(temp[1]) < (PA_restrictions + radii)))
            {
              SImageType::IndexType index3D;
              index3D = roiImage->TransformPhysicalPointToIndex(currentPointLocation);
              roiImage->SetPixel(index3D, imInterp->Evaluate(currentPointLocation));
              roiMask->SetPixel(index3D, 1);
            }
          }
        }
      }
    }
  }
  ////////
  // Now we need to normalize only bounding region inside the roiImage
  ///////
  using BinaryImageToLabelMapFilterType = itk::BinaryImageToLabelMapFilter<SImageType>;
  const BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
  binaryImageToLabelMapFilter->SetInput(roiMask);
  binaryImageToLabelMapFilter->Update();

  using LabelMapToLabelImageFilterType =
    itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, SImageType>;
  const LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
  labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
  labelMapToLabelImageFilter->Update();

  using LabelStatisticsImageFilterType = itk::LabelStatisticsImageFilter<SImageType, SImageType>;
  const LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetLabelInput(labelMapToLabelImageFilter->GetOutput());
  labelStatisticsImageFilter->SetInput(roiImage);
  labelStatisticsImageFilter->Update();

  if (labelStatisticsImageFilter->GetNumberOfLabels() != 1)
  {
    itkGenericExceptionMacro(<< "The bounding box mask should be connected for " << mapID << ".");
  }

  using LabelPixelType = LabelStatisticsImageFilterType::LabelPixelType;
  const LabelPixelType labelValue = labelStatisticsImageFilter->GetValidLabelValues()[0];
  const double         ROImean = labelStatisticsImageFilter->GetMean(labelValue);
  const double         ROIvar = labelStatisticsImageFilter->GetVariance(labelValue);
  const unsigned long  ROIcount = labelStatisticsImageFilter->GetCount(labelValue);

  // The area inside the bounding box is normalized using the mean and variance statistics
  using SubtractImageFilterType = itk::SubtractImageFilter<SImageType, SImageType, FImageType3D>;
  const SubtractImageFilterType::Pointer subtractConstantFromImageFilter = SubtractImageFilterType::New();
  subtractConstantFromImageFilter->SetInput(roiImage);
  subtractConstantFromImageFilter->SetConstant2(ROImean);
  subtractConstantFromImageFilter->Update();

  if (std::sqrt(ROIcount * ROIvar) < std::numeric_limits<double>::epsilon())
  {
    if (globalImagedebugLevel > 8)
    {
      const std::string bbmask(this->m_ResultsDir + "/bbmaks_" + itksys::SystemTools::GetFilenameName(mapID) +
                               ".nii.gz");
      itkUtil::WriteImage<SImageType>(roiMask, bbmask);
    }
    std::cerr << "WARNING: search region has no variance (uniform values), or i size 0: " << mapID << "." << std::endl;
    std::cerr << "The detection has probably large error!" << std::endl;
    std::cerr << "Zero norm for bounding area for " << mapID << ".\n"
              << "ROI Pixel Count: " << ROIcount << " ROIvariance: " << ROIvar << std::endl;

    return InitialGuessPoint;
  }
  const double normInv = 1 / (std::sqrt(ROIcount * ROIvar));

  using MultiplyImageFilterType = itk::MultiplyImageFilter<FImageType3D, FImageType3D, FImageType3D>;
  const MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(subtractConstantFromImageFilter->GetOutput());
  multiplyImageFilter->SetConstant(normInv);

  const FImageType3D::Pointer normalizedRoiImage = multiplyImageFilter->GetOutput();
  /////////////// End of normalization of roiImage //////////////
  // Now each landmark template should be converted to a moving template image
  //
  const FImageType3D::Pointer lmkTemplateImage = FImageType3D::New();
  lmkTemplateImage->SetOrigin(roiImage->GetOrigin());
  lmkTemplateImage->SetSpacing(roiImage->GetSpacing());
  lmkTemplateImage->SetDirection(roiImage->GetDirection());
  FImageType3D::SizeType mi_size;
  mi_size[0] = 2 * height + 1;
  mi_size[1] = 2 * radii + 1;
  mi_size[2] = 2 * radii + 1;
  FImageType3D::IndexType mi_start;
  mi_start.Fill(0);
  const FImageType3D::RegionType mi_region(mi_start, mi_size);
  lmkTemplateImage->SetRegions(mi_region);
  lmkTemplateImage->Allocate();

  // Since each landmark template is a cylinder, a template mask is needed.
  //
  const SImageType::Pointer templateMask = SImageType::New();
  templateMask->CopyInformation(lmkTemplateImage);
  templateMask->SetRegions(lmkTemplateImage->GetLargestPossibleRegion());
  templateMask->Allocate();

  // Fill the template moving image based on the vector index locations
  //  and template mean values for different angle rotations
  //
  double                cc_rotation_max = 0.0;
  SImageType::PointType TransformedGuessPoint = InitialGuessPoint;
  for (unsigned int curr_rotationAngle_index = 0; curr_rotationAngle_index < TemplateMean.size();
       curr_rotationAngle_index++)
  {
    lmkTemplateImage->FillBuffer(0);
    templateMask->FillBuffer(0);
    // iterate over mean values for the current rotation angle
    auto mean_iter = TemplateMean[curr_rotationAngle_index].begin();
    // Fill the lmk template image using the mean values
    for (auto it = model.begin(); it != model.end(); ++it, ++mean_iter)
    {
      FImageType3D::IndexType pixelIndex;
      pixelIndex[0] = (*it)[0] + height;
      pixelIndex[1] = (*it)[1] + radii;
      pixelIndex[2] = (*it)[2] + radii;
      lmkTemplateImage->SetPixel(pixelIndex, *mean_iter);
      templateMask->SetPixel(pixelIndex, 1);
    }

    // Finally NCC is calculated in frequency domain
    //
    using CorrelationFilterType =
      itk::MaskedFFTNormalizedCorrelationImageFilter<FImageType3D, FImageType3D, SImageType>;
    const CorrelationFilterType::Pointer correlationFilter = CorrelationFilterType::New();
    correlationFilter->SetFixedImage(normalizedRoiImage);
    correlationFilter->SetFixedImageMask(roiMask);
    correlationFilter->SetMovingImage(lmkTemplateImage);
    correlationFilter->SetMovingImageMask(templateMask);
    correlationFilter->SetRequiredFractionOfOverlappingPixels(1);
    correlationFilter->Update();
    if (globalImagedebugLevel > 8)
    {

      LandmarksMapType msp_lmks_algo_found; // named points in EMSP space
      {
        SImageType::PointType centerOfSearchAreaPoint;
        for (unsigned int i = 0; i < centerOfSearchAreaPoint.size(); i++)
        {
          centerOfSearchAreaPoint[i] = CenterOfSearchArea[i];
        }
        msp_lmks_algo_found["CenterOfSearchArea"] = centerOfSearchAreaPoint;
      }
      msp_lmks_algo_found["LPS_BEGIN"] = LPS_BEGIN;
      msp_lmks_algo_found["LPS_END"] = LPS_END;
      const std::string pointsName(this->m_ResultsDir + "/NCCOutput_" + itksys::SystemTools::GetFilenameName(mapID) +
                                   "_" + local_to_string(curr_rotationAngle_index) + "_markers.fcsv");
      WriteITKtoSlicer3Lmk(pointsName, msp_lmks_algo_found);

      // Print the name of the file that is the master volume, the fixed image sample should
      // be aligned with this!
      const std::string volumeMSPName(this->m_ResultsDir + "/NCCOutput_" + itksys::SystemTools::GetFilenameName(mapID) +
                                      "_" + local_to_string(curr_rotationAngle_index) + "_volumeMSP.nii.gz");
      itkUtil::WriteImage<SImageType>(volumeMSP, volumeMSPName);


      const std::string mask_LRName(this->m_ResultsDir + "/NCCOutput_" + itksys::SystemTools::GetFilenameName(mapID) +
                                    "_" + local_to_string(curr_rotationAngle_index) + "_mask_LRName.nii.gz");
      itkUtil::WriteImage<SImageType>(mask_LR, mask_LRName);

      const std::string ncc_output_name_fixed(this->m_ResultsDir + "/NCCOutput_" +
                                              itksys::SystemTools::GetFilenameName(mapID) + "_" +
                                              local_to_string(curr_rotationAngle_index) + "_fixed.nii.gz");
      itkUtil::WriteImage<FImageType3D>(normalizedRoiImage, ncc_output_name_fixed);

      const std::string ncc_output_name_fixedmask(this->m_ResultsDir + "/NCCOutput_" +
                                                  itksys::SystemTools::GetFilenameName(mapID) + "_" +
                                                  local_to_string(curr_rotationAngle_index) + "_fixedmask.nii.gz");
      itkUtil::WriteImage<SImageType>(roiMask, ncc_output_name_fixedmask);

      const std::string ncc_output_name_moving(this->m_ResultsDir + "/NCCOutput_" +
                                               itksys::SystemTools::GetFilenameName(mapID) + "_" +
                                               local_to_string(curr_rotationAngle_index) + "_moving.nii.gz");
      itkUtil::WriteImage<FImageType3D>(lmkTemplateImage, ncc_output_name_moving);

      const std::string ncc_output_name_movingmask(this->m_ResultsDir + "/NCCOutput_" +
                                                   itksys::SystemTools::GetFilenameName(mapID) + "_" +
                                                   local_to_string(curr_rotationAngle_index) + "_movingmask.nii.gz");
      itkUtil::WriteImage<SImageType>(templateMask, ncc_output_name_movingmask);


      const std::string ncc_output_name(this->m_ResultsDir + "/NCCOutput_" +
                                        itksys::SystemTools::GetFilenameName(mapID) + "_" +
                                        local_to_string(curr_rotationAngle_index) + ".nii.gz");
      itkUtil::WriteImage<FImageType3D>(correlationFilter->GetOutput(), ncc_output_name);
    }

    // Maximum NCC for current rotation angle
    using MinimumMaximumImageCalculatorType = itk::MinimumMaximumImageCalculator<FImageType3D>;
    const MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter =
      MinimumMaximumImageCalculatorType::New();
    minimumMaximumImageCalculatorFilter->SetImage(correlationFilter->GetOutput());
    minimumMaximumImageCalculatorFilter->Compute();
    const double cc = minimumMaximumImageCalculatorFilter->GetMaximum();
    if (cc > cc_rotation_max)
    {
      cc_rotation_max = cc;
      // Where maximum happens
      const FImageType3D::IndexType maximumCorrelationPatchCenter =
        minimumMaximumImageCalculatorFilter->GetIndexOfMaximum();
      correlationFilter->GetOutput()->TransformIndexToPhysicalPoint(maximumCorrelationPatchCenter,
                                                                    TransformedGuessPoint);
    }
  }
  cc_Max = cc_rotation_max;

  if (LMC::globalverboseFlag)
  {
    std::cout << "cc max: " << cc_Max << std::endl;
    std::cout << mapID << " Final in MSP aligned physical space: " << TransformedGuessPoint << std::endl;
  }
  return TransformedGuessPoint;
}

void
landmarksConstellationDetector::EulerToVersorRigid(VersorRigidTransformType::Pointer &        result,
                                                   const Euler3DTransformType::ConstPointer & eulerRigid)
{
  if (result.IsNotNull() && eulerRigid.IsNotNull())
  {
    result->SetFixedParameters(eulerRigid->GetFixedParameters());
    itk::Versor<double>               versorRotation;
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized =
      itk::Orthogonalize3DRotationMatrix(eulerRigid->GetMatrix());
    versorRotation.Set(CleanedOrthogonalized);
    result->SetRotation(versorRotation);
    result->SetTranslation(eulerRigid->GetTranslation());
  }
  else
  {
    itkGenericExceptionMacro(<< "Error missing Pointer data, assigning "
                             << "Euler3DTransformPointer to VersorRigid3DTransformPointer." << std::endl);
  }
}


void
landmarksConstellationDetector::DoResampleInPlace(const SImageType::ConstPointer &           inputImg,
                                                  const Euler3DTransformType::ConstPointer & rigidTx,
                                                  SImageType::Pointer &                      inPlaceResampledImg)
{
  VersorRigidTransformType::Pointer versorRigidTx = VersorRigidTransformType::New();
  EulerToVersorRigid(versorRigidTx, rigidTx.GetPointer());

  using ResampleIPFilterType = itk::ResampleInPlaceImageFilter<SImageType, SImageType>;
  using ResampleIPFilterPointer = ResampleIPFilterType::Pointer;
  const ResampleIPFilterPointer InPlaceResampler = ResampleIPFilterType::New();
  InPlaceResampler->SetInputImage(inputImg);
  InPlaceResampler->SetRigidTransform(versorRigidTx.GetPointer());
  InPlaceResampler->Update();

  inPlaceResampledImg = InPlaceResampler->GetOutput();
}

void
landmarksConstellationDetector::LinearEstimation(LandmarksMapType &               msp_lmks_linearly_estimated,
                                                 const std::vector<std::string> & processingList,
                                                 unsigned                         numBasePoints)
{
  const unsigned int dim = msp_lmks_linearly_estimated[processingList[0]].GetPointDimension();

  if (processingList.size() <= numBasePoints)
  {
    std::cout << "No EPCA landmarks to be estimated." << std::endl;
    return;
  }
  const std::string &   newPointName = processingList[processingList.size() - 1];
  SImageType::PointType newPoint;
  newPoint.Fill(0);

  // Construct Xi_t
  VectorType Xi_t;
  Xi_t.set_size(dim * (processingList.size() - 2));
  for (unsigned int k = 1; k <= processingList.size() - 2; ++k)
  {
    for (unsigned int d = 0; d < dim; ++d)
    {
      Xi_t((k - 1) * dim + d) = msp_lmks_linearly_estimated[processingList[k]][d] -
                                msp_lmks_linearly_estimated[processingList[0]][d] - this->m_LlsMeans[newPointName][d];
    }
  }
  SImageType::PointType::VectorType tmp;
  tmp.SetVnlVector(this->m_LlsMatrices[newPointName] * Xi_t);
  newPoint = msp_lmks_linearly_estimated[processingList[0]] + tmp;
  msp_lmks_linearly_estimated[newPointName] = newPoint;
}

SImageType::PointType::VectorType
landmarksConstellationDetector::FindVectorFromPointAndVectors(SImageType::PointType::VectorType BA,
                                                              SImageType::PointType::VectorType BAMean,
                                                              SImageType::PointType::VectorType BCMean,
                                                              int                               sign)
{
  SImageType::PointType::VectorType BC;

  double cosTheta = NAN; // cosine of the angle from BA to BC

  cosTheta = BAMean * BCMean / BAMean.GetNorm() / BCMean.GetNorm();

  // Start searching on MSP
  BC[0] = 0;
  const double a = BA * BA;
  const double b = -2. * BA.GetNorm() * BCMean.GetNorm() * BA[2] * cosTheta;
  const double c = BCMean * BCMean * (BA * BA * cosTheta * cosTheta - BA[1] * BA[1]);
  const double delta = b * b - 4 * a * c;
  if (delta < 0)
  {
    itkGenericExceptionMacro(<< "Failed to solve a 2rd-order equation!");
  }
  else if (sign == 1 || sign == -1)
  {
    BC[2] = -(b - sign * sqrt(delta)) / 2. / a;
  }
  else
  {
    itkGenericExceptionMacro(<< "Bad parameter! sign = 1 or sign = -1 please");
  }
  BC[1] = (BA.GetNorm() * BCMean.GetNorm() * cosTheta - BA[2] * BC[2]) / BA[1];
  return BC;
}


void
WriteManualFixFiles(const std::string &      EMSP_Fiducial_file_name,
                    SImageType * const       mspVolume,
                    const std::string &      resultDir,
                    const LandmarksMapType & errorLmks,
                    const std::string &      failureMessage,
                    const bool               throwException)
{ // ADD MetaData for EMSP_FCSV_FILENAME
  itk::MetaDataDictionary & dict = mspVolume->GetMetaDataDictionary();
  const char * const        metaDataEMSP_FCSVName = "EMSP_FCSV_FILENAME";
  itk::EncapsulateMetaData<std::string>(dict, metaDataEMSP_FCSVName, EMSP_Fiducial_file_name.c_str());

  // write EMSP aligned image
  itkUtil::WriteImage<SImageType>(mspVolume, resultDir + "/EMSP.nrrd");

  if (!errorLmks.empty())
  {
    WriteITKtoSlicer3Lmk(resultDir + "/" + EMSP_Fiducial_file_name, errorLmks);
  }
  if (throwException)
  {
    itkGenericExceptionMacro(<< failureMessage);
  }
}
