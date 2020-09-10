//
// Created by johnsonhj on 7/25/20.
// A program to process an entire set of T1, T2, PD images at once
// and create a deface mask suitable for anonymizing data

//
#include <iostream>
#include <list>
#include <itkSmoothingRecursiveGaussianImageFilter.h>

#include <BRAINSDefaceCLP.h>
#include <itkIndexRange.h>

#include <itkIO.h>
#include <Slicer3LandmarkIO.h>
#include <BRAINSIntensityTransform.h>
#include <itkBRAINSROIAutoImageFilter.h>
#include "itkMaskImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "BRAINSCommonLib.h"
#include "itkUnaryGeneratorImageFilter.h"


/*
 * This program takes in a landmark file with at least LE, and RE points defined, and series of ACPC aligned data
 * (i.e. AC at ~(0,0,0) and PC at ~(0,X,0) to generate "defaced" or "face-obscured" data as a de-identification
 * process.
 *
 * A coded mask image is created in a standardized 0.5mm space that is 512x512x512 with the center of the image at
 * 0,0,0.
 *
 */


int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  // usage check
  if (!inputMask.empty() && noMaskApplication)
  {
    std::cout << "A mask has been provided, but the noMaskApplication flag has been set.\nThis program will do nothing "
                 "in this configuration.";
    exit(-1);
  }

  // STEP 1: Using all the anatomical images for this session, generate a label mask defining
  // regions that are 0= in all the images
  //                  1= part of the face region
  //                  2= 80mm below the ac point --> altered to remove this region by resampling
  //                  3= Outside FOV for at least 1 image
  using InternalImageType = itk::Image<float, 3>;
  using MaskImageType = itk::Image<unsigned char, 3>;

  const auto lmks = ReadSlicer3toITKLmk(inputLandmarks);

  MaskImageType::PointType AC_pnt = lmks.at("AC");
  MaskImageType::PointType RE_pnt = lmks.at("RE");
  MaskImageType::PointType LE_pnt = lmks.at("LE");

  constexpr size_t LR = 0;
  constexpr size_t PA = 1;
  constexpr size_t SI = 2;

  // First add images used for mask generation
  std::vector<InternalImageType::Pointer> img_list;
  std::vector<std::string>                img_filenames;
  for (const auto & im_fn : inputVolume)
  {
    img_list.emplace_back(itkUtil::ReadImage<InternalImageType>(im_fn));
    img_filenames.emplace_back(im_fn);
  }
  inputVolume.clear();

  constexpr int valid_inside_pixel = 0;
  constexpr int face_rm = 1;
  constexpr int below_ac = 2;
  constexpr int outside_fov = 3;
  constexpr int auto_roi_background = 4;
  constexpr int eye_boxes_code = 5;
  // Find all the out of FOV spaces fro the mask image
  MaskImageType::PointType     maskpnt;
  InternalImageType::IndexType imgindex;

  using FadeMapType = itk::Image<float, 3>;

  MaskImageType::Pointer mask_labels;


  if (!inputMask.empty())
  { // Read pre-generated mask
    mask_labels = itkUtil::ReadImage<MaskImageType>(inputMask);
  }
  else
  { // Generate a new mask

    // Make a 1/2 mm masked volume
    mask_labels = MaskImageType::New();
    {
      const double                               spacing[3] = { 0.5, 0.5, 0.5 };
      const MaskImageType::RegionType::IndexType starting_index{ 0, 0, 0 };
      const MaskImageType::RegionType::SizeType  img_size{ 512, 512, 352 };
      // 512 - 160 (80mm @ 0.5 spacing) = 352
      // Empirically chosen origin voxel based on standard center of brain  WRONG:   ->  -128.0, -114.0, -132.0
      const double origin[3] = { AC_pnt[LR] - 128.0, AC_pnt[PA] - 128.0, AC_pnt[SI] - 80.0 };
      // -132 + 80  = 52 removing 80mm inferior
      // -114 - 40 = -154 adding 40mm anterior

      MaskImageType::RegionType region;
      region.SetSize(img_size);
      region.SetIndex(starting_index);
      mask_labels->SetRegions(region);

      MaskImageType::DirectionType id;
      id.SetIdentity();
      mask_labels->SetDirection(id);

      mask_labels->SetSpacing(spacing);
      mask_labels->SetOrigin(origin);
      mask_labels->Allocate(true);
    }

    // ROIAuto declaration
    // unused typedef below
    //  using ROIAutoType = itk::BRAINSROIAutoImageFilter<InternalImageType, MaskImageType>;
    // Unused variable below
    //  constexpr double             max_smoothing_size = 0.0; // Try 10.0 for very safe margins

    for (auto & curr_img : img_list)
    {
      //    MaskImageType::Pointer roi = [](InternalImageType::Pointer current_image) -> MaskImageType::Pointer {
      //      ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
      //      ROIFilter->SetClosingSize(20);                                  // close 10mm holes
      //      ROIFilter->SetDilateSize(static_cast<int>(max_smoothing_size)); // 55 mm ~ IPD ~ tip of the nose
      //
      //      // pre-compute ROIAuto
      //      ROIFilter->SetInput(current_image);
      //      ROIFilter->Update();
      //
      //      return ROIFilter->GetOutput();
      //    }(curr_img);
      //    itk::NearestNeighborInterpolateImageFunction<MaskImageType>::Pointer roiAutoInterpolator =
      //      itk::NearestNeighborInterpolateImageFunction<MaskImageType>::New();
      //    roiAutoInterpolator->SetInputImage(roi);
      //    if (debugLevel >= 2)
      //    {
      //      itkUtil::WriteImage<MaskImageType>(roi, outputDirectory + "/roi.nii.gz");
      //    }
      {
        itk::ImageRegionIteratorWithIndex<MaskImageType> mit(mask_labels, mask_labels->GetLargestPossibleRegion());
        constexpr double                                 approx_eye_radius = 9; // Size of eyes
        while (!mit.IsAtEnd())
        {
          mask_labels->TransformIndexToPhysicalPoint<double>(mit.GetIndex(), maskpnt);
          const bool isInside = curr_img->TransformPhysicalPointToIndex(maskpnt, imgindex);
          if (!isInside)
          {
            mit.Set(outside_fov);
          }
          if (mit.Value() == 0) // Only change if not yet set.
          {
            if (maskpnt[SI] < AC_pnt[SI] - 80.0) // removing 80 mm inferior
            {
              mit.Set(below_ac);
            }
            // blur nose region
            else if ((maskpnt[PA] < RE_pnt[PA] || maskpnt[PA] < LE_pnt[PA]) &&
                     (maskpnt[SI] < (RE_pnt[SI]) || (maskpnt[SI] < LE_pnt[SI])))
            {
              mit.Set(face_rm);
            }
            // blur Cheak region
            else if ((maskpnt[PA] < RE_pnt[PA] || maskpnt[PA] < LE_pnt[PA]) &&
                     (maskpnt[SI] < (RE_pnt[SI]) || (maskpnt[SI] < LE_pnt[SI])))
            {
              mit.Set(face_rm);
            }

            else if ((
                       // Now remove RE eye boxes
                       (maskpnt[PA] < (RE_pnt[PA] + approx_eye_radius)) && // anterior to back of eye
                       (maskpnt[SI] < (RE_pnt[SI] + approx_eye_radius)) && // in si eye region
                       (maskpnt[LR] < (RE_pnt[LR] + approx_eye_radius))    // lateral to eye
                       ) ||
                     (
                       // Now remove LE eye boxes
                       (maskpnt[PA] < (RE_pnt[PA] + approx_eye_radius)) && // anterior to back of eye
                       (maskpnt[SI] < (RE_pnt[SI] + approx_eye_radius)) && // in si eye region
                       (maskpnt[LR] > (LE_pnt[LR] - approx_eye_radius)))   // lateral to eye
            )
            {
              mit.Set(eye_boxes_code);
            }
            // Now try to remove eyebrows/forehead
            else if ((maskpnt[PA] < (RE_pnt[PA] - approx_eye_radius * 3) ||
                      maskpnt[PA] < (LE_pnt[PA] - approx_eye_radius * 3))
                     //                   && (maskpnt[SI] < (RE_pnt[SI] + approx_eye_radius * 5) ||
                     //                    maskpnt[SI] < (LE_pnt[SI] + approx_eye_radius * 5))
            )
            {
              mit.Set(face_rm);
            }
            //          else if (roiAutoInterpolator->Evaluate(maskpnt) != 1)
            //          {
            //            mit.Set(auto_roi_background);
            //          }
          }
          ++mit;
        }
      }
    }
  }

  // lambdas for the other two regions based on the mask values
  using MaskLabelToDistanceMapSeedFilter = typename itk::UnaryGeneratorImageFilter<MaskImageType, MaskImageType>;
  MaskLabelToDistanceMapSeedFilter::Pointer maskLabelToDistanceMapSeedFilter = MaskLabelToDistanceMapSeedFilter::New();
  maskLabelToDistanceMapSeedFilter->SetFunctor([&](const typename MaskImageType::PixelType input_mask_value) ->
                                               typename MaskImageType ::PixelType {
                                                 if (input_mask_value == valid_inside_pixel)
                                                 {
                                                   return static_cast<typename MaskImageType::PixelType>(0);
                                                 }
                                                 return static_cast<typename MaskImageType::PixelType>(1);
                                               });

  maskLabelToDistanceMapSeedFilter->SetInput(mask_labels);
  maskLabelToDistanceMapSeedFilter->Update();
  MaskImageType::Pointer binaryDistanceMapSeed = maskLabelToDistanceMapSeedFilter->GetOutput();

  using MaskLabelToNotFaceRegionFilter = typename itk::UnaryGeneratorImageFilter<MaskImageType, FadeMapType>;
  MaskLabelToNotFaceRegionFilter::Pointer maskLabelToNotFaceRegionFilter = MaskLabelToNotFaceRegionFilter::New();
  maskLabelToNotFaceRegionFilter->SetFunctor([&](const typename MaskImageType::PixelType input_mask_value) ->
                                             typename FadeMapType::PixelType {
                                               if (input_mask_value == valid_inside_pixel)
                                               {
                                                 return static_cast<typename MaskImageType::PixelType>(1.0);
                                               }
                                               return static_cast<typename MaskImageType::PixelType>(0.0);
                                             });

  maskLabelToNotFaceRegionFilter->SetInput(mask_labels);
  maskLabelToNotFaceRegionFilter->Update();
  FadeMapType::Pointer not_face_region = maskLabelToNotFaceRegionFilter->GetOutput();

  if (debugLevel >= 5)
  {
    std::cout << "Writing output distanceMapSeed" << std::endl;
    itkUtil::WriteImage<MaskImageType>(binaryDistanceMapSeed, outputDirectory + "/dist_map_seed.nii.gz");
  }

  // STEP 2: Deface the image values
  if (!noMaskApplication)
  {
    itk::NearestNeighborInterpolateImageFunction<MaskImageType>::Pointer maskInterpolator =
      itk::NearestNeighborInterpolateImageFunction<MaskImageType>::New();
    maskInterpolator->SetInputImage(mask_labels);

    // pad the distance map
    MaskImageType::SizeType lowerPadBound;
    lowerPadBound.Fill(80); // 40 mm @ 0.5 mm spacing
    MaskImageType::SizeType upperPadBound;
    upperPadBound.Fill(80); // 40 mm @ 0.5 mm spacing

    itk::ConstantPadImageFilter<MaskImageType, MaskImageType>::Pointer padImageFilter =
      itk::ConstantPadImageFilter<MaskImageType, MaskImageType>::New();
    padImageFilter->SetInput(binaryDistanceMapSeed);
    padImageFilter->SetPadLowerBound(lowerPadBound);
    padImageFilter->SetPadUpperBound(upperPadBound);
    padImageFilter->SetConstant(face_rm);
    padImageFilter->Update();
    if (debugLevel >= 5)
    {
      itkUtil::WriteImage<MaskImageType>(padImageFilter->GetOutput(), outputDirectory + "/padded.nii.gz");
    }
    // compute the distance map
    itk::SignedMaurerDistanceMapImageFilter<MaskImageType, InternalImageType>::Pointer signedDistanceMap =
      itk::SignedMaurerDistanceMapImageFilter<MaskImageType, InternalImageType>::New();
    signedDistanceMap->SetInput(padImageFilter->GetOutput());
    signedDistanceMap->SetInsideIsPositive(true);
    signedDistanceMap->Update();

    using DistBlurFilter = itk::SmoothingRecursiveGaussianImageFilter<InternalImageType, InternalImageType>;
    DistBlurFilter::Pointer blur = DistBlurFilter::New();
    blur->SetInput(signedDistanceMap->GetOutput());
    blur->SetSigma(2.25);
    blur->Update();
    InternalImageType::Pointer distanceMap = blur->GetOutput();

    for (itk::ImageRegionIterator<InternalImageType> diit(distanceMap, distanceMap->GetLargestPossibleRegion());
         !diit.IsAtEnd();
         ++diit)
    {
      const auto &                           refvalue = diit.Value();
      constexpr InternalImageType::PixelType allowed_blurring_area = 0.0; // Allow blurring if distance from edge
      if (refvalue < allowed_blurring_area)
      {
        diit.Set(0.0);
      }
      else
      {
        const InternalImageType::PixelType absrefvalue = std::abs(refvalue);
        const InternalImageType::PixelType distout = absrefvalue < 10.0 ? absrefvalue : 10.0;
        diit.Set(distout);
      }
    }


    if (debugLevel >= 2)
    {
      itkUtil::WriteImage<InternalImageType>(distanceMap, outputDirectory + "/dist_img.nii.gz");
    }
    itk::LinearInterpolateImageFunction<InternalImageType>::Pointer distanceMapInterpolator =
      itk::LinearInterpolateImageFunction<InternalImageType>::New();
    distanceMapInterpolator->SetInputImage(distanceMap);

    // blur configuration
    //  const int    numSigmas = 2;
    //  const double sigmas[numSigmas] = { 1.0, 1.1};
    constexpr size_t numSigmas = 13;
    // constexpr size_t lastSigmaIndex = numSigmas - 1;
    constexpr double sigmas[numSigmas] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };

    // Second add images that are also used for defacing
    for (const auto & passive_im_fn : passiveVolume)
    {
      img_list.emplace_back(itkUtil::ReadImage<InternalImageType>(passive_im_fn));
      img_filenames.emplace_back(passive_im_fn);
    }
    passiveVolume.clear();

    for (auto & curr_img : img_list)
    {
      // pre-compute different levels of gaussian blur
      InternalImageType::Pointer                                      blurredImages[numSigmas];
      itk::LinearInterpolateImageFunction<InternalImageType>::Pointer blurredImageInterpolators[numSigmas];
      using BlurFilter = itk::SmoothingRecursiveGaussianImageFilter<FadeMapType, FadeMapType>;
      for (size_t i = 0; i < numSigmas; ++i)
      {
        BlurFilter::Pointer blur = BlurFilter::New();
        blur->SetInput(curr_img);
        blur->SetSigma(sigmas[i]);
        blur->Update();
        blurredImages[i] = blur->GetOutput();
        blurredImageInterpolators[i] = itk::LinearInterpolateImageFunction<InternalImageType>::New();
        blurredImageInterpolators[i]->SetInputImage(blurredImages[i]);
        if (debugLevel >= 5)
        {
          itkUtil::WriteImage<InternalImageType>(blur->GetOutput(),
                                                 outputDirectory + "/blur_img_" + std::to_string(i) + ".nii.gz");
        }
      }

      const bool                                           use_zeros_face = (defaceMode == "zero");
      InternalImageType::PointType                         imgpnt;
      itk::ImageRegionIteratorWithIndex<InternalImageType> iit(curr_img, curr_img->GetLargestPossibleRegion());
      while (!iit.IsAtEnd())
      {
        const auto & curr_index{ iit.GetIndex() };
        curr_img->TransformIndexToPhysicalPoint<double>(curr_index, imgpnt);
        if (maskInterpolator->IsInsideBuffer(imgpnt) && distanceMapInterpolator->IsInsideBuffer(imgpnt))
        {
          const MaskImageType::PixelType mask_value = maskInterpolator->Evaluate(imgpnt);
          if (mask_value == valid_inside_pixel)
          {
            // Pass
          }
          else if (mask_value == face_rm || mask_value == auto_roi_background || mask_value == eye_boxes_code)
          {
            const InternalImageType::PixelType distanceValue = distanceMapInterpolator->Evaluate(imgpnt);
            // determine the correct blur value
            int              sigmaIndex = 0;
            constexpr double sigma_distance_ratio = 1.0; // Factor of sigma smoothing to distance ratio
            for (size_t i = 0; i < numSigmas; ++i)
            {
              if (sigmas[i] * sigma_distance_ratio < distanceValue)
              {
                sigmaIndex = i;
              }
            }
            if (blurredImageInterpolators[sigmaIndex]->IsInsideBuffer(imgpnt))
            {
              const InternalImageType::PixelType blurredValue{
                use_zeros_face
                  ? 0.0F
                  : static_cast<InternalImageType::PixelType>(blurredImageInterpolators[sigmaIndex]->Evaluate(imgpnt))
              };
              iit.Set(blurredValue);
              // iit.Set(sigmas[sigmaIndex]* 100);
            }
            else
            {
              iit.Set(0);
            }
          }
          else
          {
            iit.Set(0);
          }
        }
        else
        {
          iit.Set(0);
        }
        ++iit;
      }
    }



    // STEP 3 Generate intensity normalized images from the clippings.
    for (size_t i = 0; i < img_filenames.size(); ++i)
    {
      const auto input_fn = itksys::SystemTools::GetFilenameName(img_filenames[i]);

      const std::vector<std::string> output_fn_components{ outputDirectory, "/", input_fn };
      const std::string              output_fn = itksys::SystemTools::JoinPath(output_fn_components);
      using OutputImageType = itk::Image<unsigned short, 3>;
      const bool clip = !no_clip;
      const bool relative = !no_relative;

      auto intensity_normalized_output = brains_intensity_normalize_quantiles<InternalImageType, OutputImageType>(
        img_list[i], lowerPercentile, upperPercentile, lowerOutputIntensity, upperOutputIntensity, clip, relative);
      std::cout << "Writing output filename: " << output_fn << std::endl;
      itkUtil::WriteImage<OutputImageType>(intensity_normalized_output, output_fn);
    }
  }

  // Write out mask as the last element after used for images
  // If outputMask has path elements, then do not use outputDirectory.
  if (!outputMask.empty() && inputMask.empty())
  {
    if (outputMask.find('/') == std::string::npos)
    {
      const std::vector<std::string> output_fn_components{ outputDirectory, "/", outputMask };
      const std::string              output_fn = itksys::SystemTools::JoinPath(output_fn_components);
      outputMask = itksys::SystemTools::JoinPath(output_fn_components);
    }
    std::cout << "Writing output mask filename: " << outputMask << std::endl;
    itkUtil::WriteImage<MaskImageType>(mask_labels, outputMask);
  }
  return EXIT_SUCCESS;
}
