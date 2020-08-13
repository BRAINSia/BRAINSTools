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
#include "itkMaskImageFilter.h"

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
  PARSE_ARGS; // <-- this is not really used yet TODO: Make argument parsing SEM compatible.


  // STEP 1: Using all the anatomical images for this session, generate a label mask defining
  // regions that are 0= in all the images
  //                  1= part of the face region
  //                  2= 80mm below the ac point --> altered to remove this region by resampling
  //                  3= Outside FOV for at least 1 image
  using InternalImageType = itk::Image<float, 3>;

  std::vector<InternalImageType::Pointer> img_list;
  for (const auto & im_fn : inputVolume)
  {
    img_list.emplace_back(itkUtil::ReadImage<InternalImageType>(im_fn));
  }

  using MaskImageType = itk::Image<unsigned char, 3>;

  // Make a 1/2 mm masked volume
  MaskImageType::Pointer mask_labels = MaskImageType::New();
  {
    const MaskImageType::RegionType::IndexType starting_index{ 0, 0, 0 };
    const MaskImageType::RegionType::SizeType  img_size{ 512, 512, 352 }; // 512 - 160 (80mm @ 0.5 spacing)
    MaskImageType::RegionType                  region;
    region.SetSize(img_size);
    region.SetIndex(starting_index);

    mask_labels->SetRegions(region);

    MaskImageType::DirectionType id;
    id.SetIdentity();
    mask_labels->SetDirection(id);

    const double spacing[3] = { 0.5, 0.5, 0.5 };
    mask_labels->SetSpacing(spacing);


    // Empirically chosen origin voxel based on standard center of brain
    // -128.0, -114.0, -132.0
    const double origin[3] = { -128.0, -114.0, -52.0 }; // -132 + 80 removing 80mm inferior
    mask_labels->SetOrigin(origin);
    mask_labels->Allocate(true);
  }
  using FadeMapType = itk::Image<float, 3>;
  FadeMapType::Pointer not_face_region = FadeMapType::New();
  not_face_region->CopyInformation(mask_labels);
  not_face_region->SetRegions(mask_labels->GetLargestPossibleRegion());
  not_face_region->Allocate();
  not_face_region->FillBuffer(1.0);

  const auto               lmks = ReadSlicer3toITKLmk(inputLandmarks);
//  MaskImageType::PointType AC_pnt = lmks.at("AC");
  MaskImageType::PointType RE_pnt = lmks.at("RE");
  MaskImageType::PointType LE_pnt = lmks.at("LE");

  // constexpr size_t LR=0;
  constexpr size_t PA = 1;
  constexpr size_t SI = 2;

  constexpr int valid_inside_pixel = 0;
  constexpr int face_rm = 1;
  // constexpr int below_ac = 2;
  constexpr int outside_fov = 3;
  // Find all the out of FOV spaces fro the mask image
  MaskImageType::PointType     maskpnt;
  InternalImageType::IndexType imgindex;
  for (auto & curr_img : img_list)
  {
    {
      itk::ImageRegionIteratorWithIndex<MaskImageType> mit(mask_labels, mask_labels->GetLargestPossibleRegion());
      itk::ImageRegionIteratorWithIndex<FadeMapType> fit(not_face_region, not_face_region->GetLargestPossibleRegion());
      while (!mit.IsAtEnd())
      {
        mask_labels->TransformIndexToPhysicalPoint<double>(mit.GetIndex(), maskpnt);

        const bool isInside = curr_img->TransformPhysicalPointToIndex(maskpnt, imgindex);
        if (!isInside)
        {
          mit.Set(outside_fov);
          fit.Set(0.0);
        }
//        else if (maskpnt[SI] < AC_pnt[SI] - 80.0) // removing 80 mm inferior
//        {
//          mit.Set(below_ac);
//          fit.Set(0.0);
//        }
        else if ((maskpnt[PA] < RE_pnt[PA] || maskpnt[PA] < LE_pnt[PA]) &&
                 (maskpnt[SI] < RE_pnt[SI] || maskpnt[SI] < LE_pnt[SI]))
        {
          mit.Set(face_rm);
          fit.Set(0.0);
        }
        ++mit;
        ++fit;
      }
    }
  }

  // If outputMask has path elements, then do not use outputDirectory.
  if (outputMask.find('/') == std::string::npos)
  {
    const std::vector<std::string> output_fn_components{ outputDirectory, "/", outputMask };
    const std::string              output_fn = itksys::SystemTools::JoinPath(output_fn_components);
    ;
    outputMask = itksys::SystemTools::JoinPath(output_fn_components);
  }
  std::cout << "Writing output mask filename: " << outputMask << std::endl;
  itkUtil::WriteImage<MaskImageType>(mask_labels, outputMask);


  // STEP 2: Deface the image values
  const bool           doBlur = (defaceMode == "blur");
  FadeMapType::Pointer blur_image = not_face_region;
  if (doBlur)
  {
    // This creates a gradual fading rather than a sudden edge.
    using BlurFilter = itk::SmoothingRecursiveGaussianImageFilter<FadeMapType, FadeMapType>;
    BlurFilter::Pointer blur = BlurFilter::New();
    blur->SetInput(not_face_region);
    blur->SetSigma(3);                    // 3mm smoothing
    blur->SetNormalizeAcrossScale(false); // If true, negative values can result
    blur->Update();
    blur_image = blur->GetOutput();
  }
  // DEBUG
  // itkUtil::WriteImage<FadeMapType>(blur_image,"/tmp/blur_img.nii.gz");

  itk::NearestNeighborInterpolateImageFunction<MaskImageType>::Pointer interp =
    itk::NearestNeighborInterpolateImageFunction<MaskImageType>::New();
  interp->SetInputImage(mask_labels);
  InternalImageType::PointType imgpnt;
  for (auto & curr_img : img_list)
  {
    itk::ImageRegionIteratorWithIndex<InternalImageType> iit(curr_img, curr_img->GetLargestPossibleRegion());
    while (!iit.IsAtEnd())
    {
      const auto & curr_index{ iit.GetIndex() };
      curr_img->TransformIndexToPhysicalPoint<double>(curr_index, imgpnt);
      if (interp->IsInsideBuffer(imgpnt))
      {
        const MaskImageType::PixelType mask_value = interp->Evaluate(imgpnt);
        if (mask_value == valid_inside_pixel)
        {
          // Pass
        }
        else if (mask_value == face_rm)
        {
          const auto blur_image_value = blur_image->GetPixel(curr_index);
          const auto curr_value = static_cast<FadeMapType::PixelType>(iit.Get());
          const auto new_value = static_cast<InternalImageType::PixelType>(blur_image_value * curr_value);
          iit.Set(new_value);
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
  for (size_t i = 0; i < inputVolume.size(); ++i)
  {
    const auto input_fn = itksys::SystemTools::GetFilenameName(inputVolume[i]);

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
  return EXIT_SUCCESS;
}
