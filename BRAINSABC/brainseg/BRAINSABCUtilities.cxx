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
#include "BRAINSABCUtilities.h"
/*****************************
 * Now call the instantiations
 */
#include "BRAINSABCUtilities.hxx"
#include "LLSBiasCorrector.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkImageMaskSpatialObject.h"

template std::vector<FloatImageType::Pointer>
DuplicateImageList<FloatImageType>(const std::vector<FloatImageType::Pointer> &);

template std::vector<ShortImageType::Pointer>
DuplicateImageList<ShortImageType>(const std::vector<ShortImageType::Pointer> &);

template void
NormalizeProbListInPlace<FloatImageType>(std::vector<FloatImageType::Pointer> &);

template void
ZeroNegativeValuesInPlace<FloatImageType>(std::vector<FloatImageType::Pointer> &);

MapOfFloatImageVectors
ResampleImageListToFirstKeyImage(const std::string &            resamplerInterpolatorType,
                                 const MapOfFloatImageVectors & inputImageMap)
{
  muLogMacro(<< "Resampling input image map to the first key image." << std::endl);

  FloatImageType::ConstPointer KeyImageFirstRead = GetMapVectorFirstElement(inputImageMap).GetPointer();

  // Clear image list
  MapOfFloatImageVectors outputImageMap;

  // Resample the other images
  for (const auto & elem : inputImageMap)
  {
    auto         currImageIter = elem.second.begin();
    unsigned int i(0);
    while (currImageIter != elem.second.end())
    {
      FloatImageType::Pointer tmp = ResampleImageWithIdentityTransform<FloatImageType>(
        resamplerInterpolatorType, 0, (*currImageIter).GetPointer(), KeyImageFirstRead.GetPointer());
      // Add the image
      outputImageMap[elem.first].push_back(tmp);
      ++currImageIter;
      ++i;
    }
  }
  return outputImageMap;
}

ByteImageType::Pointer
intersectFOV(const MapOfMaskImageVectors & FOVMap)
{
  for (const auto & inputImageMapIter : FOVMap)
  {
    auto currModalIter = inputImageMapIter.second.begin();
    while (currModalIter != inputImageMapIter.second.end())
    {
    }
  }
  return nullptr;
}

// ByteImagePointer globalFeildOfView::Pointer = intersectFOV(intraSubjectRegisteredFOVMap);


ByteImageType::Pointer
ResampleToFirstImageList(const std::string &            resamplerInterpolatorType,
                         const MapOfFloatImageVectors & inputImageMap,
                         const MapOfTransformLists &    intraSubjectTransforms,
                         MapOfFloatImageVectors &       outputImageMap)
{
  outputImageMap.clear();
  muLogMacro(<< "ResampleToFirstImageList..." << std::endl);
  /*
   * This function, first, transforms all inputImageMap to the space of the first image of the map
   * using rigid transforms (intraSubjectTransforms) and Resampling InPlace interpolation.
   * Then, it resamples all images within one modality to the voxel lattice of the fist image of that modality channel
   * using resamplerInterpolatorType and Identity transform.
   */
  MapOfFloatImageVectors resampledInPlaceMap;
  MapOfFloatImageVectors resampledToFirstImageMap;
  PrintMapOfImageVectors(inputImageMap);

  // ResampleInPlace all images to the physical space of the first image
  //
  for (auto inputImageMapIter = inputImageMap.begin(); inputImageMapIter != inputImageMap.end(); ++inputImageMapIter)
  {
    auto         currModalIter = inputImageMapIter->second.begin();
    unsigned int i(0);
    auto         xfrmIt = intraSubjectTransforms.at(inputImageMapIter->first).begin();
    while (currModalIter != inputImageMapIter->second.end())
    {
      muLogMacro(<< "Resampling input image " << inputImageMapIter->first << " #" << i
                 << " to the physical space of the first image." << std::endl);
      using ResampleIPFilterType = itk::ResampleInPlaceImageFilter<FloatImageType, FloatImageType>;
      using ResampleIPFilterPointer = ResampleIPFilterType::Pointer;

      using VersorRigid3DTransformType = itk::VersorRigid3DTransform<double>;
      const VersorRigid3DTransformType::ConstPointer tempRigidTransform =
        dynamic_cast<VersorRigid3DTransformType const *>((*xfrmIt).GetPointer());
      if (tempRigidTransform.IsNull())
      {
        std::cerr << "Error in type conversion. " << __FILE__ << __LINE__ << std::endl;
        std::cerr << "ResampleInPlace is only allowed with rigid transform type." << std::endl;
        throw;
      }

      ResampleIPFilterPointer resampleIPFilter = ResampleIPFilterType::New();
      resampleIPFilter->SetInputImage((*currModalIter));
      resampleIPFilter->SetRigidTransform(tempRigidTransform);
      resampleIPFilter->Update();
      FloatImageType::Pointer tmp = resampleIPFilter->GetOutput();

      // Add the image
      resampledInPlaceMap[inputImageMapIter->first].push_back(tmp);
      ++currModalIter;
      ++xfrmIt;
      ++i;
    }
  }

  // Resample each intra subject image to the first image of its modality
  //
  muLogMacro(<< "Resampling each intra-subject image to the first image of the primary modality using "
             << resamplerInterpolatorType << " interpolation." << std::endl);
  // Set the outsideFOV value to 100x the max floating point value.  This allows for some numerical computations to
  // occur in the resampling, but also allows identifcation of very large values.
  // constexpr FloatImageType::PixelType outsideFOVCode = std::numeric_limits<FloatImageType::PixelType>::max()/100.0;

  // Use the first image of the vector from the ordered modality preference as the key image for all others.
  FloatImageType::Pointer allModalityKeySubjectImage = resampledInPlaceMap.cbegin()->second.cbegin()->GetPointer();
  ByteImageType::Pointer  intraSubjectFOVIntersectionMask = ByteImageType::New();
  intraSubjectFOVIntersectionMask->CopyInformation(allModalityKeySubjectImage.GetPointer());
  intraSubjectFOVIntersectionMask->SetRegions(allModalityKeySubjectImage.GetPointer()->GetLargestPossibleRegion());
  intraSubjectFOVIntersectionMask->Allocate();
  intraSubjectFOVIntersectionMask->FillBuffer(1); // Initialize so that entire FOV is included.


  for (auto & modalityMapEntry : resampledInPlaceMap)
  {
    auto currModalIter = modalityMapEntry.second.begin();
    while (currModalIter != modalityMapEntry.second.end())
    {
      using MinMaxType = itk::MinimumMaximumImageCalculator<FloatImageType>;
      MinMaxType::Pointer minmaxcalc = MinMaxType::New();
      minmaxcalc->SetImage((*currModalIter).GetPointer());
      minmaxcalc->ComputeMinimum();
      minmaxcalc->Compute();

      const FloatImageType::PixelType     minValue = minmaxcalc->GetMinimum();
      constexpr FloatImageType::PixelType min_background_offset_value = 8.00;
      const FloatImageType::PixelType     background_threshold_value = minValue - min_background_offset_value;
      // Set background offset so background fill is smaller than
      FloatImageType::Pointer resampledToFirstImage =
        ResampleImageWithIdentityTransform<FloatImageType>(resamplerInterpolatorType,
                                                           background_threshold_value,
                                                           (*currModalIter).GetPointer(),
                                                           allModalityKeySubjectImage.GetPointer());

      using OtsuThresholdType = itk::OtsuThresholdImageFilter<FloatImageType, ByteImageType>;
      OtsuThresholdType::Pointer otsufilter = OtsuThresholdType::New();
      otsufilter->SetInput(resampledToFirstImage);
      otsufilter->SetInsideValue(0);
      otsufilter->SetOutsideValue(1);
      otsufilter->Update();

      using BinaryThresholdType = itk::BinaryThresholdImageFilter<FloatImageType, ByteImageType>;
      BinaryThresholdType::Pointer bthresh = BinaryThresholdType::New();
      bthresh->SetInput(resampledToFirstImage);
      // Use 1/2 the otsu threshold to get expanded background
      bthresh->SetLowerThreshold(otsufilter->GetThreshold() * 0.5f);
      bthresh->Update();

      using MaskType = itk::ImageMaskSpatialObject<FloatImageType::ImageDimension>;

      MaskType::Pointer spatialObjectMask = MaskType::New();
      spatialObjectMask->SetImage(bthresh->GetOutput());
      spatialObjectMask->Update();
      const itk::SpatialObject<FloatImageType::ImageDimension>::BoundingBoxType * bb =
        spatialObjectMask->GetMyBoundingBoxInWorldSpace();

      // Zero the mask region outside FOV and also the intensities with outside FOV code
      using InternalIteratorType = itk::ImageRegionConstIterator<FloatImageType>;
      InternalIteratorType resampledToFirstIt(resampledToFirstImage, resampledToFirstImage->GetLargestPossibleRegion());

      using MaskIteratorType = itk::ImageRegionIteratorWithIndex<ByteImageType>;
      MaskIteratorType maskIt(intraSubjectFOVIntersectionMask,
                              intraSubjectFOVIntersectionMask->GetLargestPossibleRegion());
      maskIt.GoToBegin();
      resampledToFirstIt.GoToBegin();
      // Any value larger that an unsigned 32 bit integer is too large and assumed .
      // constexpr FloatImageType::PixelType halfoutsideFOVCode = 8388608; // 2^23 The mantissa of a floating point
      // number.
      while (!maskIt.IsAtEnd())
      {
        using PT = MaskType::PointType;
        const PT thisPoint = resampledToFirstImage->TransformIndexToPhysicalPoint<float>(maskIt.GetIndex());
        if (resampledToFirstIt.Value() <= background_threshold_value ||
            (!bb->IsInside(thisPoint))) // Voxel came from outside
        // the original FOV during registration, so invalidate it.
        {
          maskIt.Set(0); // Set it as an invalid voxel in
        }
        ++maskIt;
        ++resampledToFirstIt;
      }
      ++currModalIter;
      resampledToFirstImageMap[modalityMapEntry.first].push_back(resampledToFirstImage);
    }
  }

  //  //Fill holes filter
  //  using HoleFillType = itk::BinaryFillholeImageFilter<ByteImageType>;
  //  HoleFillType::Pointer hole_fill_filter = HoleFillType::New();
  //  hole_fill_filter->SetInput(intraSubjectFOVIntersectionMask);
  //  hole_fill_filter->FullyConnectedOn();
  //  hole_fill_filter->SetForegroundValue(1);
  //  hole_fill_filter->Update();
  //  ByteImageType::Pointer intraSubjectFOVIntersectionMaskFilled = hole_fill_filter->GetOutput();
  //  intraSubjectFOVIntersectionMask=nullptr;

  // Now zero out all values not in the FOV mask for all resampled to first images
  for (auto & modalityToFirstMapEntry : resampledToFirstImageMap)
  {
    auto currModalToFirstIter = modalityToFirstMapEntry.second.begin();
    while (currModalToFirstIter != modalityToFirstMapEntry.second.end())
    {
      // Zero the mask region outside FOV and also the intensities with outside FOV code
      using InternalIteratorType = itk::ImageRegionIterator<FloatImageType>;
      InternalIteratorType resampledToFirstIt(*currModalToFirstIter,
                                              (*currModalToFirstIter)->GetLargestPossibleRegion());

      using MaskIteratorType = itk::ImageRegionConstIterator<ByteImageType>;
      MaskIteratorType maskIt(intraSubjectFOVIntersectionMask,
                              intraSubjectFOVIntersectionMask->GetLargestPossibleRegion());
      maskIt.GoToBegin();
      resampledToFirstIt.GoToBegin();
      while (!maskIt.IsAtEnd())
      {
        if (maskIt.Value() < 1)
        {
          // Voxel came from outside
          // the original FOV during registration, so invalidate it.
          resampledToFirstIt.Set(0);
        }
        ++maskIt;
        ++resampledToFirstIt;
      }
      outputImageMap[modalityToFirstMapEntry.first].push_back(*currModalToFirstIter);
      ++currModalToFirstIter;
    }
  }
  return intraSubjectFOVIntersectionMask;
}
