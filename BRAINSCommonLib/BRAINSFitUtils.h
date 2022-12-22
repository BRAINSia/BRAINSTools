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
#ifndef __BRAINSFitUtils_h
#define __BRAINSFitUtils_h

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageDuplicator.h"

#include "itkScaleVersor3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkBSplineTransform.h"
#include "itkBRAINSROIAutoImageFilter.h"

#include "algorithm"
#include "math.h"
#include <iostream>
#include <vcl_compiler.h>

/**
 * This file contains utility functions that are common to a few of the
 *BRAINSFit Programs.
 */

/**
 * Boilerplate conversion to get a safe reference to the internal Image stored in a
 * ImageMaskSpatialObjectType
 */
extern itk::Image<unsigned char, 3>::ConstPointer
ExtractConstPointerToImageMaskFromImageSpatialObject(const SpatialObjectType::ConstPointer & inputSpatialObject);

extern itk::ImageMaskSpatialObject<3>::ConstPointer
ConvertMaskImageToSpatialMask(const itk::Image<unsigned char, 3>::ConstPointer & inputImage);

template <typename TransformType, unsigned int VImageDimension>
void
DoCenteredTransformMaskClipping(ImageMaskPointer &              fixedMask,
                                ImageMaskPointer &              movingMask,
                                typename TransformType::Pointer transform,
                                double                          maskInferiorCutOffFromCenter)
{
  if (fixedMask.IsNull() || movingMask.IsNull())
  {
    return;
  }
  if (maskInferiorCutOffFromCenter >= 1000.0)
  {
    return;
  }
  std::cerr << "maskInferiorCutOffFromCenter is " << maskInferiorCutOffFromCenter << std::endl;

  typename TransformType::InputPointType   rotationCenter = transform->GetCenter();
  typename TransformType::OutputVectorType translationVector = transform->GetTranslation();

  using MaskImageType = itk::Image<unsigned char, 3>;

  typename MaskImageType::PointType fixedCenter;
  typename MaskImageType::PointType movingCenter;
  for (unsigned int i = 0; i < VImageDimension; ++i)
  {
    fixedCenter[i] = rotationCenter[i];
    movingCenter[i] = translationVector[i] - rotationCenter[i];
  }

  typename MaskImageType::Pointer fixedMaskImage = nullptr;
  {
    const typename MaskImageType::ConstPointer tempOutputFixedVolumeROI =
      ExtractConstPointerToImageMaskFromImageSpatialObject(fixedMask.GetPointer());
    using DuplicatorType = itk::ImageDuplicator<MaskImageType>;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(tempOutputFixedVolumeROI);
    duplicator->Update();
    fixedMaskImage = duplicator->GetOutput();
  }
  typename MaskImageType::Pointer movingMaskImage = nullptr;
  {
    const typename MaskImageType::ConstPointer tempOutputMovingVolumeROI =
      ExtractConstPointerToImageMaskFromImageSpatialObject(movingMask.GetPointer());
    using DuplicatorType = itk::ImageDuplicator<MaskImageType>;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(tempOutputMovingVolumeROI);
    duplicator->Update();
    movingMaskImage = duplicator->GetOutput();
  }

  typename MaskImageType::PointType fixedInferior = fixedCenter;
  typename MaskImageType::PointType movingInferior = movingCenter;

  fixedInferior[2] -= maskInferiorCutOffFromCenter;  // negative because
                                                     // Superior is large in
                                                     // magnitude.
  movingInferior[2] -= maskInferiorCutOffFromCenter; // ITK works in an LPS
                                                     // system.

  //  Here we will set the appropriate parts of the f/m MaskImages to zeros....
  typename MaskImageType::PixelType zero = 0;
  typename MaskImageType::PointType location;
  using MaskIteratorType = itk::ImageRegionIteratorWithIndex<MaskImageType>;

  MaskIteratorType fixedIter(fixedMaskImage, fixedMaskImage->GetLargestPossibleRegion());
  fixedIter.GoToBegin();
  while (!fixedIter.IsAtEnd())
  {
    fixedMaskImage->TransformIndexToPhysicalPoint(fixedIter.GetIndex(), location);
    if (location[2] < fixedInferior[2])
    {
      fixedIter.Set(zero);
    }
    ++fixedIter;
  }

  MaskIteratorType movingIter(movingMaskImage, movingMaskImage->GetLargestPossibleRegion());
  movingIter.GoToBegin();
  while (!movingIter.IsAtEnd())
  {
    movingMaskImage->TransformIndexToPhysicalPoint(movingIter.GetIndex(), location);
    if (location[2] < movingInferior[2])
    {
      movingIter.Set(zero);
    }
    ++movingIter;
  }

  using ImageMaskSpatialObjectType = itk::ImageMaskSpatialObject<MaskImageType::ImageDimension>;

  typename ImageMaskSpatialObjectType::Pointer fixedMaskSpatialObject = ImageMaskSpatialObjectType::New();
  fixedMaskSpatialObject->SetImage(fixedMaskImage);
  fixedMaskSpatialObject->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()
  fixedMask = fixedMaskSpatialObject.GetPointer();

  typename ImageMaskSpatialObjectType::Pointer movingMaskSpatialObject = ImageMaskSpatialObjectType::New();
  movingMaskSpatialObject->SetImage(movingMaskImage);
  movingMaskSpatialObject->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()
  movingMask = movingMaskSpatialObject.GetPointer();
}

/**
 * The ComputeRobustMinMaxMean function will identify the 10th and 90th quantile
 * an extropolate an expected min/max values (i.e the 0th, and 100th quantile) from
 * that estimate.  The value returned for the true min/max values are
 * the min(estimate 0thquantile, actualmin)/max(estimated 100thquantile, actualmax);
 *
 * This is a robust way to remove small numbers of high spike noise from artificially
 * expanding the range of the image.  It helps to ensure that the middle 80% of the image
 * occupies a sufficient part of the range computed.
 *
 */
template <typename TInputImage, typename TMaskImage>
void
ComputeRobustMinMaxMean(const float Qalpha, // Remove 1% from computations by setting Qalpha=0.005
                        typename TInputImage::ConstPointer image,
                        typename TMaskImage::ConstPointer  mask,
                        float &                            minValue, // INFO:  Make this into
                                          // itk::NumericTraits<TInputImage::PixelType>::RealType;
                        float & maxValue, // INFO:  Make this into
                                          // itk::NumericTraits<TInputImage::PixelType>::RealType;
                        float & meanValue // INFO:  Make this into
                                          // itk::NumericTraits<TInputImage::PixelType>::RealType;
)
{
  // This is a more stable way of determining the range of values that the image
  // has.
  // By eliminating possible "bright or dark" noise in the image.
  minValue = std::numeric_limits<float>::max();
  maxValue = std::numeric_limits<float>::min();
  const auto num_pixels = image->GetBufferedRegion().GetNumberOfPixels();
  std::vector<typename TInputImage::PixelType> fixedList;
  fixedList.reserve(num_pixels);
  {
    itk::ImageRegionConstIteratorWithIndex<TInputImage> fi(image, image->GetBufferedRegion());
    while (!fi.IsAtEnd())
    {
      typename TInputImage::PointType physicalPoint;
      image->TransformIndexToPhysicalPoint(fi.GetIndex(), physicalPoint);

      bool inCaluationRegion = true;
      if (mask.IsNotNull() && (!mask->IsInsideInWorldSpace(physicalPoint)))
      // A null mask implies entire space is to be used.
      {
        inCaluationRegion = false;
      }
      if (inCaluationRegion)
      {
        const typename TInputImage::PixelType currValue = fi.Get();
        minValue = std::min(minValue, currValue);
        maxValue = std::max(maxValue, currValue);
        fixedList.push_back(currValue);
      }
      ++fi;
    }
  }
  std::sort(fixedList.begin(), fixedList.end());

  // Make sure that center 1/2 (25%-75%) of intensity values spans center 1/2 of
  // histogram
  // Compute extended line through these two Quantile points
  const float LQx = Qalpha;
  const float HQx = 1.0 - Qalpha;
  const auto list_size = fixedList.size();
  const auto index_lowest_quantile_value = static_cast<size_t>(list_size * LQx);
  const float fixedLQy = fixedList[index_lowest_quantile_value];
  const auto index_highest_quantile_value = static_cast<size_t>(list_size * HQx);
  const float fixedHQy = fixedList[index_highest_quantile_value];
  const float fixedQSlope = (fixedHQy - fixedLQy) / (HQx - LQx);
  const float fixedZeroQy = fixedLQy - fixedQSlope * LQx;
  const float fixedOneQy = fixedQSlope * 1.0 + fixedZeroQy;

  std::cout << " Quantile Fixed Points " << Qalpha << ": " << fixedLQy << " " << fixedHQy << " slope: " << fixedQSlope
            << std::endl;
  std::cout << " Quantile Range"
            << ": " << fixedZeroQy << " " << fixedOneQy << std::endl;
  std::cout << "PreFix Range: " << minValue << " " << maxValue << std::endl;

  // Now set to the range of values based on linear extrapolation of the
  // quantiles
  minValue = std::max(fixedZeroQy, minValue);
  maxValue = std::min(fixedOneQy, maxValue);
  std::cout << "PostFix Range: " << minValue << " " << maxValue << std::endl;

  { // For all voxels in valid range, compute the mean.
    // This assumes that other values are noise values, and that their values
    // have no meaning.
    float  sum = 0.0F;
    size_t counter = 0;
    for (typename std::vector<typename TInputImage::PixelType>::const_iterator it = fixedList.begin();
         it != fixedList.end();
         ++it)
    {
      const auto value = static_cast<float>(*it);
      if (value <= maxValue && value >= minValue)
      {
        sum += value;
        ++counter;
      }
    }
    meanValue = sum / static_cast<float>(counter);
  }
}

/**
 *This function produces an image where the very high and low tails of the image are clamped
 */
template <typename TInputImage, typename TMaskSpatialObject>
typename TInputImage::Pointer
ClampNoisyTailsOfImage(const float                               m_RemoveIntensityOutliers,
                       typename TInputImage::ConstPointer        InputImage,
                       typename TMaskSpatialObject::ConstPointer mask)
{
  using DuplicatorType = itk::ImageDuplicator<TInputImage>;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(InputImage);
  duplicator->Update();
  typename TInputImage::Pointer image = duplicator->GetOutput();

  float min = NAN;
  float max = NAN;
  float mean = NAN;
  ComputeRobustMinMaxMean<TInputImage, TMaskSpatialObject>(
    m_RemoveIntensityOutliers, image.GetPointer(), mask.GetPointer(), min, max, mean);
  itk::ImageRegionIterator<TInputImage> fi(image, image->GetBufferedRegion());
  while (!fi.IsAtEnd())
  {
    if (fi.Value() > max)
    {
      fi.Set(max);
    }
    if (fi.Value() < min)
    {
      fi.Set(min);
    }
    ++fi;
  }

  return image;
}

#endif // __BRAINSFITUTILS_h
