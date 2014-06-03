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

#include "vcl_algorithm.h"

/**
  * This file contains utility functions that are common to a few of the
  *BRAINSFit Programs.
  */

static const unsigned int BFNSSpaceDimension = 3;
static const unsigned int BFNSplineOrder = 3;
typedef double CoordinateRepType;
typedef itk::BSplineTransform<
    CoordinateRepType,
    BFNSSpaceDimension,
    BFNSplineOrder>                                        BSplineTransformType;

typedef itk::VersorRigid3DTransform<double>              VersorRigid3DTransformType;
typedef itk::ScaleVersor3DTransform<double>              ScaleVersor3DTransformType;
typedef itk::ScaleSkewVersor3DTransform<double>          ScaleSkewVersor3DTransformType;
typedef itk::AffineTransform<double, BFNSSpaceDimension> AffineTransformType;

typedef itk::SpatialObject<3>                                      SpatialObjectType;
typedef SpatialObjectType::Pointer                                 ImageMaskPointer;
typedef itk::Image<unsigned char, 3>                               MaskImageType;
typedef itk::ImageMaskSpatialObject<MaskImageType::ImageDimension> ImageMaskSpatialObjectType;

/**
 * Boilerplate conversion to get a safe reference to the internal Image stored in a
 * ImageMaskSpatialObjectType
 */
extern MaskImageType::ConstPointer ExtractConstPointerToImageMaskFromImageSpatialObject(
  SpatialObjectType::ConstPointer inputSpatialObject);

extern SpatialObjectType::ConstPointer ConvertMaskImageToSpatialMask(
  MaskImageType::ConstPointer inputImage );

template <class TransformType, unsigned int VImageDimension>
void DoCenteredTransformMaskClipping(
  ImageMaskPointer & fixedMask,
  ImageMaskPointer & movingMask,
  typename TransformType::Pointer transform,
  double maskInferiorCutOffFromCenter)
{
  if( fixedMask.IsNull()  ||  movingMask.IsNull() )
    {
    return;
    }
  if( maskInferiorCutOffFromCenter >= 1000.0 )
    {
    return;
    }
  std::cerr << "maskInferiorCutOffFromCenter is " << maskInferiorCutOffFromCenter << std::endl;

  typedef unsigned char PixelType;

  typename TransformType::InputPointType rotationCenter = transform->GetCenter();
  typename TransformType::OutputVectorType translationVector = transform->GetTranslation();

  typename MaskImageType::PointType fixedCenter;
  typename MaskImageType::PointType movingCenter;
  for( unsigned int i = 0; i < VImageDimension; ++i )
    {
    fixedCenter[i]  = rotationCenter[i];
    movingCenter[i] = translationVector[i] - rotationCenter[i];
    }

  typename MaskImageType::Pointer fixedMaskImage  = NULL;
    {
    const typename MaskImageType::ConstPointer tempOutputFixedVolumeROI =
      ExtractConstPointerToImageMaskFromImageSpatialObject(fixedMask.GetPointer() );
    typedef itk::ImageDuplicator<MaskImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(tempOutputFixedVolumeROI);
    duplicator->Update();
    fixedMaskImage = duplicator->GetModifiableOutput();
    }
  typename MaskImageType::Pointer movingMaskImage = NULL;
    {
    const typename MaskImageType::ConstPointer tempOutputMovingVolumeROI =
      ExtractConstPointerToImageMaskFromImageSpatialObject(movingMask.GetPointer() );
    typedef itk::ImageDuplicator<MaskImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(tempOutputMovingVolumeROI);
    duplicator->Update();
    movingMaskImage = duplicator->GetModifiableOutput();
    }

  typename MaskImageType::PointType fixedInferior  = fixedCenter;
  typename MaskImageType::PointType movingInferior = movingCenter;

  fixedInferior[2] -= maskInferiorCutOffFromCenter;   // negative because
                                                      // Superior is large in
                                                      // magnitude.
  movingInferior[2] -= maskInferiorCutOffFromCenter;  // ITK works in an LPS
                                                      // system.

  //  Here we will set the appropriate parts of the f/m MaskImages to zeros....
  typename MaskImageType::PixelType zero = 0;
  typename MaskImageType::PointType location;
  typedef itk::ImageRegionIteratorWithIndex<MaskImageType> MaskIteratorType;

  MaskIteratorType fixedIter( fixedMaskImage, fixedMaskImage->GetLargestPossibleRegion() );
  fixedIter.GoToBegin();
  while( !fixedIter.IsAtEnd() )
    {
    fixedMaskImage->TransformIndexToPhysicalPoint(fixedIter.GetIndex(), location);
    if( location[2] < fixedInferior[2] )
      {
      fixedIter.Set(zero);
      }
    ++fixedIter;
    }

  MaskIteratorType movingIter( movingMaskImage, movingMaskImage->GetLargestPossibleRegion() );
  movingIter.GoToBegin();
  while( !movingIter.IsAtEnd() )
    {
    movingMaskImage->TransformIndexToPhysicalPoint(movingIter.GetIndex(), location);
    if( location[2] < movingInferior[2] )
      {
      movingIter.Set(zero);
      }
    ++movingIter;
    }

  typename ImageMaskSpatialObjectType::Pointer  fixedMaskSpatialObject = ImageMaskSpatialObjectType::New();
  fixedMaskSpatialObject->SetImage(fixedMaskImage);
  fixedMaskSpatialObject->ComputeObjectToWorldTransform();
  fixedMask = fixedMaskSpatialObject.GetPointer();

  typename ImageMaskSpatialObjectType::Pointer  movingMaskSpatialObject = ImageMaskSpatialObjectType::New();
  movingMaskSpatialObject->SetImage(movingMaskImage);
  movingMaskSpatialObject->ComputeObjectToWorldTransform();
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
template <class TInputImage, class TMaskImage>
void ComputeRobustMinMaxMean(
  const float Qalpha, // Remove 1% from computations by setting Qalpha=0.005
  typename TInputImage::ConstPointer image,
  typename TMaskImage::ConstPointer mask,
  float & minValue, // TODO:  Make this into
                    // itk::NumericTraits<TInputImage::PixelType>::RealType;
  float & maxValue, // TODO:  Make this into
                    // itk::NumericTraits<TInputImage::PixelType>::RealType;
  float & meanValue // TODO:  Make this into
                    // itk::NumericTraits<TInputImage::PixelType>::RealType;
  )
{
  // This is a more stable way of determining the range of values that the image
  // has.
  // By eliminating possible "bright or dark" noise in the image.
  minValue = vcl_numeric_limits<float>::max();
  maxValue = vcl_numeric_limits<float>::min();
  std::vector<typename TInputImage::PixelType> fixedList(image->GetBufferedRegion().GetNumberOfPixels() );
    {
    itk::ImageRegionConstIteratorWithIndex<TInputImage> fi(image, image->GetBufferedRegion() );
    while( !fi.IsAtEnd() )
      {
      typename TInputImage::PointType physicalPoint;
      image->TransformIndexToPhysicalPoint(fi.GetIndex(), physicalPoint);

      bool inCaluationRegion = true;
      if( mask.IsNotNull() &&  ( !mask->IsInside(physicalPoint) ) )  // A null
                                                                     // mask
                                                                     // implies
                                                                     // entire
                                                                     // space is
                                                                     // to be
                                                                     // used.
        {
        inCaluationRegion = false;
        }
      if( inCaluationRegion )
        {
        const typename TInputImage::PixelType currValue = fi.Get();
        minValue = vcl_min(minValue, currValue);
        maxValue = vcl_max(maxValue, currValue);
        fixedList.push_back(currValue);
        }
      ++fi;
      }
    }
  std::sort(fixedList.begin(), fixedList.end() );

  // Make sure that center 1/2 (25%-75%) of intensity values spans center 1/2 of
  // histogram
  // Compute extended line through these two Quantile points
  const float LQx = Qalpha;
  const float HQx = 1.0 - Qalpha;

  const float fixedLQy = fixedList[static_cast<size_t>(fixedList.size() * LQx)];
  const float fixedHQy = fixedList[static_cast<size_t>(fixedList.size() * HQx)];
  const float fixedQSlope = (fixedHQy - fixedLQy) / (HQx - LQx);
  const float fixedZeroQy = fixedLQy - fixedQSlope * LQx;
  const float fixedOneQy = fixedQSlope * 1.0 + fixedZeroQy;

  std::cout << " Quantile Fixed Points " << Qalpha << ": " << fixedLQy << " " << fixedHQy << " slope: "
            << fixedQSlope <<  std::endl;
  std::cout << " Quantile Range" << ": " << fixedZeroQy << " " << fixedOneQy << std::endl;
  std::cout << "PreFix Range: " << minValue   << " " << maxValue << std::endl;

  // Now set to the range of values based on linear extrapolation of the
  // quantiles
  minValue = vcl_max(fixedZeroQy, minValue);
  maxValue = vcl_min(fixedOneQy, maxValue);
  std::cout << "PostFix Range: " << minValue   << " " << maxValue << std::endl;

    { // For all voxels in valid range, compute the mean.
      // This assumes that other values are noise values, and that their values
      // have no meaning.
    float  sum = 0.0F;
    size_t counter = 0;
    for( typename std::vector<typename TInputImage::PixelType>::const_iterator it = fixedList.begin();
         it != fixedList.end(); ++it )
      {
      const float value = static_cast<float>(*it);
      if( value <= maxValue && value >= minValue )
        {
        sum += value;
        counter++;
        }
      }
    meanValue = sum / static_cast<float>(counter);
    }
}

/**
 *This function produces an image where the very high and low tails of the image are clamped
 */
template <class TInputImage, class TMaskSpatialObject>
typename TInputImage::Pointer ClampNoisyTailsOfImage(
  const float m_RemoveIntensityOutliers,
  typename TInputImage::ConstPointer InputImage,
  typename TMaskSpatialObject::ConstPointer mask)
{
  typedef itk::ImageDuplicator<TInputImage> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(InputImage);
  duplicator->Update();
  typename TInputImage::Pointer image = duplicator->GetModifiableOutput();

  float min;
  float max;
  float mean;
  ComputeRobustMinMaxMean<TInputImage, TMaskSpatialObject>(
    m_RemoveIntensityOutliers,
    image.GetPointer(),
    mask.GetPointer(), min, max, mean);
  itk::ImageRegionIterator<TInputImage> fi(image, image->GetBufferedRegion() );
  while( !fi.IsAtEnd() )
    {
    if( fi.Value() > max )
      {
      fi.Set(max);
      }
    if( fi.Value() < min )
      {
      fi.Set(min);
      }
    ++fi;
    }

  return image;
}

#endif // __BRAINSFITUTILS_h
