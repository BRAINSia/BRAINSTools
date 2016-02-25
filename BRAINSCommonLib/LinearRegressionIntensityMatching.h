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

// When averaging multiple images together, it is implicitly assumed that
// the images occupy the same dynamic range, and that each anatomical
// intensity region is represented by the same intensity values.  That is
// to say that it is assumed that a given intensity value has the same
// intepretation accross all of the averaged images.
//
// In practice, the above assumptions are not always true.  In order
// to improve the averaging interpretation, we will minimize the
// difference between the two data set through linear regression
// before averaging.
//
// The intent is that the the averaging filter will recieve
//
// AveFilter->SetInput1(ReferenceImage)
// AVeFilter->SetInput2(  LinearRegressionIntensityMatching(ReferenceImage, MaskImage, RescaleToReferenceDynamicRange) )
//
// For this to work, it is requiered that ReferenceImage, MaskImage, RescaleToReferenceDynamicRange all have the same
// physical space definitions and voxel space layout.
//
// The following code does a simple linear regression of the
// RescaleToReferenceDynamicRange image (x) to the original ReferenceImage (y)
// so that the intensity transfer function new_out = slope * x + intercept
// can be applied such that new_out is approximately the same dynamic range
// as the orginal input image.
#ifndef LinearRegressionIntensityMatching_h
#define LinearRegressionIntensityMatching_h

#include <itkImage.h>
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

template <typename TRealImage,
          typename TBinaryImage>
typename TRealImage::Pointer LinearRegressionIntensityMatching(
    const typename TRealImage::Pointer & ReferenceImage,
    typename TBinaryImage::Pointer & MaskImage,
    const typename TRealImage::Pointer & RescaleToReferenceDynamicRange)
{
  const bool doEntireImage = MaskImage.IsNull();  // If the mask image is null, assume entire image
  //Get average of values in Rescale Image

  if( ReferenceImage->GetOrigin() != RescaleToReferenceDynamicRange->GetOrigin() ||
      ReferenceImage->GetOrigin() != MaskImage->GetOrigin() )
    {
    std::cout<< "Image data spacing mismatch ::\n"
             << " - ReferenceImage:: " << ReferenceImage->GetOrigin() << std::endl
             << " - RescaleToReferenceDynamicRange:: " << RescaleToReferenceDynamicRange->GetOrigin() << std::endl
             << " - MaskImage:: " << MaskImage->GetOrigin() << std::endl;
    }
  itk::ImageRegionIterator<TRealImage> ItRescaledImage( RescaleToReferenceDynamicRange,
    RescaleToReferenceDynamicRange->GetRequestedRegion() );
  typedef typename TRealImage::PixelType RescaleRealType;
  RescaleRealType avgOut = 0.0;
  size_t nOut = 0;
  for( ItRescaledImage.GoToBegin(); !ItRescaledImage.IsAtEnd(); ++ItRescaledImage )
    {
    if( doEntireImage ||  ( MaskImage->GetPixel(ItRescaledImage.GetIndex()) != 0 ) )
      {
      avgOut += ItRescaledImage.Get();
      ++nOut;
      }
    }
  avgOut /= static_cast<RescaleRealType>(nOut);

  //Get average of values in Reference image
  itk::ImageRegionConstIterator<TRealImage> ItRefImg( ReferenceImage,
    ReferenceImage->GetRequestedRegion() );
  RescaleRealType avgIn = 0.0;
  size_t nIn = 0;
  for( ItRefImg.GoToBegin(); !ItRefImg.IsAtEnd(); ++ItRefImg )
    {
    // If in mask
      {
      avgIn += ItRefImg.Get();
      ++nIn;
      }
    }
  avgIn /= static_cast<RescaleRealType>(nIn);

  assert(nIn == nOut);

  RescaleRealType numerator = 0.0;
  RescaleRealType denominator = 0.0;
  for( ItRefImg.GoToBegin(), ItRescaledImage.GoToBegin(); !ItRefImg.IsAtEnd() &&
       !ItRescaledImage.IsAtEnd(); ++ItRefImg,++ItRescaledImage )
    {
    // If in mask
      {
      const RescaleRealType out_sub_avg = ItRescaledImage.Get() - avgOut;
      const RescaleRealType in_sub_avg  = ItRefImg.Get() -avgIn;
      numerator += out_sub_avg*in_sub_avg;
      denominator += out_sub_avg*out_sub_avg;
      }
    }
  if( denominator < 1e-10 )
    {
    denominator = 1.0;
    }
  const RescaleRealType slope = numerator/denominator;
  const RescaleRealType intercept = avgIn - slope*avgOut;
  //std::cout << "Slope = " << slope << " Intercept= " << intercept << std::endl;
  for(ItRescaledImage.GoToBegin(); !ItRescaledImage.IsAtEnd(); ++ItRescaledImage)
    {
    ItRescaledImage.Set( ItRescaledImage.Get()*slope + intercept);
    }

  return RescaleToReferenceDynamicRange;

}

#endif // LinearRegressionIntensityMatching_h
