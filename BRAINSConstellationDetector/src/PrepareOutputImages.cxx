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

#include "itkIO.h"
#include <BRAINSFitHelper.h>

#include "PrepareOutputImages.h"
#include "itkOrthogonalize3DRotationMatrix.h"
#include "ChopImageBelowLowerBound.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"

namespace itk
{

  SImageType::PointType GetNamedPointFromLandmarkList(const LandmarksMapType & landmarks, const std::string & NamedPoint)
    {
    LandmarksMapType::const_iterator itpair = landmarks.find(NamedPoint);

    if( itpair == landmarks.end() )
      {
      std::cout << "ERROR:  " << NamedPoint << " not found in list." << std::endl;
      return SImageType::PointType();
      }
    return itpair->second;
    }

  void PrepareOutputImages(SImageType::Pointer & lOutputResampledImage,
    SImageType::Pointer & lOutputImage,
    SImageType::Pointer & lOutputUntransformedClippedVolume,
    SImageType::ConstPointer lImageToBeResampled,
    VersorTransformType::ConstPointer lVersorTransform,
    const double lACLowerBound,
    const short int BackgroundFillValue,
    const std::string & lInterpolationMode,
    const bool lCutOutHeadInOutputVolume,
    const double lOtsuPercentileThreshold
  )
    {
    typedef ResampleInPlaceImageFilter<SImageType, SImageType> ResampleIPFilterType;
    typedef ResampleIPFilterType::Pointer                      ResampleIPFilterPointer;

    const double PhysicalLowerBound = /* ACy when zero-centered is ... */ 0.0 - lACLowerBound;
      {
      const SImageType * constImage( lImageToBeResampled.GetPointer() );

      ResampleIPFilterPointer resampleIPFilter = ResampleIPFilterType::New();
      resampleIPFilter->SetInputImage( constImage );
      resampleIPFilter->SetRigidTransform( lVersorTransform.GetPointer() );
      resampleIPFilter->Update();
      lOutputImage = resampleIPFilter->GetOutput();
      }

      {
      lOutputResampledImage = TransformResample<SImageType, SImageType>(
        lImageToBeResampled.GetPointer(),
        MakeIsoTropicReferenceImage().GetPointer(),
        BackgroundFillValue,
        GetInterpolatorFromString<SImageType>( lInterpolationMode).GetPointer(),
        lVersorTransform.GetPointer() );
      }

      {
      // ======================== Start
      // ======================== Start
      //  HACK -- chopping based on AcLowerBound
      //  This is ugly code that could be re-written much simpler.
      //
      typedef itk::ImageRegionIteratorWithIndex<SImageType> IteratorType;
      const double thousand = 1000.0;   // we need a DOUBLE constant, not a
      // FLOAT constant, for exact switch
      // comparisons.
      if( lACLowerBound < thousand )
        {
        // First Process the OutputResampledImage
        std::cout << "Chopping image below physical location: " << PhysicalLowerBound << "." << std::endl;
        ChopImageBelowLowerBound<SImageType>(lOutputResampledImage, BackgroundFillValue, PhysicalLowerBound);
        ChopImageBelowLowerBound<SImageType>(lOutputImage, BackgroundFillValue, PhysicalLowerBound);

        // Second Create a mask for inverse resampling to orignal space
        SImageType::Pointer ZeroOneImage = SImageType::New();
        ZeroOneImage->CopyInformation(lOutputResampledImage);
        ZeroOneImage->SetRegions( lOutputResampledImage->GetLargestPossibleRegion() );
        ZeroOneImage->Allocate();
        ZeroOneImage->FillBuffer(1);
        ChopImageBelowLowerBound<SImageType>(ZeroOneImage, BackgroundFillValue, PhysicalLowerBound);

        if( lCutOutHeadInOutputVolume )  // Restrict mask to head
          // tissue region if necessary
          {
          //  No double opportunity when generating both kinds of images.
          const unsigned int closingSize = 7;
          typedef itk::LargestForegroundFilledMaskImageFilter<SImageType> LFFMaskFilterType;
          LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
          LFF->SetInput(lOutputResampledImage);
          LFF->SetOtsuPercentileThreshold(lOtsuPercentileThreshold);
          LFF->SetClosingSize(closingSize);
          LFF->Update();
          SImageType::Pointer HeadOutlineMaskImage =  LFF->GetOutput();

          IteratorType ItZeroOneImage( ZeroOneImage, ZeroOneImage->GetRequestedRegion() );
          ItZeroOneImage.GoToBegin();
          IteratorType ItOutputResampledImage( lOutputResampledImage,
            lOutputResampledImage->GetRequestedRegion() );
          ItOutputResampledImage.GoToBegin();
          IteratorType ItHead( HeadOutlineMaskImage, HeadOutlineMaskImage->GetLargestPossibleRegion() );
          ItHead.GoToBegin();
          while( !ItHead.IsAtEnd() )
            {
            if( ItHead.Get() == 0 )
              {
              ItOutputResampledImage.Set(0);
              ItZeroOneImage.Set(0);
              }
            ++ItZeroOneImage;
            ++ItOutputResampledImage;
            ++ItHead;
            }
          }
        // Map the ZeroOne image through the inverse zero-centered transform
        // to make the clipping factor image:
        typedef itk::NearestNeighborInterpolateImageFunction<SImageType, double> NearestNeighborInterpolatorType;
        NearestNeighborInterpolatorType::Pointer interpolator = NearestNeighborInterpolatorType::New();
        typedef itk::ResampleImageFilter<SImageType, SImageType> ResampleFilterType;
        ResampleFilterType::Pointer ResampleFilter = ResampleFilterType::New();
        ResampleFilter->SetInput(ZeroOneImage);
        ResampleFilter->SetInterpolator(interpolator);
        ResampleFilter->SetDefaultPixelValue(0);
        ResampleFilter->SetOutputParametersFromImage(lImageToBeResampled);
          {
          VersorTransformType::Pointer lInvVersorTransform = VersorTransformType::New();
          const SImageType::PointType centerPoint = lVersorTransform->GetCenter();
          lInvVersorTransform->SetCenter(centerPoint);
          lInvVersorTransform->SetIdentity();
          lVersorTransform->GetInverse(lInvVersorTransform);
          ResampleFilter->SetTransform(lInvVersorTransform);
          }
        ResampleFilter->Update();
        SImageType::Pointer lClippingFactorImage = ResampleFilter->GetOutput();

        // Multiply the raw input image by the clipping factor image:
        typedef itk::MultiplyImageFilter<SImageType, SImageType> MultiplyFilterType;
        MultiplyFilterType::Pointer MultiplyFilter = MultiplyFilterType::New();
        MultiplyFilter->SetInput1(lImageToBeResampled);
        MultiplyFilter->SetInput2(lClippingFactorImage);
        MultiplyFilter->Update();
        lOutputUntransformedClippedVolume = MultiplyFilter->GetOutput();
        }
      // ======================== Stop
      // ======================== Stop
      }
    }

  void ApplyInverseOfTransformToLandmarks(
    VersorTransformType::ConstPointer lVersorTransform,
    const LandmarksMapType & inputLmks,
    LandmarksMapType & outputLmks
  )
    {
    VersorTransformType::Pointer lInvVersorTransform = VersorTransformType::New();
      {
      const SImageType::PointType centerPoint = lVersorTransform->GetCenter();
      lInvVersorTransform->SetCenter(centerPoint);
      lInvVersorTransform->SetIdentity();
      lVersorTransform->GetInverse(lInvVersorTransform);
      }

    outputLmks.clear();
    for( LandmarksMapType::const_iterator lit = inputLmks.begin();
      lit != inputLmks.end(); ++lit )
      {
      outputLmks[lit->first] =
        lInvVersorTransform->TransformPoint
        ( GetNamedPointFromLandmarkList(inputLmks, lit->first) );
      }
    }
}
