#include "PrepareOutputImages.h"
#include "ChopImageBelowLowerBound.h"

namespace itk
{
  void PrepareOutputImages(SImageType::Pointer & lOutputResampledImage,
    SImageType::Pointer & lOutputImage,
    SImageType::Pointer & lOutputUntransformedClippedVolume,
    SImageType::Pointer & lImageToBeResampled,
    VersorTransformType::Pointer & lInvVersorTransform,
    const double lACLowerBound,
    const double BackgroundFillValue,
    const double PhysicalLowerBound,
    const bool lCutOutHeadInOutputVolume,
    const double lOtsuPercentileThreshold
  )
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
      // don't forget to do SetRegions here!
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

      ResampleFilter->SetTransform(lInvVersorTransform);
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
