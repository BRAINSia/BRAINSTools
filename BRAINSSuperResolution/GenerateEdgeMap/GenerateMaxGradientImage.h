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
 * Author: Ali Ghayoor
 * at SINAPSE Lab,
 * The University of Iowa 2016
 */

#ifndef __GenerateMaxGradientImage_h
#define __GenerateMaxGradientImage_h

#include "itkGradientMagnitudeImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkTimeProbe.h"
#include <vector>

// Class to print the proper exception message
class EmptyVectorException
{
public:
  EmptyVectorException(const char* pStr = "The list of input images was empty.  Nothing to find maximum!") :
    pMessage(pStr)
    {
    }

  const char * what() const
  {
  return pMessage;
  }

private:
  const char * pMessage;
};

// Auxiliary function to find the maximum of the rescaled gradient image list
template <typename TImage>
typename TImage::Pointer
MaxOfImageList(const std::vector<typename TImage::Pointer> & inputImageList) // inputImageList = rescaledGradientImageList
{
  typedef itk::MaximumImageFilter<TImage, TImage, TImage>   MaximumFilterType;

  if( inputImageList.empty() )
    {
    // No images, something went wrong.
    throw EmptyVectorException();
    }
  if( inputImageList.size() == 1 )
    {
    // Only one image, no need to find maximum image.
    return inputImageList[0];
    }

  // Initialize the maximum image with the first image in the list
  typename TImage::Pointer maxImage = inputImageList[0];

  for(unsigned int i = 1; i < inputImageList.size(); ++i)
     {
     typename MaximumFilterType::Pointer myMax = MaximumFilterType::New();
     myMax->SetInput1( maxImage );
     myMax->SetInput2( inputImageList[i] );
     try
       {
       myMax->Update();
       }
     catch( itk::ExceptionObject & exp )
       {
       std::cerr << "ExceptionObject with Iterator" << std::endl;
       std::cerr << exp << std::endl;
       }
     maxImage = myMax->GetOutput();
     }

  return maxImage;
}

/*
 * Main function to generate the maximum gradient image.
 * Weak edges (whose gradients are less than %50 quantile) are set to minOutputRange (1).
 * Strong edges (whose gradients are more than %95 quantile) are set to maxOutputRange(255).
 *
 *** NOTE ***
 * All input images and the input mask must be in the same voxel space
 ************
 */
template <class InputImageType, class OutputImageType, class MaskImageType>
typename OutputImageType::Pointer
GenerateMaxGradientImage(const std::vector<typename InputImageType::Pointer> & inputImages,
                         const float LowerPercentileMatching, // Map lower quantile and below to minOutputRange
                         const float UpperPercentileMatching, // Map upper quantile and above to maxOutputRange
                         const unsigned int minOutputRange,   // epsilon
                         const unsigned int maxOutputRange,
                         typename MaskImageType::Pointer mask = NULL)
{
  std::cout << "Generating maximum gradient image..." << std::endl;
  std::cout << "[LowerQuantile UpperQuantile] = [" << LowerPercentileMatching << " " << UpperPercentileMatching << "]" << std::endl;
  std::cout << "[minOutputRange maxOutputRange] [" << minOutputRange << " " << maxOutputRange << "]" << std::endl;

  typedef itk::GradientMagnitudeImageFilter<InputImageType, InputImageType>   GradientFilterType;
  typedef itk::MinimumMaximumImageCalculator<InputImageType>                  MinMaxCalculatorType;
  typedef itk::LabelStatisticsImageFilter<InputImageType, MaskImageType>      LabelStatisticsImageFilter;
  typedef typename itk::IntensityWindowingImageFilter<InputImageType,
                                                      OutputImageType>        WindowRescalerType;
  typedef std::vector<typename OutputImageType::Pointer>                      RescaledImageGradientVectorType;

  itk::TimeProbe MaxGradientImageTimer;
  MaxGradientImageTimer.Start();

  const unsigned int numberOfImageModalities =
    inputImages.size(); // number of modality images

  // list of rescaled gradient magnitude images
  RescaledImageGradientVectorType rescaledGradientImageList( numberOfImageModalities );

  const typename MaskImageType::PixelType maskInteriorLabel = 1;
  typename MaskImageType::Pointer internalMask;
  if( mask.IsNull() )
    {
    internalMask = MaskImageType::New();
    internalMask->CopyInformation( inputImages[0] );
    internalMask->SetRegions( inputImages[0]->GetLargestPossibleRegion() );
    internalMask->Allocate();
    internalMask->FillBuffer( maskInteriorLabel );
    }
  else
    {
    internalMask = mask;
    }

  for( size_t i = 0; i < numberOfImageModalities; i++ )
     {
     typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
     gradientFilter->SetInput( inputImages[i] );
     gradientFilter->Update();

     typename MinMaxCalculatorType::Pointer myMinMax = MinMaxCalculatorType::New();
     myMinMax->SetImage( gradientFilter->GetOutput() );
     myMinMax->Compute();
     typename InputImageType::PixelType imgMin = myMinMax->GetMinimum();
     typename InputImageType::PixelType imgMax = myMinMax->GetMaximum();

     int numBins = vnl_math_rnd(imgMax - imgMin + 1);
     if( numBins < 256 )
       {
       numBins = 256;
       }

     typename LabelStatisticsImageFilter::Pointer maskedStatistics = LabelStatisticsImageFilter::New();
     maskedStatistics->SetInput( gradientFilter->GetOutput() );
     maskedStatistics->SetLabelInput( internalMask );
     maskedStatistics->UseHistogramsOn();
     maskedStatistics->SetHistogramParameters(numBins, imgMin, imgMax);
     maskedStatistics->Update();
     typename LabelStatisticsImageFilter::HistogramType::Pointer hist =
       maskedStatistics->GetHistogram( maskInteriorLabel );
     if( hist.IsNull() )
       {
       itkGenericExceptionMacro("histogram had no value for label "
                                << maskInteriorLabel);
       }

     /*
     // maximum is defined such that (Max - Epsilon)/2 = Q(%85)
     unsigned int maxOutputRange = 2 * hist->Quantile(0, 0.85F) + minOutputRange;
     if( maxOutputRange > 255 )
       {
       maxOutputRange = 255;
       }
     */

     typename WindowRescalerType::Pointer intensityMapper = WindowRescalerType::New();
     intensityMapper->SetInput( maskedStatistics->GetOutput() );
     intensityMapper->SetOutputMinimum( minOutputRange );
     intensityMapper->SetOutputMaximum( maxOutputRange );
     intensityMapper->SetWindowMinimum( hist->Quantile(0, LowerPercentileMatching) );
     intensityMapper->SetWindowMaximum( hist->Quantile(0, UpperPercentileMatching) );
     intensityMapper->Update();

     rescaledGradientImageList[i] = intensityMapper->GetOutput();
     }

  typename OutputImageType::Pointer maxGradientImage = MaxOfImageList<OutputImageType>(rescaledGradientImageList);
/*
  // Another rescaling was needed if we were computing summed gradient image.
  typedef itk::IntensityWindowingImageFilter<OutputImageType,
                                             OutputImageType>           RescaleFilterType;
  typename RescaleFilterType::Pointer outputRescaler = RescaleFilterType::New();
  outputRescaler->SetOutputMinimum(0);
  outputRescaler->SetOutputMaximum(255);
  outputRescaler->SetInput( maxGradientImage );
  outputRescaler->Update();
*/
  MaxGradientImageTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime = MaxGradientImageTimer.GetTotal();
  std::cout << "Generating maximum gradient image took " << elapsedTime
            << " " << MaxGradientImageTimer.GetUnit() << "." << std::endl;

  //return outputRescaler->GetOutput();
  return maxGradientImage;
}

#endif // __GenerateMaxGradientImage_h
