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
#ifndef __DemonsPreprocessor_hxx
#define __DemonsPreprocessor_hxx

#include "DemonsPreprocessor.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkBOBFFilter.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIO.h"
#include "itkMedianImageFilter.h"

#include "itkMultiResolutionPDEDeformableRegistration.h"
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkWarpImageFilter.h"
#include <sstream>
#include "DebugImageWrite.h"
namespace itk
{
template <typename TInputImage, typename TOutputImage>
DemonsPreprocessor<TInputImage, TOutputImage>::DemonsPreprocessor()
  : m_NumberOfHistogramLevels(256)
  , m_NumberOfMatchPoints(1)
  , m_FixedImageMinimum(NumericTraits<InputPixelType>::NonpositiveMin())
  , m_MovingImageMinimum(NumericTraits<InputPixelType>::NonpositiveMin())
  , m_FixedBinaryVolume("none")
  , m_MovingBinaryVolume("none")
  , m_Lower(NumericTraits<PixelType>::NonpositiveMin())
  , m_Upper(NumericTraits<PixelType>::max())
  , m_DefaultPixelValue(NumericTraits<PixelType>::OneValue())
  , m_OutDebug(false)
  , m_UseHistogramMatching(0)
{
  for (unsigned i = 0; i < TInputImage::ImageDimension; ++i)
  {
    m_Seed[i] = 0;
    m_MedianFilterSize[i] = 0;
  }
  m_Radius.Fill(1);
}

template <typename TInputImage, typename TOutputImage>
void
DemonsPreprocessor<TInputImage, TOutputImage>::Execute()
{
  if (m_MedianFilterSize[0] > 0 || m_MedianFilterSize[1] > 0 || m_MedianFilterSize[2] > 0)
  {
    using MedianImageFilterType = typename itk::MedianImageFilter<TInputImage, TInputImage>;
    typename MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
    medianFilter->SetRadius(m_MedianFilterSize);
    medianFilter->SetInput(m_InputFixedImage);
    medianFilter->Update();
    m_InputFixedImage = medianFilter->GetOutput();
    DebugOutput(InputImageType, m_InputFixedImage);

    //
    // reinitialize
    medianFilter = MedianImageFilterType::New();
    medianFilter->SetRadius(m_MedianFilterSize);
    medianFilter->SetInput(m_InputMovingImage);
    medianFilter->Update();
    m_InputMovingImage = medianFilter->GetOutput();
    DebugOutput(InputImageType, m_InputMovingImage);
  }
  { // Create UnNormalized...Images
    {
      this->m_UnNormalizedFixedImage = itkUtil::PreserveCast<TInputImage, TOutputImage>(this->m_InputFixedImage);
    }
    {
      m_UnNormalizedMovingImage = itkUtil::PreserveCast<TInputImage, TOutputImage>(this->m_InputMovingImage);
    }
  }

  //  m_OutputFixedImage =
  // itkUtil::CopyImage<TOutputImage>(m_UnNormalizedFixedImage);

  m_OutputMovingImage = itkUtil::CopyImage<TOutputImage>(m_UnNormalizedMovingImage);

  if (this->GetUseHistogramMatching())
  {
    using HistogramMatchingFilterType = HistogramMatchingImageFilter<OutputImageType, OutputImageType>;
    typename HistogramMatchingFilterType::Pointer histogramfilter = HistogramMatchingFilterType::New();
    if (this->GetOutDebug())
    {
      std::cout << "Performing Histogram Matching \n";
    }
    if ((std::numeric_limits<typename OutputImageType::PixelType>::max() -
         std::numeric_limits<typename OutputImageType::PixelType>::min()) < m_NumberOfHistogramLevels)
    {
      std::cout << "The intensity of range is less than Histogram levels!!" << std::endl;
    }
    histogramfilter->SetInput(m_UnNormalizedMovingImage);
    histogramfilter->SetReferenceImage(m_UnNormalizedFixedImage);

    histogramfilter->SetNumberOfHistogramLevels(m_NumberOfHistogramLevels);
    histogramfilter->SetNumberOfMatchPoints(m_NumberOfMatchPoints);
    histogramfilter->ThresholdAtMeanIntensityOn();
    histogramfilter->Update();
    m_OutputMovingImage = histogramfilter->GetOutput();
    DebugOutput(OutputImageType, m_OutputMovingImage);
    // +DANGER: ALIASING:  m_OutputMovingImage EQ m_UnNormalizedMovingImage by
    // design.
    // Create a copy just to be safe.  This s probably a waste of memory.
    // m_OutputMovingImage =
    // itkUtil::CopyImage<TOutputImage>(m_UnNormalizedMovingImage);
  }

  m_OutputFixedImage = itkUtil::CopyImage<TOutputImage>(m_UnNormalizedFixedImage);

  if (this->GetOutDebug())
  {
    std::cout << "Writing Histogram equalized image" << std::endl;
    itkUtil::WriteImage<TOutputImage>(m_OutputFixedImage, "HistogramModifiedFixedImage.nii.gz");
    std::cout << "Writing UnormalizedMovingImage equalized image" << std::endl;
    itkUtil::WriteImage<TOutputImage>(m_UnNormalizedMovingImage, "HistogramReferenceMovingImage.nii.gz");
  }
  // Make BOBF Images if specified
  if (this->m_FixedBinaryVolume != std::string("none"))
  {
    if (this->GetOutDebug())
    {
      std::cout << "Making BOBF \n";
      std::cout << "PRE Fixed Origin" << m_OutputFixedImage->GetOrigin() << std::endl;
    }
    m_OutputFixedImage = this->MakeBOBFImage(m_OutputFixedImage, m_FixedBinaryVolume);
    DebugOutput(OutputImageType, m_OutputFixedImage);
    if (this->GetOutDebug())
    {
      std::cout << "Fixed Origin" << m_OutputFixedImage->GetOrigin() << std::endl;
      std::cout << "PRE Moving Origin" << m_OutputMovingImage->GetOrigin() << std::endl;
    }
    m_OutputMovingImage = this->MakeBOBFImage(m_OutputMovingImage, m_MovingBinaryVolume);
    DebugOutput(OutputImageType, m_OutputMovingImage);
    if (this->GetOutDebug())
    {
      std::cout << "Moving Origin" << m_OutputMovingImage->GetOrigin() << std::endl;
      std::cout << "Writing Brain Only Background Filled Moving image" << std::endl;
      itkUtil::WriteImage<TOutputImage>(m_OutputMovingImage, "BOBF_Moving.nii.gz");
      itkUtil::WriteImage<TOutputImage>(m_OutputFixedImage, "BOBF_Fixed.nii.gz");
    }
  }
  m_InputMovingImage = nullptr;
  m_InputFixedImage = nullptr;
}

/*This function takes in a brain image and a whole brain mask and strips the
 * skull of the image. It uses the BOBF filter to perform the skull
 * stripping.*/

template <typename TInputImage, typename TOutputImage>
typename DemonsPreprocessor<TInputImage, TOutputImage>::OutputImagePointer
DemonsPreprocessor<TInputImage, TOutputImage>::MakeBOBFImage(OutputImagePointer input, std::string MaskName)
{
  OutputImagePointer Mask = itkUtil::ReadImage<OutputImageType>(MaskName);

  if ((m_UnNormalizedFixedImage->GetLargestPossibleRegion().GetSize() != Mask->GetLargestPossibleRegion().GetSize()) ||
      (m_UnNormalizedFixedImage->GetSpacing() != Mask->GetSpacing()))
  {
    if (this->GetOutDebug())
    {
      std::cout << "Writing Resampled Output image" << std::endl;
      itkUtil::WriteImage<TOutputImage>(Mask, "Resampled.mask");
    }
  }

  using BOBFFilterType = BOBFFilter<OutputImageType, OutputImageType>;
  typename BOBFFilterType::Pointer BOBFfilter = BOBFFilterType::New();
  if (this->GetOutDebug())
  {
    std::cout << "Making Brain only Background filled image with the following parameters. " << std::endl;
    std::cout << "Lower Threshold:  " << m_Lower << std::endl;
    std::cout << "Upper Threshold:  " << m_Upper << std::endl;
    std::cout << "Neighborhood:  " << m_Radius << std::endl;
    std::cout << "Background fill Value:  " << m_DefaultPixelValue << std::endl;
    std::cout << "Seed :  " << m_Seed << std::endl;
  }

  BOBFfilter->SetLower(m_Lower);
  BOBFfilter->SetUpper(m_Upper);
  BOBFfilter->SetRadius(m_Radius);
  BOBFfilter->SetReplaceValue(m_DefaultPixelValue);
  BOBFfilter->SetSeed(m_Seed);
  BOBFfilter->SetInputImage(input);
  BOBFfilter->SetInputMask(Mask);
  try
  {
    BOBFfilter->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << err << std::endl;
    throw;
  }

  OutputImagePointer output = BOBFfilter->GetOutput();
  return output;
}
} // namespace itk

#endif // _DemonsPreprocessor_hxx
