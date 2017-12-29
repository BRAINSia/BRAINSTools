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
/*=========================================================================
 Author: $Author: krm15 $  // Author of last commit
 Version: $Rev: 585 $  // Revision of last commit
 Date: $Date: 2009-08-20 21:25:19 -0400 (Thu, 20 Aug 2009) $  // Date of last commit
 =========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 =========================================================================*/

#ifndef __itkHoughTransformRadialVotingImageFilter_hxx
#define __itkHoughTransformRadialVotingImageFilter_hxx

#include "itkHoughTransformRadialVotingImageFilter.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::HoughTransformRadialVotingImageFilter() :
  m_MinimumRadius(0),   // by default
  m_MaximumRadius(10),  // by default
  m_Threshold(0),
  m_GradientThreshold(0),
  m_OutputThreshold(0.0),
  m_SigmaGradient(1),   // Scale of the DoG filter
  m_Variance(1),
  m_VotingRadiusRatio(0.5),
  m_SphereRadiusRatio(1),
  m_SamplingRatio(1.0),
  m_RadiusImage(nullptr),
  m_AccumulatorImage(nullptr),
  m_SpheresList(),
  m_NumberOfSpheres(1),
  m_OldModifiedTime(0),
  m_NbOfThreads(1),
  m_AllSeedsProcessed(false),
  m_HoughEyeDetectorMode(0)
{
  // Limiting the Hough Filter to run with single thread.
  // In a single threaded status, the increment of run time for this filter is not tremendous in compare with
  // the total run time of the BCD
  this->SetNumberOfThreads(1);
}

template <class TInputImage, class TOutputImage>
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::~HoughTransformRadialVotingImageFilter()
{
}

template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::SetRadius(InputCoordType radius)
{
  this->SetMinimumRadius(radius);
  this->SetMaximumRadius(radius);
}

template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);

  output->SetRequestedRegionToLargestPossibleRegion();
}

template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  if( this->GetInput() )
    {
    InputImagePointer image =
      const_cast<InputImageType *>( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}

template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  // Get the input and output pointers
  const InputImageConstPointer inputImage = this->GetInput();

  m_AccumulatorImage = InternalImageType::New();
  m_AccumulatorImage->CopyInformation(inputImage);
  m_AccumulatorImage->SetRegions( inputImage->GetLargestPossibleRegion() );

  m_AccumulatorImage->Allocate();
  m_AccumulatorImage->FillBuffer(0);

  m_RadiusImage = InternalImageType::New();
  m_RadiusImage->CopyInformation( inputImage );
  m_RadiusImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  m_RadiusImage->Allocate();
  m_RadiusImage->FillBuffer(0);
}

template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  ComputeMeanRadiusImage();

  // Copy the typecast m_AccumulatorImage to Output image
  InternalIteratorType iIt( m_AccumulatorImage, this->GetOutput()->GetRequestedRegion() );
  OutputIteratorType   oIt( this->GetOutput(), this->GetOutput()->GetRequestedRegion() );

  iIt.GoToBegin();
  oIt.GoToBegin();
  while( !iIt.IsAtEnd() )
    {
    assert( !oIt.IsAtEnd() );
    oIt.Set( static_cast<OutputPixelType>( iIt.Get() ) );
    ++iIt;
    ++oIt;
    }
}

template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & windowRegion,
  ThreadIdType /* threadId -- Not used */)

{
  // Get the input and output pointers
  const InputImageConstPointer inputImage = this->GetInput();
  const InputSpacingType       spacing = inputImage->GetSpacing();

  DoGFunctionPointer DoGFunction = DoGFunctionType::New();
  DoGFunction->SetUseImageSpacing(true);

  DoGFunction->SetInputImage(inputImage);
  DoGFunction->SetSigma(m_SigmaGradient);

  GaussianFunctionPointer GaussianFunction = GaussianFunctionType::New();

  const InputCoordType averageRadius = 0.5 * ( m_MinimumRadius + m_MaximumRadius );
  const InputCoordType averageRadius2 = averageRadius * averageRadius;

  ImageRegionConstIteratorWithIndex<InputImageType>
  image_it(inputImage, windowRegion);
  image_it.GoToBegin();

  const unsigned int sampling = static_cast<unsigned int>( 1. / m_SamplingRatio );
  unsigned int       counter = 1;

  InternalRegionType region;
  while( !image_it.IsAtEnd() )
    {
    if( image_it.Get() > m_Threshold )
      {
      const Index<ImageDimension> index = image_it.GetIndex();
      DoGVectorType               grad = DoGFunction->EvaluateAtIndex(index);

      // if the gradient is not flat
      typename DoGVectorType::ValueType norm2 = grad.GetSquaredNorm();

      if( norm2 > m_GradientThreshold )
        {
        if( counter % sampling == 0 )
          {
          // Normalization
          if( norm2 != 0 )
            {
            const typename DoGVectorType::ValueType inv_norm = 1.0 / std::sqrt(norm2);
            for( unsigned int i = 0; i < ImageDimension; i++ )
              {
              grad[i] *= inv_norm;
              }
            }
          Index<ImageDimension> center;
            {
            InternalIndexType start;
            InternalSizeType  size;
            for( unsigned int i = 0; i < ImageDimension; i++ )
              {
              // for T1, T2 images
              if( m_HoughEyeDetectorMode == 1 )
                {
                center[i] = index[i] + static_cast
                  <InternalIndexValueType>( averageRadius * grad[i] / spacing[i] );
                }
              else
                { // for PD image
                center[i] = index[i] - static_cast
                  <InternalIndexValueType>( averageRadius * grad[i] / spacing[i] );
                }

              const InputCoordType rad = m_VotingRadiusRatio * m_MinimumRadius / spacing[i];
              start[i] = center[i] - static_cast<InternalIndexValueType>( rad );
              size[i] = 1 + 2 * static_cast<InternalSizeValueType>( rad );
              }
            region.SetSize(size);
            region.SetIndex(start);
            }

          if( inputImage->GetRequestedRegion().IsInside(region) )
            {
            ImageRegionIteratorWithIndex<InternalImageType>
            It1(m_AccumulatorImage, region);

            ImageRegionIterator<InternalImageType> It2(m_RadiusImage, region);

            It1.GoToBegin();
            It2.GoToBegin();
            while( !It1.IsAtEnd() )
              {
              assert( !It2.IsAtEnd() );
              const Index<ImageDimension> indexAtVote = It1.GetIndex();
              InputCoordType              distance = 0;
              double                      d = 0;
              for( unsigned int i = 0; i < ImageDimension; i++ )
                {
                d += vnl_math_sqr(
                    static_cast<double>( indexAtVote[i] - center[i] ) * spacing[i]);
                distance += vnl_math_sqr(
                    static_cast<InputCoordType>( indexAtVote[i] - index[i] ) * spacing[i]);
                }
              d = std::sqrt(d);
              distance = std::sqrt(distance);

              // Apply a normal distribution weight;
              const double weight = GaussianFunction->EvaluatePDF(d, 0, averageRadius2);
              double       distweight(distance * weight);

              It1.Set(It1.Get() + weight);
              It2.Set(It2.Get() + distweight);

              ++It1;
              ++It2;
              }
            }
          } // end counter
        counter++;
        } // end gradient threshold
      }   // end intensity threshold
    ++image_it;
    }
}

template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>::ComputeMeanRadiusImage()
{
  const InputImageConstPointer inputImage = this->GetInput(0);
  const InputRegionType        windowRegion = inputImage->GetRequestedRegion();

  ImageRegionConstIterator<InternalImageType> acc_it(m_AccumulatorImage, windowRegion);
  ImageRegionIterator<InternalImageType>      radius_it(m_RadiusImage, windowRegion);

  acc_it.GoToBegin();
  radius_it.GoToBegin();
  while( !acc_it.IsAtEnd() )
    {
    assert(!radius_it.IsAtEnd() );
    if( acc_it.Get() > 0 )
      {
      radius_it.Set( radius_it.Get() / acc_it.Get() );
      }
    ++acc_it;
    ++radius_it;
    }
}

/** Get the list of circles. This recomputes the circles */
template <class TInputImage, class TOutputImage>
typename HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>::SpheresListType
& HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::GetSpheres()
  {
  // if the filter has not been updated
  if( this->GetMTime() == m_OldModifiedTime )
    {
    return m_SpheresList;
    }

  m_SpheresList.clear();

  /** Blur the accumulator in order to find the maximum */
  GaussianFilterPointer gaussianFilter = GaussianFilterType::New();
  gaussianFilter->SetInput( this->GetOutput(0) );   // the output is the
                                                    // accumulator image
  gaussianFilter->SetVariance(m_Variance);
  gaussianFilter->Update();

  const InternalImagePointer postProcessImage = gaussianFilter->GetOutput();
  const InternalSpacingType  spacing = postProcessImage->GetSpacing();
  const InternalSizeType     size = postProcessImage->GetRequestedRegion().GetSize();

  MinMaxCalculatorPointer minMaxCalculator = MinMaxCalculatorType::New();
  minMaxCalculator->SetImage(postProcessImage);
  minMaxCalculator->ComputeMaximum();
  const InternalPixelType pmax = minMaxCalculator->GetMaximum();

  unsigned int circles = 0;

  // Find maxima
  do
    {
    minMaxCalculator->ComputeMaximum();
    const InternalPixelType max = minMaxCalculator->GetMaximum();
    const InternalIndexType idx = minMaxCalculator->GetIndexOfMaximum();

    if( max < m_OutputThreshold * pmax )
      {
      assert( 0 == 1 ); // This should never happen.  max == pmax in all cases
      break;
      }

    SphereVectorType center;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      center[i] = idx[i];
      }

    // Create a Line Spatial Object
    SpherePointer Sphere = SphereType::New();
    Sphere->SetId(circles);
    Sphere->SetRadius( m_RadiusImage->GetPixel(idx) );
    Sphere->GetObjectToParentTransform()->SetOffset(center);
    Sphere->ComputeBoundingBox();

    m_SpheresList.push_back(Sphere);
    // Remove a black disc from the hough space domain
    InternalSizeType  sizeOfROI;
    InternalIndexType start;
    InternalIndexType end;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      const InternalIndexValueType rad =
        static_cast<InternalIndexValueType>(
          m_SphereRadiusRatio * Sphere->GetRadius()[i] / spacing[i] );

      if( idx[i] > static_cast<InternalIndexValueType>( rad ) )
        {
        start[i] = idx[i] - static_cast<InternalIndexValueType>( rad );
        }
      else
        {
        start[i] = 0;
        }

      if( idx[i] + rad < static_cast<InternalIndexValueType>( size[i] - 1 ) )
        {
        end[i] = idx[i] + static_cast<InternalIndexValueType>( rad );
        }
      else
        {
        end[i] = postProcessImage->GetLargestPossibleRegion().GetSize()[i] - 1;
        }

      sizeOfROI[i] = end[i] - start[i] + 1;
      }

    InternalRegionType region;
    region.SetSize(sizeOfROI);
    region.SetIndex(start);

    ImageRegionIterator<InternalImageType> It(postProcessImage, region);

    It.GoToBegin();
    while( !It.IsAtEnd() )
      {
      It.Set(0);
      ++It;
      }

    ++circles;
    }
  while( circles < m_NumberOfSpheres );

  m_OldModifiedTime = this->GetMTime();

  return m_SpheresList;
  }

/** Print Self information */
template <class TInputImage, class TOutputImage>
void
HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Threshold: " << m_Threshold << std::endl;
  os << "Gradient Threshold: " << m_GradientThreshold << std::endl;
  os << "Minimum Radius:  " << m_MinimumRadius << std::endl;
  os << "Maximum Radius: " << m_MaximumRadius << std::endl;
  os << "Derivative Scale : " << m_SigmaGradient << std::endl;
  os << "Sphere Radius Ratio: " << m_SphereRadiusRatio << std::endl;
  os << "Voting Radius Ratio: " << m_VotingRadiusRatio << std::endl;
  os << "Accumulator blur variance: " << m_Variance << std::endl;
  os << "Number Of Spheres: " << m_NumberOfSpheres << std::endl;
  os << "Output Threshold : " << m_OutputThreshold << std::endl;
  os << "Sampling Ratio: " << m_SamplingRatio << std::endl;
  os << "NbOfThreads: " << m_NbOfThreads << std::endl;
  os << "All Seeds Processed: " << m_AllSeedsProcessed << std::endl;
  os << "HoughEyeDetectorMode: " << m_HoughEyeDetectorMode << std::endl;

  os << "Radius Image Information : " << m_RadiusImage << std::endl;
}
} // end namespace

#endif
