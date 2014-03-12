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
#include "itkShotNoiseImageFilter.h"
#include "itkThreadSafeMersenneTwisterRandomVariateGenerator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"
#include "itkNormalVariateGenerator.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
ShotNoiseImageFilter<TInputImage, TOutputImage>
::ShotNoiseImageFilter()
{
  m_Scale = 1.0;
}

template <class TInputImage, class TOutputImage>
void
ShotNoiseImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
                        ThreadIdType threadId)
{
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput(0);

  // create a random generator per thread
  typename Statistics::ThreadSafeMersenneTwisterRandomVariateGenerator::Pointer rand =
    Statistics::ThreadSafeMersenneTwisterRandomVariateGenerator::New();
  rand->Initialize();
  typename Statistics::NormalVariateGenerator::Pointer randn = Statistics::NormalVariateGenerator::New();

  // Define the portion of the input to walk for this thread, using
  // the CallCopyOutputRegionToInputRegion method allows for the input
  // and output images to be different dimensions
  InputImageRegionType inputRegionForThread;
  this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);

  // Define the iterators
  ImageRegionConstIterator<TInputImage> inputIt(inputPtr, inputRegionForThread);
  ImageRegionIterator<TOutputImage>     outputIt(outputPtr, outputRegionForThread);

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels() );

  inputIt.GoToBegin();
  outputIt.GoToBegin();

  while( !inputIt.IsAtEnd() )
    {
    double in = m_Scale * inputIt.Get();
    if( in < 50 )
      {
      double L = vcl_exp( -in );
      long   k = 0;
      double p = 1.0;

      do
        {
        k += 1;
        p *= rand->GetVariate();
        }
      while( p > L );

      // clip the output to the actual supported range
      outputIt.Set( (OutputImagePixelType) std::min( (double)NumericTraits<OutputImagePixelType>::max(),
                                                     (k - 1) / m_Scale ) );
      }
    else
      {
      double out = in + vcl_sqrt( in ) * randn->GetVariate();
      outputIt.Set( (OutputImagePixelType) std::min( (double)NumericTraits<OutputImagePixelType>::max(), out
                                                     / m_Scale ) );
      }
    ++inputIt;
    ++outputIt;
    progress.CompletedPixel();  // potential exception thrown here
    }
}

template <class TInputImage, class TOutputImage>
void
ShotNoiseImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os,
            Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Scale: "
     << static_cast<typename NumericTraits<double>::PrintType>(this->GetScale() )
     << std::endl;
}
} /* namespace itk */
