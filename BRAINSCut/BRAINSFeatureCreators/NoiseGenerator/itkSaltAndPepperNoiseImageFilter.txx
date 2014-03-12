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
#include "itkSaltAndPepperNoiseImageFilter.h"
#include "itkThreadSafeMersenneTwisterRandomVariateGenerator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
SaltAndPepperNoiseImageFilter<TInputImage, TOutputImage>
::SaltAndPepperNoiseImageFilter()
{
  m_Probability = 0.01;
}

template <class TInputImage, class TOutputImage>
void
SaltAndPepperNoiseImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
                        ThreadIdType threadId)
{
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput(0);

  // create a random generator per thread
  typename Statistics::ThreadSafeMersenneTwisterRandomVariateGenerator::Pointer rand =
    Statistics::ThreadSafeMersenneTwisterRandomVariateGenerator::New();
  rand->Initialize();

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
    if( rand->GetVariate() < m_Probability )
      {
      if( rand->GetVariate() < 0.5 )
        {
        // salt
        outputIt.Set( NumericTraits<OutputImagePixelType>::max() );
        }
      else
        {
        // pepper
        outputIt.Set( NumericTraits<OutputImagePixelType>::NonpositiveMin() );
        }
      }
    else
      {
      // keep the data unchanged
      outputIt.Set( (OutputImagePixelType) inputIt.Get() );
      }
    ++inputIt;
    ++outputIt;
    progress.CompletedPixel();  // potential exception thrown here
    }
}

template <class TInputImage, class TOutputImage>
void
SaltAndPepperNoiseImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os,
            Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Probability: "
     << static_cast<typename NumericTraits<double>::PrintType>(this->GetProbability() )
     << std::endl;
}
} /* namespace itk */
