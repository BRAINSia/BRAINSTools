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
#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"
#include "itkNormalVariateGenerator.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
AdditiveGaussianNoiseImageFilter<TInputImage, TOutputImage>
::AdditiveGaussianNoiseImageFilter()
{
  m_Mean = 0.0;
  m_StandardDeviation = 1.0;
}

template <class TInputImage, class TOutputImage>
void
AdditiveGaussianNoiseImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
                        ThreadIdType threadId)
{
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput(0);

  // create a random generator per thread
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
    double out = inputIt.Get() + m_Mean + m_StandardDeviation * randn->GetVariate();
    out = std::min( (double)NumericTraits<OutputImagePixelType>::max(), out );
    out = std::max( (double)NumericTraits<OutputImagePixelType>::NonpositiveMin(), out );
    outputIt.Set( (OutputImagePixelType) out  );
    ++inputIt;
    ++outputIt;
    progress.CompletedPixel();  // potential exception thrown here
    }
}

template <class TInputImage, class TOutputImage>
void
AdditiveGaussianNoiseImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os,
            Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Mean: "
     << static_cast<typename NumericTraits<double>::PrintType>(m_Mean)
     << std::endl;
  os << indent << "StandardDeviation: "
     << static_cast<typename NumericTraits<double>::PrintType>(m_StandardDeviation)
     << std::endl;
}
} /* namespace itk */
