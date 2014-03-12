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
#include "itkSpeckleNoiseImageFilter.h"
#include "itkThreadSafeMersenneTwisterRandomVariateGenerator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
SpeckleNoiseImageFilter<TInputImage, TOutputImage>
::SpeckleNoiseImageFilter()
{
  m_StandardDeviation = 1.0;
}

template <class TInputImage, class TOutputImage>
void
SpeckleNoiseImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
                        ThreadIdType threadId)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  InputImageConstPointer inputPtr = this->GetInput();
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  OutputImagePointer outputPtr = this->GetOutput(0);

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

  // choose the value of the gamma distribution so that the mean is 1 and the variance depend
  // on m_StandardDeviation
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  double theta = m_StandardDeviation * m_StandardDeviation;
  double k = 1 / theta;
  double floork = itk::Math::Floor<double>( k );
  double delta = k - floork;
  double v0 = itk::Math::e / ( itk::Math::e + delta );

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  while( !inputIt.IsAtEnd() )
    {
    // first generate the gamma distributed random variable
    // ref http://en.wikipedia.org/wiki/Gamma_distribution#Generating_gamma-distributed_random_variables
    double xi;
    double nu;

    do
      {
      double v1 = 1.0 - rand->GetVariateWithOpenUpperRange(); // open *lower* range -- (0,1]
      double v2 = 1.0 - rand->GetVariateWithOpenUpperRange();
      double v3 = 1.0 - rand->GetVariateWithOpenUpperRange();
      if( v1 <= v0 )
        {
        xi = vcl_pow( v2, 1 / delta );
        nu = v3 * vcl_pow( xi, delta - 1.0 );
        }
      else
        {
        xi = 1.0 - vcl_log( v2 );
        nu = v3 * vcl_exp( -xi );
        }
      }
    while( nu > vcl_exp( -xi ) * vcl_pow( xi, delta - 1.0 ) );

    double gamma = xi;
    for( int i = 0; i < floork; i++ )
      {
      gamma -= vcl_log( 1.0 - rand->GetVariateWithOpenUpperRange() );
      }
    gamma *= theta;
    // ok, so now apply multiplicative noise
    double out = gamma * inputIt.Get();
    // and clip the value to fit in the range (saturation)
    out = std::min( (double)NumericTraits<OutputImagePixelType>::max(), out );
    out = std::max( (double)NumericTraits<OutputImagePixelType>::NonpositiveMin(), out );
    outputIt.Set( (OutputImagePixelType) out  );
    ++inputIt;
    ++outputIt;
    progress.CompletedPixel();  // potential exception thrown here
    }

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
}

template <class TInputImage, class TOutputImage>
void
SpeckleNoiseImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os,
            Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "StandardDeviation: "
     << static_cast<typename NumericTraits<double>::PrintType>(this->GetStandardDeviation() )
     << std::endl;
}
} /* namespace itk */
