/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef __itkBlendImageFilter_hxx
#define __itkBlendImageFilter_hxx

#include "itkBlendImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace itk
{
/**
 * Constructor
 */
template < typename TInputImage, typename TOutputImage >
BlendImageFilter< TInputImage, TOutputImage >::BlendImageFilter()
  : m_Blend1( 1.0 )
  , m_Blend2( 1.0 )
{}

/**
 * GenerateData Performs the accumulation
 */
template < typename TInputImage, typename TOutputImage >
void
BlendImageFilter< TInputImage, TOutputImage >::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId )
{
  // We use dynamic_cast since inputs are stored as DataObjects.  The
  // ImageToImageFilter::GetInput(int) always returns a pointer to a
  // TInputImage so it cannot be used for the second input.
  InputImagePointer  inputPtr1 = dynamic_cast< TInputImage * >( this->ProcessObject::GetInput( 0 ) );
  InputImagePointer  inputPtr2 = dynamic_cast< TInputImage * >( this->ProcessObject::GetInput( 1 ) );
  OutputImagePointer outputPtr = dynamic_cast< TOutputImage * >( this->GetOutput( 0 ) );

  ImageRegionIterator< TInputImage > inputIt1( inputPtr1, outputRegionForThread );
  ImageRegionIterator< TInputImage > inputIt2( inputPtr2, outputRegionForThread );

  ImageRegionIterator< TOutputImage > outputIt( outputPtr, outputRegionForThread );

  ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  inputIt1.GoToBegin();
  inputIt2.GoToBegin();
  outputIt.GoToBegin();

  while ( !inputIt1.IsAtEnd() )
  {
    double acc = static_cast< double >( inputIt1.Get() ) * this->m_Blend1;
    acc += static_cast< double >( inputIt2.Get() ) * this->m_Blend2;
    outputIt.Set( static_cast< typename TOutputImage::PixelType >( acc ) );
    ++inputIt2;
    ++inputIt1;
    ++outputIt;
    progress.CompletedPixel(); // potential exception thrown here
  }
}

template < typename TInputImage, typename TOutputImage >
void
BlendImageFilter< TInputImage, TOutputImage >::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Blend1: " << this->m_Blend1 << std::endl;
  os << indent << "Blend2: " << this->m_Blend2 << std::endl;
}
} // end namespace itk

#endif
