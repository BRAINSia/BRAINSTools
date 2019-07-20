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
#ifndef __itkBOBFFilter_hxx
#define __itkBOBFFilter_hxx

#include "itkBOBFFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include <vector>
#include "itkNeighborhoodConnectedImageFilter.h"

namespace itk
{
/*
 *
 */
template < typename TInputImage, typename TOutputImage >
BOBFFilter< TInputImage, TOutputImage >::BOBFFilter()
  : m_Lower( NumericTraits< InputPixelType >::NonpositiveMin() )
  , m_Upper( NumericTraits< InputPixelType >::max() )
  , m_ReplaceValue( NumericTraits< OutputPixelType >::OneValue() )
{
  this->SetNumberOfRequiredInputs( 2 );
  m_Seed.Fill( 0 );
  m_Radius.Fill( 1 );
}

/*
 *
 */
template < typename TInputImage, typename TOutputImage >
void
BOBFFilter< TInputImage, TOutputImage >::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

/*
 *
 */
template < typename TInputImage, typename TOutputImage >
void
BOBFFilter< TInputImage, TOutputImage >::SetInputMask(
  const typename BOBFFilter< TInputImage, TOutputImage >::InputImageType * image )
{
  this->SetNthInput( 1, const_cast< TInputImage * >( image ) );
}

template < typename TInputImage, typename TOutputImage >
const typename BOBFFilter< TInputImage, TOutputImage >::InputImageType *
BOBFFilter< TInputImage, TOutputImage >::GetInputMask()
{
  return static_cast< const TInputImage * >( this->ProcessObject::GetInput( 1 ) );
}

template < typename TInputImage, typename TOutputImage >
void
BOBFFilter< TInputImage, TOutputImage >::GenerateData()
{
  OutputImagePointer     OutputPtr = this->GetOutput();
  InputImageConstPointer InputImage = this->GetInputImage();
  InputImageConstPointer InputMask = this->GetInputMask();

  /*Allocate the output*/
  OutputPtr->SetRequestedRegion( InputImage->GetRequestedRegion() );
  OutputPtr->SetBufferedRegion( InputImage->GetBufferedRegion() );
  OutputPtr->SetLargestPossibleRegion( InputImage->GetLargestPossibleRegion() );
  OutputPtr->CopyInformation( InputImage );
  OutputPtr->Allocate();

  using InputIterator = ImageRegionConstIterator< TInputImage >;
  using OutputIterator = ImageRegionIterator< TOutputImage >;

  OutputIterator outItr( OutputPtr, OutputPtr->GetLargestPossibleRegion() );

  InputIterator ImgItr( InputImage, InputImage->GetLargestPossibleRegion() );

  InputIterator MskItr( InputMask, InputMask->GetLargestPossibleRegion() );
  for ( ImgItr.GoToBegin(), MskItr.GoToBegin(), outItr.GoToBegin(); !ImgItr.IsAtEnd(); ++ImgItr, ++MskItr, ++outItr )
  {
    if ( MskItr.Get() == 0 )
    {
      outItr.Set( m_ReplaceValue );
    }
    else
    {
      outItr.Set( ImgItr.Get() );
    }
  }
}
} // end namespace itk

#endif // _itkBOBFFilter_hxx
