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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef _itkAverageImageFilter_hxx
#define _itkAverageImageFilter_hxx

#include "itkAverageImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkNumericTraits.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
void
AverageImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

template< typename TInputImage, typename TOutputImage >
void
AverageImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread,
                        ThreadIdType itkNotUsed(threadId))
{
  typedef ImageRegionConstIterator< TInputImage > IteratorType;
  typedef ImageRegionIterator< TOutputImage> OutIteratorType;

  typename TOutputImage::Pointer output = this->GetOutput();

  // Record the number of input images.
  const unsigned int numberOfInputFiles = this->GetNumberOfInputs();

  //  create and initialize all input image iterators
  IteratorType *it = new IteratorType[numberOfInputFiles];
  for ( unsigned int i = 0; i < numberOfInputFiles; ++i)
    {
    it[i] = IteratorType( this->GetInput( i ),
                          outputRegionForThread );
    }

  OutIteratorType out = OutIteratorType( output, outputRegionForThread );
  for( out.GoToBegin(); !out.IsAtEnd(); ++out )
    {
    typename NumericTraits<OutputPixelType>::RealType sum = it[0].Get();
    ++it[0];

    for( unsigned int i = 1; i < numberOfInputFiles; ++i)
      {
      sum += it[i].Get();
      ++(it[i]);
      }

    sum *= (1.0 / numberOfInputFiles);
    // instead of casting to output type, we should support rounding if the
    // output type is an integer type. Unfortunately, this does not seem to
    // be supported by the NumericTraits so far.
    out.Set( (OutputPixelType) sum );
    }

  delete[] it;
}

} // end namespace itk

#endif
