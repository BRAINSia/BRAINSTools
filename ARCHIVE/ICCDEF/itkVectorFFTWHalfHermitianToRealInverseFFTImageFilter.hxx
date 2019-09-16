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
 *  itkVectorFFTComplexConjugateToRealImageFilter.hxx
 *  iccdefRegistrationNew
 *
 *  Created by Yongqiang Zhao on 5/6/09.
 *  Copyright 2009 UI. All rights reserved.
 *
 */

#ifndef __itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter_hxx
#define __itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter_hxx

#include "itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter.h"
#include "itkFFTWHalfHermitianToRealInverseFFTImageFilter.h"
#include <iostream>
#include "itkIndent.h"
#include "itkMetaDataObject.h"
#include "itkImageRegionIterator.h"

namespace itk
{
template <typename TPixel, unsigned int Dimension>
void
VectorFFTWHalfHermitianToRealInverseFFTImageFilter<TPixel, Dimension>::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  //
  // If this implementation returns a full result
  // instead of a 'half-complex' matrix, then none of this
  // is necessary
  if (this->FullMatrix())
  {
    return;
  }

  // get pointers to the input and output
  typename TInputImageType::ConstPointer inputPtr = this->GetInput();
  typename TOutputImageType::Pointer     outputPtr = this->GetOutput();

  if (!inputPtr || !outputPtr)
  {
    return;
  }

  //
  // This is all based on the same function in itk::ShrinkImageFilter
  // ShrinkImageFilter also modifies the image spacing, but spacing
  // has no meaning in the result of an FFT. For an IFFT, since the
  // spacing is propagated to the complex result, we can use the spacing
  // from the input to propagate back to the output.
  unsigned int                                i;
  const typename TInputImageType::SizeType &  inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TInputImageType::IndexType & inputStartIndex = inputPtr->GetLargestPossibleRegion().GetIndex();

  typename TOutputImageType::SizeType  outputSize;
  typename TOutputImageType::IndexType outputStartIndex;

  //
  // in 4.3.4 of the FFT documentation, they indicate the size of
  // of a real-to-complex FFT is N * N ... + (N /2+1)
  //                              1   2        d
  // complex numbers.
  // going from complex to real, you know the output is at least
  // twice the size in the last dimension as the input, but it might
  // be 2*size+1.  Consequently, the output of the FFT:R2C operation
  //
  MetaDataDictionary & InputDic = const_cast<MetaDataDictionary &>(inputPtr->GetMetaDataDictionary());

  using SizeScalarType = typename TInputImageType::SizeType::SizeValueType;

  SizeScalarType x = 0;

  outputSize[0] = (inputSize[0] - 1) * 2;
  if (this->ActualXDimensionIsOdd())
  {
    outputSize[0]++;
  }
  // backwards compatible/deprecated version
  if (ExposeMetaData<SizeScalarType>(InputDic, std::string("FFT_Actual_RealImage_Size"), x))
  {
    outputSize[0] = x;
  }

  outputStartIndex[0] = inputStartIndex[0];
  for (i = 1; i < TOutputImageType::ImageDimension; i++)
  {
    outputSize[i] = inputSize[i];
    outputStartIndex[i] = inputStartIndex[i];
  }
  typename TOutputImageType::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(outputSize);
  outputLargestPossibleRegion.SetIndex(outputStartIndex);

  outputPtr->SetRegions(outputLargestPossibleRegion);
  //  outputPtr->SetBufferedRegion(outputLargestPossibleRegion);
  //  outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
  //  outputPtr->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TPixel, unsigned int Dimension>
void
VectorFFTWHalfHermitianToRealInverseFFTImageFilter<TPixel, Dimension>::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  typename TInputImageType::Pointer inputPtr = const_cast<TInputImageType *>(this->GetInput());
  if (inputPtr)
  {
    inputPtr->SetRequestedRegionToLargestPossibleRegion();
  }
}

template <typename TPixel, unsigned int VDimension>
void
VectorFFTWHalfHermitianToRealInverseFFTImageFilter<TPixel, VDimension>::GenerateData()
{
  // get pointers to the input and output
  typename TInputImageType::ConstPointer inputPtr = this->GetInput();
  typename TOutputImageType::Pointer     outputPtr = this->GetOutput();

  if (!inputPtr || !outputPtr)
  {
    return;
  }

  // allocate output buffer memory
  outputPtr->SetBufferedRegion(outputPtr->GetRequestedRegion());
  outputPtr->Allocate();

  const typename TInputImageType::SizeType &  outputSize = outputPtr->GetLargestPossibleRegion().GetSize();
  const typename TOutputImageType::SizeType & inputSize = inputPtr->GetLargestPossibleRegion().GetSize();

  // figure out sizes
  // size of input and output aren't the same which is handled in the superclass,
  // sort of.
  // the input size and output size only differ in the fastest moving dimension
  unsigned int total_outputSize = 1;
  unsigned int total_inputSize = 1;
  for (unsigned i = 0; i < VDimension; i++)
  {
    total_outputSize *= outputSize[i];
    total_inputSize *= inputSize[i];
  }

  if (this->m_PlanComputed) // if we've already computed a plan
  {
    // if the image sizes aren't the same,
    // we have to compute the plan again
    if (this->m_LastImageSize != total_outputSize)
    {
      delete[] this->m_InputBuffer;
      delete[] this->m_OutputBuffer;
      fftwf_destroy_plan(this->m_Plan);
      this->m_PlanComputed = false;
    }
  }
  // either plan never computed, or need to re-compute
  if (!this->m_PlanComputed)
  {
    // if we've never computed the plan, or we need to redo it
    this->m_InputBuffer = new fftwf_complex[total_inputSize * 3];
    this->m_OutputBuffer = new float[total_outputSize * 3];
    this->m_LastImageSize = total_outputSize;

    int * sizes = new int[VDimension];
    for (unsigned int i = 0; i < VDimension; i++)
    {
      sizes[(VDimension - 1) - i] = outputSize[i];
    }

    this->m_Plan = fftwf_plan_many_dft_c2r(VDimension,
                                           sizes,
                                           3,
                                           this->m_InputBuffer,
                                           nullptr,
                                           3,
                                           1,
                                           this->m_OutputBuffer,
                                           nullptr,
                                           3,
                                           1,
                                           FFTW_MEASURE | FFTW_DESTROY_INPUT);
    delete[] sizes;
    this->m_PlanComputed = true;
  }
  // copy the input, because it may be destroyed by computing the plan
  memcpy(this->m_InputBuffer, inputPtr->GetBufferPointer(), total_inputSize * sizeof(fftwf_complex) * 3);
  fftwf_execute(this->m_Plan);

  // copy the output
  memcpy(outputPtr->GetBufferPointer(), this->m_OutputBuffer, total_outputSize * sizeof(float) * 3);

  using IteratorType = ImageRegionIterator<TOutputImageType>;

  IteratorType it(outputPtr, outputPtr->GetLargestPossibleRegion());

  while (!it.IsAtEnd())
  {
    it.Set(it.Value() / total_outputSize);
    ++it;
  }
}

template <typename TPixel, unsigned int VDimension>
bool
VectorFFTWHalfHermitianToRealInverseFFTImageFilter<TPixel, VDimension>::FullMatrix()
{
  return false;
}
} // namespace itk
#endif // _itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter_hxx
