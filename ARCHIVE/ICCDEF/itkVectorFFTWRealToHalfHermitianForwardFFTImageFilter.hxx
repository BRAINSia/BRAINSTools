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
 *  itkVectorVectorFFTWRealToHalfHermitianForwardFFTImageFilter.hxx
 *  iccdefRegistrationNew
 *
 *  Created by Yongqiang Zhao on 5/6/09.
 *  Copyright 2009 UI. All rights reserved.
 *
 */

#ifndef __itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter_hxx
#define __itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter_hxx

#include "itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter.h"
#include "itkFFTWRealToHalfHermitianForwardFFTImageFilter.h"
#include <iostream>
#include "itkIndent.h"
#include "itkMetaDataObject.h"

namespace itk
{
/** INFO:  There should be compile time type checks so that
           if only ITK_USE_FFTWF is defined, then only floats are valid.
           and if ITK_USE_FFTWD is defined, then only doubles are valid.
*/

template <typename TPixel, unsigned int VDimension>
void
VectorFFTWRealToHalfHermitianForwardFFTImageFilter<TPixel, VDimension>::GenerateOutputInformation()
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
  // has no meaning in the result of an FFT.
  unsigned int                                i;
  const typename TInputImageType::SizeType &  inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TInputImageType::IndexType & inputStartIndex = inputPtr->GetLargestPossibleRegion().GetIndex();

  typename TOutputImageType::SizeType  outputSize;
  typename TOutputImageType::IndexType outputStartIndex;

  //
  // in 4.3.4 of the FFTW documentation, they indicate the size of
  // of a real-to-complex FFT is N * N ... + (N /2+1)
  //                              1   2        d
  // complex numbers.
  // static_cast prob. not necessary but want to make sure integer
  // division is used.
  outputSize[0] = static_cast<unsigned int>(inputSize[0]) / 2 + 1;
  outputStartIndex[0] = inputStartIndex[0];
  for (i = 1; i < TOutputImageType::ImageDimension; i++)
  {
    outputSize[i] = inputSize[i];
    outputStartIndex[i] = inputStartIndex[i];
  }
  //
  // the halving of the input size hides the actual size of the input.
  // to get the same size image out of the IFFT, need to send it as
  // Metadata.
  using SizeScalarType = typename TOutputImageType::SizeType::SizeValueType;
  itk::MetaDataDictionary & OutputDic = outputPtr->GetMetaDataDictionary();
  itk::EncapsulateMetaData<SizeScalarType>(OutputDic, std::string("FFT_Actual_RealImage_Size"), inputSize[0]);
  typename TOutputImageType::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(outputSize);
  outputLargestPossibleRegion.SetIndex(outputStartIndex);

  outputPtr->SetRegions(outputLargestPossibleRegion);
}

template <typename TPixel, unsigned int VDimension>
void
VectorFFTWRealToHalfHermitianForwardFFTImageFilter<TPixel, VDimension>::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  typename TInputImageType::Pointer input = const_cast<TInputImageType *>(this->GetInput());

  if (!input)
  {
    return;
  }

  input->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TPixel, unsigned int VDimension>
void
VectorFFTWRealToHalfHermitianForwardFFTImageFilter<TPixel, VDimension>::EnlargeOutputRequestedRegion(
  DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);

  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TPixel, unsigned int VDimension>
void
VectorFFTWRealToHalfHermitianForwardFFTImageFilter<TPixel, VDimension>::GenerateData()
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

  const typename TInputImageType::SizeType &  inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TOutputImageType::SizeType & outputSize = outputPtr->GetLargestPossibleRegion().GetSize();

  // figure out sizes
  // size of input and output aren't the same which is handled in the superclass,
  // sort of.
  // the input size and output size only differ in the fastest moving dimension
  unsigned int total_inputSize = 1;
  unsigned int total_outputSize = 1;
  for (unsigned i = 0; i < VDimension; i++)
  {
    total_inputSize *= inputSize[i];
    total_outputSize *= outputSize[i];
  }

  if (this->m_PlanComputed) // if we've already computed a plan
  {
    // if the image sizes aren't the same,
    // we have to compute the plan again
    if (this->m_LastImageSize != total_inputSize)
    {
      delete[] this->m_InputBuffer;
      delete[] this->m_OutputBuffer;
      fftwf_destroy_plan(this->m_Plan);
      this->m_PlanComputed = false;
    }
  }
  if (!this->m_PlanComputed)
  {
    this->m_InputBuffer = new float[total_inputSize * 3];
    this->m_OutputBuffer = new fftwf_complex[total_outputSize * 3];
    this->m_LastImageSize = total_inputSize;

    int * sizes = new int[VDimension];
    for (unsigned int i = 0; i < VDimension; i++)
    {
      sizes[(VDimension - 1) - i] = inputSize[i];
    }

    this->m_Plan = fftwf_plan_many_dft_r2c(VDimension,
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
  memcpy(this->m_InputBuffer, inputPtr->GetBufferPointer(), total_inputSize * sizeof(float) * 3);
  fftwf_execute(this->m_Plan);
  memcpy(outputPtr->GetBufferPointer(), this->m_OutputBuffer, total_outputSize * sizeof(fftwf_complex) * 3);
}

template <typename TPixel, unsigned int VDimension>
bool
VectorFFTWRealToHalfHermitianForwardFFTImageFilter<TPixel, VDimension>::FullMatrix()
{
  return false;
}
} // namespace itk

#endif // _itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter_hxx
