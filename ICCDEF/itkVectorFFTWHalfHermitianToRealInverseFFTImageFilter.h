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
 *  itkVectorFFTComplexConjugateToRealImageFilter.h
 *  iccdefRegistrationNew
 *
 *  Created by Yongqiang Zhao on 5/6/09.
 *  Copyright 2009 UI. All rights reserved.
 *
 */

#ifndef __itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter_h
#define __itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <complex>

//
// FFTWCommon defines proxy classes based on data types
#if defined(ITK_USE_FFTWF) || defined(ITK_USE_FFTWD)
#include "fftw3.h"
#endif

namespace itk
{
template <typename TPixel, unsigned int VDimension = 3>
class VectorFFTWHalfHermitianToRealInverseFFTImageFilter :
  public ImageToImageFilter<Image<Vector<std::complex<typename TPixel::ValueType>, 3>, VDimension>,
                            Image<TPixel, VDimension> >
{
public:
  /** Standard class typedefs.*/
  typedef Image<Vector<std::complex<typename TPixel::ValueType>, 3>, VDimension> TInputImageType;
  typedef Image<TPixel, VDimension>                                              TOutputImageType;

  typedef VectorFFTWHalfHermitianToRealInverseFFTImageFilter    Self;
  typedef ImageToImageFilter<TInputImageType, TOutputImageType> Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;
  //

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorFFTWHalfHermitianToRealInverseFFTImageFilter,
               ImageToImageFilter);

  /** Image type typedef support. */
  typedef TInputImageType              ImageType;
  typedef typename ImageType::SizeType ImageSizeType;
  virtual void GenerateOutputInformation() ITK_OVERRIDE; // figure out allocation for output image

  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  //
  // these should be defined in every FFT filter class
  virtual void GenerateData() ITK_OVERRIDE;  // generates output from input

  virtual bool FullMatrix();

  void SetActualXDimensionIsOdd(bool isodd)
  {
    m_ActualXDimensionIsOdd = isodd;
  }

  void SetActualXDimensionIsOddOn()
  {
    this->SetActualXDimensionIsOdd(true);
  }

  void SetActualXDimensionIsOddOff()
  {
    this->SetActualXDimensionIsOdd(false);
  }

  bool ActualXDimensionIsOdd()
  {
    return m_ActualXDimensionIsOdd;
  }

protected:
  VectorFFTWHalfHermitianToRealInverseFFTImageFilter() : m_PlanComputed(false),
    m_LastImageSize(0),
    m_InputBuffer(0),
    m_OutputBuffer(0),
    m_ActualXDimensionIsOdd(false)
  {
  }

  virtual ~VectorFFTWHalfHermitianToRealInverseFFTImageFilter()
  {
    if( m_PlanComputed )
      {
      fftwf_destroy_plan(this->m_Plan);
      delete [] m_InputBuffer;
      delete [] m_OutputBuffer;
      }
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(VectorFFTWHalfHermitianToRealInverseFFTImageFilter);

  bool         m_PlanComputed;
  fftwf_plan   m_Plan;
  unsigned int m_LastImageSize;
  // local storage needed to keep fftw from scribbling on
  fftwf_complex * m_InputBuffer;
  float *         m_OutputBuffer;
  bool            m_ActualXDimensionIsOdd;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter.hxx"
#endif

#endif // __itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter_h
