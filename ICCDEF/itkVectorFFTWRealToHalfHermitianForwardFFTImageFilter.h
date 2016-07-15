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
 *  itkVectorVectorFFTWRealToHalfHermitianForwardFFTImageFilter.h
 *  iccdefRegistrationNew
 *
 *  Created by Yongqiang Zhao on 5/6/09.
 *  Copyright 2009 UI. All rights reserved.
 *
 */
#ifndef __itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter_h
#define __itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <complex>
#include "itkVector.h"

#if defined(ITK_USE_FFTWF) || defined(ITK_USE_FFTWD)
#include "fftw3.h"
#endif

namespace itk
{
/** /class VectorFFTWRealToHalfHermitianForwardFFTImageFilter
 * /brief
 *
 * \ingroup
 */

template <class TPixel, unsigned int VDimension = 3>
class VectorFFTWRealToHalfHermitianForwardFFTImageFilter :
  public         ImageToImageFilter<Image<TPixel, VDimension>,
                                    Image<Vector<std::complex<typename TPixel::ValueType>, 3>, VDimension> >
{
public:
  /** Standard class typedefs. */
  typedef Image<TPixel, VDimension>                                              TInputImageType;
  typedef Image<Vector<std::complex<typename TPixel::ValueType>, 3>, VDimension> TOutputImageType;

  typedef VectorFFTWRealToHalfHermitianForwardFFTImageFilter    Self;
  typedef ImageToImageFilter<TInputImageType, TOutputImageType> Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorFFTWRealToHalfHermitianForwardFFTImageFilter,
               ImageToImageFilter);
  /** Image type typedef support. */
  typedef TInputImageType              ImageType;
  typedef typename ImageType::SizeType ImageSizeType;
  virtual void GenerateOutputInformation() ITK_OVERRIDE; // figure out allocation for output image

  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  virtual void EnlargeOutputRequestedRegion(DataObject *output) ITK_OVERRIDE;

  //
  // these should be defined in every FFT filter class
  virtual void GenerateData() ITK_OVERRIDE;  // generates output from input

protected:
  VectorFFTWRealToHalfHermitianForwardFFTImageFilter() : m_PlanComputed(false),
    m_LastImageSize(0),
    m_InputBuffer(0),
    m_OutputBuffer(0)
  {
  }

  ~VectorFFTWRealToHalfHermitianForwardFFTImageFilter()
  {
    if( m_PlanComputed )
      {
      fftwf_destroy_plan(this->m_Plan);
      delete [] this->m_InputBuffer;
      delete [] this->m_OutputBuffer;
      }
  }

  virtual bool FullMatrix();

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(VectorFFTWRealToHalfHermitianForwardFFTImageFilter);

  bool         m_PlanComputed;
  fftwf_plan   m_Plan;
  unsigned int m_LastImageSize;

  float *         m_InputBuffer;
  fftwf_complex * m_OutputBuffer;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter.hxx"
#endif

#endif // __itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter_h
