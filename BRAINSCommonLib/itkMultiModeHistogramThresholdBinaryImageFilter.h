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
#ifndef __itkMultiModeHistogramThresholdBinaryImageFilter_h
#define __itkMultiModeHistogramThresholdBinaryImageFilter_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>
#include <itkArray.h>

namespace itk
{
/**
  * \author Hans J. Johnson
  *
  * This filter
  *
  */
template <class TInputImage, class TOutputImage = Image<unsigned short, TInputImage::ImageDimension> >
class MultiModeHistogramThresholdBinaryImageFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::ConstPointer InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputPixelType;

  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef MultiModeHistogramThresholdBinaryImageFilter        Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef TOutputImage                                        IntegerImageType;
  typedef typename IntegerImageType::PixelType                IntegerPixelType;

  typedef Array<double> ThresholdArrayType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiModeHistogramThresholdBinaryImageFilter, ImageToImageFilter);

  itkSetMacro(LinearQuantileThreshold, double);
  itkGetConstMacro(LinearQuantileThreshold, double);

  /** set Quantile Threshold Arrays */
  itkSetMacro(QuantileLowerThreshold, ThresholdArrayType);
  itkGetConstMacro(QuantileLowerThreshold, ThresholdArrayType);
  itkSetMacro(QuantileUpperThreshold, ThresholdArrayType);
  itkGetConstMacro(QuantileUpperThreshold, ThresholdArrayType);

  itkGetConstObjectMacro(BinaryPortionImage, IntegerImageType);
  itkSetObjectMacro(BinaryPortionImage, IntegerImageType);

  itkSetMacro(InsideValue, IntegerPixelType);
  itkGetConstMacro(InsideValue, IntegerPixelType);
  itkSetMacro(OutsideValue, IntegerPixelType);
  itkGetConstMacro(OutsideValue, IntegerPixelType);
protected:
  MultiModeHistogramThresholdBinaryImageFilter();
  ~MultiModeHistogramThresholdBinaryImageFilter() override;
  void PrintSelf(std::ostream & os, Indent indent) const override;

  void GenerateData() override;

private:
  ThresholdArrayType m_QuantileLowerThreshold;
  ThresholdArrayType m_QuantileUpperThreshold;
  double             m_LinearQuantileThreshold;

  typename IntegerImageType::Pointer m_BinaryPortionImage;

  IntegerPixelType m_InsideValue;
  IntegerPixelType m_OutsideValue;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiModeHistogramThresholdBinaryImageFilter.hxx"
#endif

#endif
