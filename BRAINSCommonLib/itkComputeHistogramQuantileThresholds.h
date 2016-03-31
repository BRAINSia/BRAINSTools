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
#ifndef __itkComputeHistogramQuantileThresholds_h
#define __itkComputeHistogramQuantileThresholds_h

#include <itkImage.h>
#include <itkNumericTraits.h>

namespace itk
{
/**
  * \class ComputeHistogramQuantileThresholds
  * \author Hans J. Johnson
  *
  * This filter just computes Histogram Quantile Thresholds.  It does not apply
  *the thresholds.
  *
  */
template <class TInputImage, class TMaskImage>
class ComputeHistogramQuantileThresholds :
  public         Object
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::ConstPointer InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputPixelType;

  typedef ComputeHistogramQuantileThresholds Self;
  typedef Object                             Superclass;
  typedef SmartPointer<Self>                 Pointer;
  typedef typename TMaskImage::PixelType     MaskPixelType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ComputeHistogramQuantileThresholds, Object);

  /** set Quantile Threshold */
  itkSetMacro(QuantileLowerThreshold, double);
  itkGetConstMacro(QuantileLowerThreshold, double);
  itkSetMacro(QuantileUpperThreshold, double);
  itkGetConstMacro(QuantileUpperThreshold, double);

  itkGetConstMacro(LowerIntensityThresholdValue, typename InputImageType::PixelType);
  itkGetConstMacro(UpperIntensityThresholdValue, typename InputImageType::PixelType);
  itkGetConstMacro(NumberOfValidHistogramsEntries, unsigned int);

  itkGetConstObjectMacro(Image, InputImageType);
  itkSetConstObjectMacro(Image, InputImageType);

  itkSetMacro(ImageMin, typename TInputImage::PixelType);
  itkGetConstMacro(ImageMin, typename TInputImage::PixelType);
  itkSetMacro(ImageMax, typename TInputImage::PixelType);
  itkGetConstMacro(ImageMax, typename TInputImage::PixelType);

  itkGetConstObjectMacro(BinaryPortionImage, TMaskImage);
  itkSetObjectMacro(BinaryPortionImage, TMaskImage);

  void Calculate();

protected:
  ComputeHistogramQuantileThresholds();
  ~ComputeHistogramQuantileThresholds();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  void ImageMinMax(InputPixelType & min, InputPixelType & max);

  InputImagePointer m_Image;
  typename TMaskImage::Pointer m_BinaryPortionImage;

  double m_QuantileLowerThreshold;
  double m_QuantileUpperThreshold;
  unsigned int m_NumberOfValidHistogramsEntries;

  typename TInputImage::PixelType m_ImageMin;
  typename TInputImage::PixelType m_ImageMax;

  typename InputImageType::PixelType m_LowerIntensityThresholdValue;
  typename InputImageType::PixelType m_UpperIntensityThresholdValue;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeHistogramQuantileThresholds.hxx"
#endif

#endif
