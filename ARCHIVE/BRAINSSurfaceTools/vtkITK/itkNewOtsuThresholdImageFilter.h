/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNewOtsuThresholdImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/16 19:59:56 $
  Version:   $Revision: 1.1.2.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNewOtsuThresholdImageFilter_h
#define __itkNewOtsuThresholdImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkFixedArray.h"

namespace itk
{

/** \class NewOtsuThresholdImageFilter
 * \brief Threshold an image using the Otsu Threshold
 *
 * This filter creates a binary thresholded image that separates an
 * image into foreground and background components. The filter
 * computes the threshold using the OtsuThresholdImageCalculator and
 * applies that theshold to the input image using the
 * BinaryThresholdImageFilter. The NunberOfHistogram bins can be set
 * for the Calculator. The InsideValue and OutsideValue can be set
 * for the BinaryThresholdImageFilter.
 *
 * \sa NewOtsuThresholdImageCalculator
 * \sa BinaryThresholdImageFilter
 * \ingroup IntensityImageFilters  Multithreaded
 */
template <typename TInputImage, typename TOutputImage>
class NewOtsuThresholdImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(NewOtsuThresholdImageFilter);

  /** Standard Self type alias */
  using Self = NewOtsuThresholdImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(NewOtsuThresholdImageFilter, ImageToImageFilter);

  /** Image pixel value type alias. */
  using InputPixelType = typename TInputImage::PixelType;
  using OutputPixelType = typename TOutputImage::PixelType;

  /** Image related type alias. */
  using InputImagePointer = typename TInputImage::Pointer;
  using OutputImagePointer = typename TOutputImage::Pointer;

  using InputSizeType = typename TInputImage::SizeType;
  using InputIndexType = typename TInputImage::IndexType;
  using InputImageRegionType = typename TInputImage::RegionType;
  using OutputSizeType = typename TOutputImage::SizeType;
  using OutputIndexType = typename TOutputImage::IndexType;
  using OutputImageRegionType = typename TOutputImage::RegionType;


  /** Image related type alias. */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Set the "outside" pixel value. The default value
   * NumericTraits<OutputPixelType>::ZeroValue(). */
  itkSetMacro(OutsideValue, OutputPixelType);

  /** Get the "outside" pixel value. */
  itkGetMacro(OutsideValue, OutputPixelType);

  /** Set the "inside" pixel value. The default value
   * NumericTraits<OutputPixelType>::max() */
  itkSetMacro(InsideValue, OutputPixelType);

  /** Get the "inside" pixel value. */
  itkGetMacro(InsideValue, OutputPixelType);

  /** Set/Get the number of histogram bins. Defaults is 128. */
  itkSetClampMacro(NumberOfHistogramBins, unsigned long, 1, NumericTraits<unsigned long>::max());
  itkGetMacro(NumberOfHistogramBins, unsigned long);

  itkSetMacro(Omega, double);
  itkGetMacro(Omega, double);

  /** Get the computed threshold. */
  itkGetMacro(Threshold, InputPixelType);


protected:
  NewOtsuThresholdImageFilter();
  ~NewOtsuThresholdImageFilter(){};
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateInputRequestedRegion() override;
  void
  GenerateData() override;

private:
  InputPixelType  m_Threshold;
  OutputPixelType m_InsideValue;
  OutputPixelType m_OutsideValue;
  unsigned long   m_NumberOfHistogramBins;
  double          m_Omega;
}; /// end of class

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkNewOtsuThresholdImageFilter.txx"
#endif

#endif
