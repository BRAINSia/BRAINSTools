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
/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: itkOtsuHistogramMatchingImageFilter.h,v $
 *  Language:  C++
 *  Date:      $Date: 2009-05-02 05:43:54 $
 *  Version:   $Revision: 1.13 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __itkOtsuHistogramMatchingImageFilter_h
#define __itkOtsuHistogramMatchingImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkHistogram.h"
#include "BRAINSFitUtils.h"
#include "BRAINSTypes.h"
#include "vnl/vnl_matrix.h"

namespace itk
{
/** \class OtsuHistogramMatchingImageFilter
 * \brief Normalize the grayscale values between two image by histogram
 * matching.
 *
 * OtsuHistogramMatchingImageFilter normalizes the grayscale values of a source
 * image based on the grayscale values of a reference image.
 * This filter uses a histogram matching technique where the histograms of the
 * two images are matched only at a specified number of quantile values.
 *
 * This filter was inspired by the HistogramMatchingImagFilter that was
 * orginally designed to normalize MR images of the same
 * MR protocol and same body part. This algorithm is a specialization that
 * takes advantage of the fact that the algorihtms work best if background
 * pixels are excluded from both the source and reference histograms.
 * The background exclusion method uses the Otsu threshold to exclude all
 * background voxles whose grayscale values are smaller than Otsu threshold.
 *
 * The source image can be set via either SetInput() or SetSourceImage().
 * The reference image can be set via SetReferenceImage().
 *
 * SetNumberOfHistogramLevels() sets the number of bins used when
 * creating histograms of the source and reference images.
 * SetNumberOfMatchPoints() governs the number of quantile values to be
 * matched.
 *
 * This filter assumes that both the source and reference are of the same
 * type and that the input and output image type have the same number of
 * dimension and have scalar pixel types.
 *
 * \ingroup IntensityImageFilters Multithreaded
 *
 */
/* THistogramMeasurement -- The precision level for which to do
 * HistogramMeasurmenets */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement = typename TInputImage::PixelType>
class OtsuHistogramMatchingImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(OtsuHistogramMatchingImageFilter);

  /** Standard class type alias. */
  using Self = OtsuHistogramMatchingImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(OtsuHistogramMatchingImageFilter);

  /** ImageDimension enumeration. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Typedef to describe the output image region type. */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Inherited type alias. */
  using InputImageType = typename Superclass::InputImageType;
  using InputImagePointer = typename Superclass::InputImagePointer;
  using InputImageConstPointer = typename Superclass::InputImageConstPointer;
  using OutputImageType = typename Superclass::OutputImageType;
  using OutputImagePointer = typename Superclass::OutputImagePointer;

  /** Pixel related type alias. */
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;

  /** Histogram related type alias. */
  using HistogramType = Statistics::Histogram<THistogramMeasurement>;
  using HistogramPointer = typename HistogramType::Pointer;

  /** Set/Get the source image. */
  void
  SetSourceImage(const InputImageType * source)
  {
    this->SetInput(source);
  }

  const InputImageType *
  GetSourceImage()
  {
    return this->GetInput();
  }

  /** Set/Get the reference image. */
  void
  SetReferenceImage(const InputImageType * reference);

  const InputImageType *
  GetReferenceImage();

  itkSetObjectMacro(SourceMask, SpatialObjectType);
  itkSetObjectMacro(ReferenceMask, SpatialObjectType);

  itkGetConstObjectMacro(SourceMask, SpatialObjectType);
  itkGetConstObjectMacro(ReferenceMask, SpatialObjectType);

  /** Set/Get the number of histogram levels used. */
  itkSetMacro(NumberOfHistogramLevels, unsigned long);
  itkGetConstMacro(NumberOfHistogramLevels, unsigned long);

  /** Set/Get the number of match points used. */
  itkSetMacro(NumberOfMatchPoints, unsigned long);
  itkGetConstMacro(NumberOfMatchPoints, unsigned long);

  /** This filter requires all of the input to be in the buffer. */
  void
  GenerateInputRequestedRegion() override;

  /** Methods to get the histograms of the source, reference, and
   * output. Objects are only valid after Update() has been called
   * on this filter. */
  itkGetConstObjectMacro(SourceHistogram, HistogramType);
  itkGetConstObjectMacro(ReferenceHistogram, HistogramType);
  itkGetConstObjectMacro(OutputHistogram, HistogramType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(IntConvertibleToInputCheck, (Concept::Convertible<int, InputPixelType>));
  itkConceptMacro(SameDimensionCheck, (Concept::SameDimension<ImageDimension, OutputImageDimension>));
  itkConceptMacro(DoubleConvertibleToInputCheck, (Concept::Convertible<double, InputPixelType>));
  itkConceptMacro(DoubleConvertibleToOutputCheck, (Concept::Convertible<double, OutputPixelType>));
  itkConceptMacro(InputConvertibleToDoubleCheck, (Concept::Convertible<InputPixelType, double>));
  itkConceptMacro(OutputConvertibleToDoubleCheck, (Concept::Convertible<OutputPixelType, double>));
  itkConceptMacro(SameTypeCheck, (Concept::SameType<InputPixelType, OutputPixelType>));
  /** End concept checking */
#endif
protected:
  /** Override VeriyInputInformation() since this filter does not expect
   * the input images to occupy the same physical space.
   *
   * \sa ProcessObject::VerifyInputInformation
   */
  void
  VerifyInputInformation() const override
  {}

  OtsuHistogramMatchingImageFilter();
  ~OtsuHistogramMatchingImageFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  BeforeThreadedGenerateData() override;

  void
  AfterThreadedGenerateData() override;

  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

  /** Compute min, max and mean of an image. */
  void
  ComputeMinMaxMean(const InputImageType *  image,
                    THistogramMeasurement & minValue,
                    THistogramMeasurement & maxValue,
                    THistogramMeasurement & meanValue);

  /** Construct a histogram from an image. */
  void
  ConstructHistogram(const InputImageType *                      image,
                     const typename SpatialObjectType::Pointer & mask,
                     HistogramType *                             histogram,
                     const THistogramMeasurement                 minValue,
                     const THistogramMeasurement                 maxValue);

private:
  unsigned long m_NumberOfHistogramLevels{ 256 };
  unsigned long m_NumberOfMatchPoints{ 1 };
  bool          m_ThresholdAtMeanIntensity;

  InputPixelType  m_SourceIntensityThreshold;
  InputPixelType  m_ReferenceIntensityThreshold;
  OutputPixelType m_OutputIntensityThreshold;

  THistogramMeasurement m_SourceMinValue;
  THistogramMeasurement m_SourceMaxValue;
  THistogramMeasurement m_SourceMeanValue;
  THistogramMeasurement m_ReferenceMinValue;
  THistogramMeasurement m_ReferenceMaxValue;
  THistogramMeasurement m_ReferenceMeanValue;
  THistogramMeasurement m_OutputMinValue;
  THistogramMeasurement m_OutputMaxValue;
  THistogramMeasurement m_OutputMeanValue;

  HistogramPointer m_SourceHistogram;
  HistogramPointer m_ReferenceHistogram;
  HistogramPointer m_OutputHistogram;

  using TableType = vnl_matrix<double>;
  TableType m_QuantileTable;

  using GradientArrayType = vnl_vector<double>;
  GradientArrayType m_Gradients;
  double            m_LowerGradient{ 0.0 };
  double            m_UpperGradient{ 0.0 };

  typename SpatialObjectType::Pointer m_SourceMask;
  typename SpatialObjectType::Pointer m_ReferenceMask;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkOtsuHistogramMatchingImageFilter.hxx"
#endif

#endif
