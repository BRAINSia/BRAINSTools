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
#ifndef __itkGtractInverseDisplacementFieldImageFilter_h
#define __itkGtractInverseDisplacementFieldImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkKernelTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
/** \class GtractInverseDisplacementFieldImageFilter
 * \brief Computes the inverse of a deformation field.
 *
 * GtractInverseDisplacementFieldImageFilter takes a deformation field as input and
 * computes the deformation field that is its inverse. If the input deformation
 * field was mapping coordinates from a space A into a space B, the output of
 * this filter will map coordinates from the space B into the space A.
 *
 * Given that both the input and output deformation field are represented as
 * discrete images with pixel type vector, the inverse will be only an
 * estimation and will probably not correspond to a perfect inverse.  The
 * precision of the inverse can be improved at the price of increasing the
 * computation time and memory consumption in this filter.
 *
 * The method used for computing the inverse deformation field is to subsample
 * the input field using a regular grid and create Kerned-Base Spline in which
 * the reference landmarks are the coordinates of the deformed point and the
 * target landmarks are the negative of the displacement vectors. The
 * kernel-base spline is then used for regularly sampling the output space and
 * recover vector values for every std::single pixel.
 *
 * The subsampling factor used for the regular grid of the input field will
 * determine the number of landmarks in the KernelBased spline and therefore it
 * will have a dramatic effect on both the precision of output deformation
 * field and the computational time required for the filter to complete the
 * estimation. A large subsampling factor will result in few landmarks in the
 * KernelBased spline, therefore on fast computation and low precision.  A
 * small subsampling factor will result in a large number of landmarks in the
 * KernelBased spline, therefore a large memory consumption, long computation
 * time and high precision for the inverse estimation.
 *
 * This filter std::expects both the input and output images to be of pixel type
 * Vector.
 *
 * \ingroup ImageToImageFilter
 */
template <typename TInputImage, typename TOutputImage>
class GtractInverseDisplacementFieldImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = GtractInverseDisplacementFieldImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GtractInverseDisplacementFieldImageFilter, ImageToImageFilter);

  /** Number of dimensions. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** Transform type alias.
   *
   * \todo Check that input and output images have the same number of
   * dimensions; this is required for consistency.  */
  using KernelTransformType = KernelTransform<double, Self::ImageDimension>;
  using KernelTransformPointerType = typename KernelTransformType::Pointer;

  /** Image size type alias. */
  using SizeType = typename OutputImageType::SizeType;

  /** Image index type alias. */
  using IndexType = typename OutputImageType::IndexType;

  /** Image pixel value type alias. */
  using OutputPixelType = typename TOutputImage::PixelType;
  using OutputPixelComponentType = typename OutputPixelType::ValueType;

  /** Typedef to describe the output image region type. */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Image spacing type alias */
  using SpacingType = typename TOutputImage::SpacingType;
  using OriginPointType = typename TOutputImage::PointType;

  /** Image direction type alias */
  using DirectionType = typename TOutputImage::DirectionType;

  /** Set the coordinate transformation.
   * Set the KernelBase spline used for resampling the deformation grid.
   * */
  itkSetObjectMacro(KernelTransform, KernelTransformType);

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro(KernelTransform, KernelTransformType);

  /** Set the size of the output image. */
  itkSetMacro(Size, SizeType);

  /** Get the size of the output image. */
  itkGetConstReferenceMacro(Size, SizeType);

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void
  SetOutputSpacing(const double * spacing);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, OriginPointType);
  virtual void
  SetOutputOrigin(const double * origin);

  /** Get the output image origin. */
  itkGetConstReferenceMacro(OutputOrigin, OriginPointType);

  /** Set the output image direction. */
  itkSetMacro(OutputDirection, DirectionType);

  /** Get the output image direction. */
  itkGetConstReferenceMacro(OutputDirection, DirectionType);

  /** Set/Get the factor used for subsampling the input deformation field.  A
   * large value in this factor will produce a fast computation of the inverse
   * field but with low precision. A small value of this factor will produce a
   * precise computation of the inverse field at the price of large memory
   * consumption and long computational time. */
  itkSetMacro(SubsamplingFactor, unsigned int);
  itkGetConstMacro(SubsamplingFactor, unsigned int);

  /** GtractInverseDisplacementFieldImageFilter produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  void
  GenerateOutputInformation() override;

  /** GtractInverseDisplacementFieldImageFilter needs a different input requested region than
   * the output requested region.  As such, GtractInverseDisplacementFieldImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  void
  GenerateInputRequestedRegion() override;

  /** Method Compute the Modified Time based on changed to the components. */
  unsigned long
  GetMTime() const override;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck, (Concept::HasNumericTraits<OutputPixelComponentType>));
  /** End concept checking */
#endif
protected:
  GtractInverseDisplacementFieldImageFilter();
  ~GtractInverseDisplacementFieldImageFilter() override {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /**
   * GenerateData() computes the internal KernelBase spline and resamples
   * the deformation field.
   */
  void
  GenerateData() override;

  /** Subsample the input deformation field and generate the
   *  landmarks for the kernel base spline
   */
  void
  PrepareKernelBaseSpline();

private:
  GtractInverseDisplacementFieldImageFilter(const Self &); // purposely not
                                                           // implemented
  void
  operator=(const Self &); // purposely not

  // implemented

  SizeType                   m_Size;            // Size of the output image
  KernelTransformPointerType m_KernelTransform; // Coordinate transform to
                                                // use
  SpacingType     m_OutputSpacing;              // output image spacing
  OriginPointType m_OutputOrigin;               // output image origin
  DirectionType   m_OutputDirection;            // output image direction

  unsigned int m_SubsamplingFactor; // factor to subsample the
                                    // input field.
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkGtractInverseDisplacementFieldImageFilter.hxx"
#endif

#endif // __itkGtractInverseDisplacementFieldImageFilter_h
