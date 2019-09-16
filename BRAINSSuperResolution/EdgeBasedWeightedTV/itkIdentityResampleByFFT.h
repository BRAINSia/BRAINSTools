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
#ifndef itkIdentityResampleByFFT_h
#define itkIdentityResampleByFFT_h

#include "itkFixedArray.h"
#include "itkImageRegionIterator.h"
#include "itkImageToImageFilter.h"
#include "itkExtrapolateImageFunction.h"
#include "itkSize.h"
#include "itkDataObjectDecorator.h"


namespace itk
{
/** \class IdentityResampleByFFT
 * \brief Resample an image via a coordinate transform
 *
 * IdentityResampleByFFT resamples an existing image through some coordinate
 * transform, interpolating via some image function.  The class is templated
 * over the types of the input and output images.
 *
 * Note that the choice of interpolator function can be important.
 * This function is set via SetInterpolator().  The default is
 * LinearInterpolateImageFunction<InputImageType,
 * TInterpolatorPrecisionType>, which
 * is reasonable for ordinary medical images.  However, some synthetic
 * images have pixels drawn from a finite prescribed set.  An example
 * would be a mask indicating the segmentation of a brain into a small
 * number of tissue types.  For such an image, one does not want to
 * interpolate between different pixel values, and so
 * NearestNeighborInterpolateImageFunction< InputImageType,
 * TCoordRep > would be a better choice.
 *
 * If an sample is taken from outside the image domain, the default behavior is
 * to use a default pixel value.  If different behavior is desired, an
 * extrapolator function can be set with SetExtrapolator().
 *
 * Output information (spacing, size and direction) for the output
 * image should be set. This information has the normal defaults of
 * unit spacing, zero origin and identity direction. Optionally, the
 * output information can be obtained from a reference image. If the
 * reference image is provided and UseReferenceImage is On, then the
 * spacing, origin and direction of the reference image will be used.
 *
 * Since this filter produces an image which is a different size than
 * its input, it needs to override several of the methods defined
 * in ProcessObject in order to properly manage the pipeline execution model.
 * In particular, this filter overrides
 * ProcessObject::GenerateInputRequestedRegion() and
 * ProcessObject::GenerateOutputInformation().
 *
 * This filter is implemented as a multithreaded filter.  It provides a
 * ThreadedGenerateData() method for its implementation.
 * \warning For multithreading, the TransformPoint method of the
 * user-designated coordinate transform must be threadsafe.
 *
 * \ingroup GeometricTransform
 * \ingroup ITKImageGrid
 *
 * \wiki
 * \wikiexample{SimpleOperations/TranslationTransform,Translate an image}
 * \wikiexample{ImageProcessing/Upsampling,Upsampling an image}
 * \wikiexample{ImageProcessing/IdentityResampleByFFT,Resample (stretch or compress) an image}
 * \endwiki
 */
template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType = double>
class IdentityResampleByFFT : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(IdentityResampleByFFT);

  /** Standard class type alias. */
  using Self = IdentityResampleByFFT;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using InputImageRegionType = typename InputImageType::RegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(IdentityResampleByFFT, ImageToImageFilter);

  /** Number of dimensions. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

  /** base type for images of the current ImageDimension */
  using ImageBaseType = ImageBase<Self::ImageDimension>;

  /** Image size type alias. */
  using SizeType = Size<Self::ImageDimension>;

  /** Image index type alias. */
  using IndexType = typename TOutputImage::IndexType;

  /** Image point type alias. */
  using PointType = typename InterpolatorType::PointType;
  // using PointType = typename TOutputImage::PointType;

  /** Image pixel value type alias. */
  using PixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;

  using PixelConvertType = DefaultConvertPixelTraits<PixelType>;

  using PixelComponentType = typename PixelConvertType::ComponentType;

  /** Input pixel continuous index typdef */
  using ContinuousInputIndexType = ContinuousIndex<TTransformPrecisionType, ImageDimension>;

  /** Typedef to describe the output image region type. */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Image spacing,origin and direction type alias */
  using SpacingType = typename TOutputImage::SpacingType;
  using OriginPointType = typename TOutputImage::PointType;
  using DirectionType = typename TOutputImage::DirectionType;

  /** Typedef the reference image type to be the ImageBase of the OutputImageType */
  using ReferenceImageBaseType = ImageBase<ImageDimension>;

  itkSetMacro(Size, SizeType);
  itkGetConstReferenceMacro(Size, SizeType);

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void
  SetOutputSpacing(const double * values);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, OriginPointType);
  virtual void
  SetOutputOrigin(const double * values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro(OutputOrigin, OriginPointType);

  /** Set the output direciton cosine matrix. */
  itkSetMacro(OutputDirection, DirectionType);
  itkGetConstReferenceMacro(OutputDirection, DirectionType);

  /** Helper method to set the output parameters based on this image */
  void
  SetOutputParametersFromImage(const ImageBaseType * image);

  /** Set the start index of the output largest possible region.
   * The default is an index of all zeros. */
  itkSetMacro(OutputStartIndex, IndexType);

  /** Get the start index of the output largest possible region. */
  itkGetConstReferenceMacro(OutputStartIndex, IndexType);

  /** Set a reference image to use to define the output information.
   *  By default, output information is specificed through the
   *  SetOutputSpacing, Origin, and Direction methods.  Alternatively,
   *  this method can be used to specify an image from which to
   *  copy the information. UseReferenceImageOn must be set to utilize the
   *  reference image. */
  itkSetInputMacro(ReferenceImage, ReferenceImageBaseType);

  /** Get the reference image that is defining the output information. */
  itkGetInputMacro(ReferenceImage, ReferenceImageBaseType);

  /** Turn on/off whether a specified reference image should be used to define
   *  the output information. */
  itkSetMacro(UseReferenceImage, bool);
  itkBooleanMacro(UseReferenceImage);
  itkGetConstMacro(UseReferenceImage, bool);

  /** IdentityResampleByFFT produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  virtual void
  GenerateOutputInformation() override;

  /** IdentityResampleByFFT needs a different input requested region than
   * the output requested region.  As such, IdentityResampleByFFT needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void
  GenerateInputRequestedRegion() override;

  /** This method is used to set the state of the filter before
   * multi-threading. */
  virtual void
  BeforeThreadedGenerateData() override;

  /** This method is used to set the state of the filter after
   * multi-threading. */
  virtual void
  AfterThreadedGenerateData() override;

  /** Method Compute the Modified Time based on changed to the components. */
  ModifiedTimeType
  GetMTime(void) const override;

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro(OutputHasNumericTraitsCheck, (Concept::HasNumericTraits<PixelComponentType>));
  // End concept checking
#endif

protected:
  IdentityResampleByFFT();
  ~IdentityResampleByFFT() {}
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Override VeriyInputInformation() since this filter's inputs do
   * not need to occoupy the same physical space.
   *
   * \sa ProcessObject::VerifyInputInformation
   */
  void
  VerifyInputInformation() const override
  {}

  /** IdentityResampleByFFT can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior
   * to calling ThreadedGenerateData().
   * ThreadedGenerateData can only write to the portion of the output image
   * specified by the parameter "outputRegionForThread"
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

  /** Default implementation for resampling that works for any
   * transformation type. */
  virtual void
  NonlinearThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkIdentityResampleByFFT.hxx"
#endif

#endif
