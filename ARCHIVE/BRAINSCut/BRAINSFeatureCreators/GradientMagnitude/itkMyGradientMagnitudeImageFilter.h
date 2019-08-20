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
#ifndef __itkMyGradientMagnitudeImageFilter_h
#define __itkMyGradientMagnitudeImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class GradientMagnitudeImageFilter
 * \brief Computes the gradient magnitude of an image region at each pixel.
 *
 * \ingroup GradientFilters
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 * \ingroup ITKImageGradient
 *
 * \wiki
 * \wikiexample{EdgesAndGradients/GradientMagnitudeImageFilter,Compute the gradient magnitude image}
 * \endwiki
 */
template < typename TInputImage, typename TOutputImage >
class GradientMagnitudeImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( GradientMagnitudeImageFilter );

  /** Standard class type alias. */
  using Self = GradientMagnitudeImageFilter;
  using Superclass = ImageToImageFilter< TInputImage, TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( GradientMagnitudeImageFilter, ImageToImageFilter );

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  using OutputPixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;
  using RealType = typename NumericTraits< InputPixelType >::RealType;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** Image type alias support */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass type alias. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  /** GradientMagnitudeImageFilter needs a larger input requested
   * region than the output requested region (larger by the kernel
   * size to calculate derivatives).  As such,
   * GradientMagnitudeImageFilter needs to provide an implementation
   * for GenerateInputRequestedRegion() in order to inform the
   * pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void
  GenerateInputRequestedRegion() throw( InvalidRequestedRegionError );

  /** Use the image spacing information in calculations. Use this option if you
   *  want derivatives in physical space. Default is UseImageSpacingOn. */
  void
  SetUseImageSpacingOn()
  {
    this->SetUseImageSpacing( true );
  }

  /** Ignore the image spacing. Use this option if you want derivatives in
      isotropic pixel space.  Default is UseImageSpacingOn. */
  void
  SetUseImageSpacingOff()
  {
    this->SetUseImageSpacing( false );
  }

  /** Set/Get whether or not the filter will use the spacing of the input
      image in its calculations */
  itkSetMacro( UseImageSpacing, bool );
  itkGetConstMacro( UseImageSpacing, bool );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputHasNumericTraitsCheck, (Concept::HasNumericTraits< InputPixelType >));
  /** End concept checking */
#endif
protected:
  GradientMagnitudeImageFilter() { m_UseImageSpacing = true; }

  virtual ~GradientMagnitudeImageFilter() {}

  /** GradientMagnitudeImageFilter can be implemented as a
   * multithreaded filter.  Therefore, this implementation provides a
   * ThreadedGenerateData() routine which is called for each
   * processing thread. The output image data is allocated
   * automatically by the superclass prior to calling
   * ThreadedGenerateData().  ThreadedGenerateData can only write to
   * the portion of the output image specified by the parameter
   * "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void
  ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId );

  void
  PrintSelf( std::ostream &, Indent ) const;

private:
  bool m_UseImageSpacing;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMyGradientMagnitudeImageFilter.hxx"
#endif

#endif
