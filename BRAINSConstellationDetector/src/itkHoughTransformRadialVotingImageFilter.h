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
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 585 $  // Revision of last commit
  Date: $Date: 2009-08-20 21:25:19 -0400 (Thu, 20 Aug 2009) $  // Date of last commit
=========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#ifndef __itkHoughTransformRadialVotingImageFilter_h
#define __itkHoughTransformRadialVotingImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkEllipseSpatialObject.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCastImageFilter.h"
#include <itkGaussianDistribution.h>
#include "itkAddImageFilter.h"
#include "itkImageRegionIterator.h"

namespace itk
{
/**
 * \class HoughTransformRadialVotingImageFilter
 * \brief Performs the Hough Transform to find circles in a 2D image.
 *
 * This filter derives from the base class ImageToImageFilter
 * The input is an image, and all pixels above some threshold are those
 * we want to consider during the process.
 *
 * This filter produces two output:
 *   1) The accumulator array, which represents probability of centers.
 *   2) The array or radii, which has the radius value at each coordinate point.
 *
 *  When the filter found a "correct" point, it computes the gradient at this
 * point and votes on a small region defined using the minimum and maximum
 * radius given by the user, and fill in the array of radii.
 *
 * \ingroup ImageFeatureExtraction
 * \todo Update the doxygen documentation!!!
 * */

template < typename TInputImage, typename TOutputImage >
class HoughTransformRadialVotingImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard "Self" type alias. */
  using Self = HoughTransformRadialVotingImageFilter;
  /** Standard "Superclass" type alias. */
  using Superclass = ImageToImageFilter< TInputImage, TOutputImage >;
  /** Smart pointer type alias support. */
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Input Image type alias */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using InputIndexType = typename InputImageType::IndexType;
  using InputIndexValueType = typename InputIndexType::IndexValueType;
  using InputPixelType = typename InputImageType::PixelType;
  using InputSizeType = typename InputImageType::SizeType;
  using InputSizeValueType = typename InputSizeType::SizeValueType;
  using InputRegionType = typename InputImageType::RegionType;
  using InputSpacingType = typename InputImageType::SpacingType;
  using InputPointType = typename InputImageType::PointType;
  using InputCoordType = typename InputPointType::CoordRepType;

  /** Output Image type alias */
  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelType = typename OutputImageType::PixelType;
  using OutputImageRegionType = typename OutputImageType::RegionType;

  using InternalImageType = Image< InputCoordType, ImageDimension >;
  using InternalImagePointer = typename InternalImageType::Pointer;
  using InternalIndexType = typename InternalImageType::IndexType;
  using InternalIndexValueType = typename InternalIndexType::IndexValueType;
  using InternalPixelType = typename InternalImageType::PixelType;
  using InternalRegionType = typename InternalImageType::RegionType;
  using InternalSizeType = typename InternalImageType::SizeType;
  using InternalSizeValueType = typename InternalSizeType::SizeValueType;
  using InternalSpacingType = typename InternalImageType::SpacingType;

  /** Sphere type alias */
  using SphereType = EllipseSpatialObject< ImageDimension >;
  using SpherePointer = typename SphereType::Pointer;
  using SphereVectorType = typename SphereType::VectorType;
  using SpheresListType = std::list< SpherePointer >;

  using InternalIteratorType = ImageRegionIterator< InternalImageType >;
  using OutputIteratorType = ImageRegionIterator< OutputImageType >;

  using GaussianFilterType = DiscreteGaussianImageFilter< OutputImageType, InternalImageType >;
  typedef typename GaussianFilterType::Pointer GaussianFilterPointer;

  using GaussianFunctionType = itk::Statistics::GaussianDistribution;
  using GaussianFunctionPointer = typename GaussianFunctionType::Pointer;

  using DoGFunctionType = GaussianDerivativeImageFunction< InputImageType, InputCoordType >;
  typedef typename DoGFunctionType::Pointer    DoGFunctionPointer;
  typedef typename DoGFunctionType::VectorType DoGVectorType;

  using MinMaxCalculatorType = MinimumMaximumImageCalculator< InternalImageType >;
  using MinMaxCalculatorPointer = typename MinMaxCalculatorType::Pointer;

  using CastFilterType = CastImageFilter< InternalImageType, OutputImageType >;
  using CastFilterPointer = typename CastFilterType::Pointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( HoughTransformRadialVotingImageFilter, ImageToImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Set both Minimum and Maximum radius values */
  void
  SetRadius( InputCoordType radius );

  /** Set the minimum radiu value the filter should look for */
  itkSetMacro( MinimumRadius, InputCoordType );

  /** Set the maximum radius value the filter should look for */
  itkSetMacro( MaximumRadius, InputCoordType );

  /** Set the threshold above which the filter should consider
   the point as a valid point */
  itkSetMacro( Threshold, double );

  /** Get the threshold value */
  itkGetConstMacro( Threshold, double );

  /** Set the threshold above which the filter should consider
      the point as a valid point */
  itkSetMacro( OutputThreshold, InternalPixelType );

  /** Get the threshold value */
  itkGetConstMacro( OutputThreshold, InternalPixelType );

  /** Set the threshold above which the filter should consider
      the point as a valid point */
  itkSetMacro( GradientThreshold, InputCoordType );

  /** Get the threshold value */
  itkGetConstMacro( GradientThreshold, InputCoordType );

  /** Get the radius image */
  itkGetConstObjectMacro( RadiusImage, InternalImageType );

  /** Get the accumulator image */
  itkGetConstObjectMacro( AccumulatorImage, InternalImageType );

  /** Set the scale of the derivative function (using DoG) */
  itkSetMacro( SigmaGradient, double );

  /** Get the scale value */
  itkGetConstMacro( SigmaGradient, double );

  /** Get the list of circles. This recomputes the circles */
  SpheresListType &
  GetSpheres();

  /** Set/Get the number of circles to extract */
  itkSetMacro( NumberOfSpheres, unsigned int );
  itkGetConstMacro( NumberOfSpheres, unsigned int );

  /** Set/Get the radius of the disc to remove from the accumulator
   *  for each circle found */
  itkSetMacro( SphereRadiusRatio, InputCoordType );
  itkGetConstMacro( SphereRadiusRatio, InputCoordType );

  itkSetMacro( VotingRadiusRatio, InputCoordType );
  itkGetConstMacro( VotingRadiusRatio, InputCoordType );

  /** Set the variance of the gaussian bluring for the accumulator */
  itkSetMacro( Variance, double );
  itkGetConstMacro( Variance, double );

  /** Set the number of threads */
  itkSetMacro( NbOfThreads, unsigned int );
  itkGetConstMacro( NbOfThreads, unsigned int );

  /** Set the number of threads */
  itkSetMacro( SamplingRatio, double );
  itkGetConstMacro( SamplingRatio, double );

  /** Set the mode of the algorithm */
  /** HoughEyeDetectorMode = 0: Detecting bright spheres in a dark environment.
   */
  /** HoughEyeDetectorMode = 1: Detecting dark spheres in a bright environment.
   */
  itkSetMacro( HoughEyeDetectorMode, int );

  /** Get the mode of the algorithm */
  itkGetConstMacro( HoughEyeDetectorMode, int );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( IntConvertibleToOutputCheck, (Concept::Convertible< int, OutputPixelType >));
  itkConceptMacro( InputGreaterThanDoubleCheck, (Concept::GreaterThanComparable< InputPixelType, double >));
  itkConceptMacro( OutputPlusIntCheck, (Concept::AdditiveOperators< OutputPixelType, int >));
  itkConceptMacro( OutputDividedByIntCheck, (Concept::DivisionOperators< OutputPixelType, int >));
  /** End concept checking */
#endif
protected:
  HoughTransformRadialVotingImageFilter();
  ~HoughTransformRadialVotingImageFilter() override;

  InputCoordType       m_MinimumRadius;
  InputCoordType       m_MaximumRadius;
  double               m_Threshold;
  InputCoordType       m_GradientThreshold;
  InternalPixelType    m_OutputThreshold;
  double               m_SigmaGradient;
  double               m_Variance;
  InputCoordType       m_VotingRadiusRatio;
  InputCoordType       m_SphereRadiusRatio;
  double               m_SamplingRatio;
  InternalImagePointer m_RadiusImage;
  InternalImagePointer m_AccumulatorImage;
  SpheresListType      m_SpheresList;
  unsigned int         m_NumberOfSpheres;
  unsigned long        m_OldModifiedTime;
  unsigned int         m_NbOfThreads;
  bool                 m_AllSeedsProcessed;

  // -- Add by Wei Lu
  int m_HoughEyeDetectorMode;

  /** Method for evaluating the implicit function over the image. */
  void
  BeforeThreadedGenerateData() override;

  void
  AfterThreadedGenerateData() override;

  void
  ThreadedGenerateData( const OutputImageRegionType & windowRegion, ThreadIdType threadId ) override;

  void
  PrintSelf( std::ostream & os, Indent indent ) const override;

  /** HoughTransformRadialVotingImageFilter needs the entire input. Therefore
   * it must provide an implementation GenerateInputRequestedRegion().
   * \sa ProcessObject::GenerateInputRequestedRegion(). */
  void
  GenerateInputRequestedRegion() override;

  /** HoughTransformRadialVotingImageFilter's produces all the output.
   * Therefore, it must provide an implementation of
   * EnlargeOutputRequestedRegion.
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  void
  EnlargeOutputRequestedRegion( DataObject * itkNotUsed( output ) ) override;

  void
  ComputeMeanRadiusImage();

private:
  HoughTransformRadialVotingImageFilter( const Self & ) {}

  void
  operator=( const Self & )
  {}
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkHoughTransformRadialVotingImageFilter.hxx"
#endif

#endif
