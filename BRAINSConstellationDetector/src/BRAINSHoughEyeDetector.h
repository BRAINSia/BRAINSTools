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
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef __BRAINSHoughEyeDetector_h
#define __BRAINSHoughEyeDetector_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkVersorRigid3DTransform.h"

#include "itkResampleInPlaceImageFilter.h"
#include "itkHoughTransformRadialVotingImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"

#include <cstring>

namespace itk
{
template < typename TInputImage, typename TOutputImage >
class BRAINSHoughEyeDetector : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard ITK type alias */
  using Self = BRAINSHoughEyeDetector;
  using Superclass = ImageToImageFilter< TInputImage, TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  static constexpr unsigned int Dimension = TInputImage::ImageDimension;

  /** Input image type alias */
  using InputImagePointer = typename TInputImage::Pointer;
  using InputImageConstPointer = typename TInputImage::ConstPointer;
  using InputIndexType = typename TInputImage::IndexType;
  using InputPixelType = typename TInputImage::PixelType;
  using InputSizeType = typename TInputImage::SizeType;
  using InputRegionType = typename TInputImage::RegionType;
  using InputSpacingType = typename TInputImage::SpacingType;
  using InputPointType = typename TInputImage::PointType;
  using InputDirectionType = typename TInputImage::DirectionType;
  using InputImageConstIterator = ImageRegionConstIterator< TInputImage >;

  /** Output image type alias */
  using OutputImagePointer = typename TOutputImage::Pointer;
  using OutputPixelType = typename TOutputImage::PixelType;
  using OutputImageRegionType = typename TOutputImage::RegionType;
  using OutputPointType = typename TOutputImage::PointType;

  using OutputImageIterator = ImageRegionIterator< TOutputImage >;
  using OutputImageIteratorWithIndex = ImageRegionIteratorWithIndex< TOutputImage >;

  /* Transform and filter type alias */
  using VersorTransformType = VersorRigid3DTransform< double >;
  using VersorVectorType = typename VersorTransformType::OutputVectorType;

  using WriterType = ImageFileWriter< TInputImage >;
  using HoughFilterType = HoughTransformRadialVotingImageFilter< TInputImage, TOutputImage >;
  using SpheresListType = typename HoughFilterType::SpheresListType;
  using SphereIterator = typename SpheresListType::const_iterator;
  using HoughFilterPointer = typename HoughFilterType::Pointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro( BRAINSHoughEyeDetector, ImageToImageFilter );

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Display */
  void
  PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Set/Get the number of circles to extract */
  itkSetMacro( NumberOfSpheres, unsigned int );

  /** Set the minimum radiu value the filter should look for */
  itkSetMacro( MinimumRadius, double );

  /** Set the maximum radius value the filter should look for */
  itkSetMacro( MaximumRadius, double );

  /** Set the scale of the derivative function (using DoG) */
  itkSetMacro( SigmaGradient, double );

  /** Set the variance of the gaussian bluring for the accumulator */
  itkSetMacro( Variance, double );

  /** Set the radius of the disc to remove from the accumulator
   *  for each circle found */
  itkSetMacro( SphereRadiusRatio, double );

  /** Set the voting radius */
  itkSetMacro( VotingRadiusRatio, double );

  /** Set the threshold above which the filter should consider
   the point as a valid point */
  itkSetMacro( Threshold, double );

  /** Set the threshold above which the filter should consider
   the point as a valid point */
  itkSetMacro( OutputThreshold, double );

  /** Set the threshold above which the filter should consider
   the point as a valid point */
  itkSetMacro( GradientThreshold, double );

  /** Set the number of threads */
  itkSetMacro( NbOfThreads, unsigned int );

  /** Set the number of threads */
  itkSetMacro( SamplingRatio, double );

  /** Set the mode of the algorithm
   * HoughEyeDetectorMode = 0: Detecting bright spheres in a dark environment
   * HoughEyeDetectorMode = 1: Detecting dark spheres in a bright environment */
  itkSetMacro( HoughEyeDetectorMode, int );

  /** Set the center of head mass of the image */
  itkSetMacro( CenterOfHeadMass, InputPointType );

  /** Set the interior radius of the shell-like RoI */
  itkSetMacro( R1, double );

  /** Set the exterior radius of the shell-like RoI */
  itkSetMacro( R2, double );

  /** Set the spread angle of the shell-like RoI */
  itkSetMacro( Theta, double );

  /** Get the left eye center coordinate */
  itkGetMacro( LE, InputPointType );

  /** Get the right eye center coordinate */
  itkGetMacro( RE, InputPointType );

  /** Set the debug output dir */
  itkSetMacro( ResultsDir, std::string );

  /** Set the write debug image level */
  itkSetMacro( WritedebuggingImagesLevel, unsigned int );

  /** Get the accumulator image */
  itkGetConstObjectMacro( AccumulatorImage, TInputImage );

  /** Get the RoI image */
  itkGetConstObjectMacro( RoIImage, TInputImage );

  /** Get the maximum output pixel value */
  itkGetConstMacro( MaxInputPixelValue, OutputPixelType );

  /** Get the minimum output pixel value */
  itkGetConstMacro( MinInputPixelValue, OutputPixelType );

  /** Get the versor transform of the detector */
  itkGetModifiableObjectMacro( VersorTransform, VersorTransformType );
  itkGetModifiableObjectMacro( InvVersorTransform, VersorTransformType );

  /** Get/Set the failure report */
  itkGetConstMacro( Failure, bool );
  itkSetMacro( Failure, bool );


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( IntConvertibleToOutputCheck, (Concept::Convertible< int, OutputPixelType >));
  itkConceptMacro( InputGreaterThanDoubleCheck, (Concept::GreaterThanComparable< InputPixelType, double >));
  itkConceptMacro( OutputPlusIntCheck, (Concept::AdditiveOperators< OutputPixelType, int >));
  itkConceptMacro( OutputDividedByIntCheck, (Concept::DivisionOperators< OutputPixelType, int >));
  /** End concept checking */
#endif
  BRAINSHoughEyeDetector( const Self & ) = delete;
  void
  operator=( const Self & ) = delete;
  virtual ~BRAINSHoughEyeDetector() = default;

protected:
  BRAINSHoughEyeDetector();
  void
  GenerateData() override;

  /** Input Parameters */
  // Pass parameters from Hough Transform Radial Voting Filter
  unsigned int   m_NumberOfSpheres;
  double         m_MinimumRadius;
  double         m_MaximumRadius;
  double         m_SigmaGradient;
  double         m_Variance;
  double         m_SphereRadiusRatio;
  double         m_VotingRadiusRatio;
  double         m_Threshold;
  double         m_OutputThreshold;
  double         m_GradientThreshold;
  unsigned int   m_NbOfThreads;
  double         m_SamplingRatio;
  int            m_HoughEyeDetectorMode;
  InputPointType m_CenterOfHeadMass;

  // Interior radius (mm), exterior radius (mm), and spread
  // angle (rad) of the shell-like RoI.
  double             m_R1;
  double             m_R2;
  double             m_Theta;
  OutputImagePointer m_OutputImage;

  // Debug settings
  std::string  m_ResultsDir;
  unsigned int m_WritedebuggingImagesLevel;

  /** Output parameters */
  OutputImagePointer m_AccumulatorImage;
  OutputImagePointer m_RoIImage;
  OutputPointType    m_LE;
  OutputPointType    m_RE;
  bool               m_Failure; // indicating whether the detector realizes the failure
  OutputPixelType    m_MaxInputPixelValue;
  OutputPixelType    m_MinInputPixelValue;

  VersorTransformType::Pointer m_VersorTransform;
  VersorTransformType::Pointer m_InvVersorTransform;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "BRAINSHoughEyeDetector.hxx"
#endif

#endif
