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
template <typename TInputImage, typename TOutputImage>
class BRAINSHoughEyeDetector : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard ITK type alias */
  using Self = BRAINSHoughEyeDetector;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

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
  using InputImageConstIterator = ImageRegionConstIterator<TInputImage>;

  /** Output image type alias */
  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelType = typename TOutputImage::PixelType;
  using OutputImageRegionType = typename TOutputImage::RegionType;
  using OutputPointType = typename TOutputImage::PointType;

  using OutputImageIterator = ImageRegionIterator<TOutputImage>;
  using OutputImageIteratorWithIndex = ImageRegionIteratorWithIndex<TOutputImage>;

  /* Transform and filter type alias */
  using VersorTransformType = VersorRigid3DTransform<double>;
  using VersorVectorType = typename VersorTransformType::OutputVectorType;

  using WriterType = ImageFileWriter<TInputImage>;
  using HoughFilterType = HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage>;
  using SpheresListType = typename HoughFilterType::SpheresListType;
  using SphereIterator = typename SpheresListType::const_iterator;
  using HoughFilterPointer = typename HoughFilterType::Pointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(BRAINSHoughEyeDetector, ImageToImageFilter);

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Display */
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Set/Get the number of circles to extract */
  itkSetMacro(NumberOfSpheres, unsigned int);

  /** Set the minimum radiu value the filter should look for */
  itkSetMacro(MinimumRadius, double);

  /** Set the maximum radius value the filter should look for */
  itkSetMacro(MaximumRadius, double);

  /** Set the scale of the derivative function (using DoG) */
  itkSetMacro(SigmaGradient, double);

  /** Set the variance of the gaussian bluring for the accumulator */
  itkSetMacro(Variance, double);

  /** Set the radius of the disc to remove from the accumulator
   *  for each circle found */
  itkSetMacro(SphereRadiusRatio, double);

  /** Set the voting radius */
  itkSetMacro(VotingRadiusRatio, double);

  /** Set the threshold above which the filter should consider
   the point as a valid point */
  itkSetMacro(Threshold, double);

  /** Set the threshold above which the filter should consider
   the point as a valid point */
  itkSetMacro(OutputThreshold, double);

  /** Set the threshold above which the filter should consider
   the point as a valid point */
  itkSetMacro(GradientThreshold, double);

  /** Set the number of threads */
  itkSetMacro(NbOfThreads, unsigned int);

  /** Set the number of threads */
  itkSetMacro(SamplingRatio, double);

  /** Set the mode of the algorithm
   * HoughEyeDetectorMode = 0: Detecting bright spheres in a dark environment
   * HoughEyeDetectorMode = 1: Detecting dark spheres in a bright environment */
  itkSetMacro(HoughEyeDetectorMode, int);

  /** Set the center of head mass of the image */
  itkSetMacro(orig_lmk_CenterOfHeadMass, InputPointType);

  /** Get the left eye center coordinate */
  itkGetMacro(orig_lmk_LE, InputPointType);

  /** Get the right eye center coordinate */
  itkGetMacro(orig_lmk_RE, InputPointType);

  /** Set the debug output dir */
  itkSetMacro(ResultsDir, std::string);

  /** Set the write debug image level */
  itkSetMacro(WritedebuggingImagesLevel, unsigned int);

  /** Get the accumulator image */
  itkGetConstObjectMacro(AccumulatorImage, TInputImage);

  /** Get the RoI image */
  itkGetConstObjectMacro(RoIImage, TInputImage);

  /** Get the maximum output pixel value */
  itkGetConstMacro(MaxInputPixelValue, OutputPixelType);

  /** Get the minimum output pixel value */
  itkGetConstMacro(MinInputPixelValue, OutputPixelType);

  /** Get the versor transform of the detector */
  itkGetModifiableObjectMacro(orig2eyeFixedTransform, VersorTransformType);

  /** Get/Set the failure report */
  itkGetConstMacro(Failure, bool);
  itkSetMacro(Failure, bool);


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(IntConvertibleToOutputCheck, (Concept::Convertible<int, OutputPixelType>));
  itkConceptMacro(InputGreaterThanDoubleCheck, (Concept::GreaterThanComparable<InputPixelType, double>));
  itkConceptMacro(OutputPlusIntCheck, (Concept::AdditiveOperators<OutputPixelType, int>));
  itkConceptMacro(OutputDividedByIntCheck, (Concept::DivisionOperators<OutputPixelType, int>));
  /** End concept checking */
#endif
  BRAINSHoughEyeDetector(const Self &) = delete;
  void
  operator=(const Self &) = delete;
  ~BRAINSHoughEyeDetector() override = default;

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
  InputPointType m_orig_lmk_CenterOfHeadMass;

  // Debug settings
  std::string  m_ResultsDir;
  unsigned int m_WritedebuggingImagesLevel;

  /** Output parameters */
  OutputImagePointer m_AccumulatorImage;
  OutputImagePointer m_RoIImage;
  OutputPointType    m_orig_lmk_LE;
  OutputPointType    m_orig_lmk_RE;
  bool               m_Failure; // indicating whether the detector realizes the failure
  OutputPixelType    m_MaxInputPixelValue;
  OutputPixelType    m_MinInputPixelValue;

  VersorTransformType::Pointer m_orig2eyeFixedTransform;

  // Eye Radius rages from 11-13 mm
  static constexpr double default_minimum_radius = 11.0; // mm
  static constexpr double default_maximum_radius = 13.0; // mm
  // After evaluating over 9000 T1 weighted space scans from the PREDICT HD project
  // 4885 scans across 50 scanners evaluated the radius enclosing the eye
  // center was 68 and max was (101.7,  64.98)
  // mean=81.467760 std=4.448509

  // Interior radius (mm), exterior radius (mm), of the shell-like RoI
  static constexpr double m_R1{ 64.98 - 1.25 * default_maximum_radius }; // was 30
  static constexpr double m_R2{ 101.7 + 1.25 * default_maximum_radius }; // was 120

private:
  OutputImagePointer
  MakeROICandiadteRegion(InputImageConstPointer image, float min_si = -52., float max_si = +52.);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "BRAINSHoughEyeDetector.hxx"
#endif

#endif
