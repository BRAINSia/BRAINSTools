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

#include "itkMultiThreader.h"

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
template <class TInputImage, class TOutputImage>
class BRAINSHoughEyeDetector :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  /** Standard ITK typedef */
  typedef BRAINSHoughEyeDetector                        Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  itkStaticConstMacro(Dimension, unsigned int, TInputImage::ImageDimension);

  /** Input image typedef */
  typedef typename TInputImage::Pointer         InputImagePointer;
  typedef typename TInputImage::ConstPointer    InputImageConstPointer;
  typedef typename TInputImage::IndexType       InputIndexType;
  typedef typename TInputImage::PixelType       InputPixelType;
  typedef typename TInputImage::SizeType        InputSizeType;
  typedef typename TInputImage::RegionType      InputRegionType;
  typedef typename TInputImage::SpacingType     InputSpacingType;
  typedef typename TInputImage::PointType       InputPointType;
  typedef typename TInputImage::DirectionType   InputDirectionType;
  typedef ImageRegionConstIterator<TInputImage> InputImageConstIterator;

  /** Output image typedef */
  typedef typename TOutputImage::Pointer    OutputImagePointer;
  typedef typename TOutputImage::PixelType  OutputPixelType;
  typedef typename TOutputImage::RegionType OutputImageRegionType;
  typedef typename TOutputImage::PointType  OutputPointType;

  typedef ImageRegionIterator<TOutputImage>          OutputImageIterator;
  typedef ImageRegionIteratorWithIndex<TOutputImage> OutputImageIteratorWithIndex;

  /* Transform and filter typedef */
  typedef VersorRigid3DTransform<double>                 VersorTransformType;
  typedef typename VersorTransformType::OutputVectorType VersorVectorType;

  typedef ImageFileWriter<TInputImage>                                     WriterType;
  typedef HoughTransformRadialVotingImageFilter<TInputImage, TOutputImage> HoughFilterType;
  typedef typename HoughFilterType::SpheresListType                        SpheresListType;
  typedef typename SpheresListType::const_iterator                         SphereIterator;
  typedef typename HoughFilterType::Pointer                                HoughFilterPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(BRAINSHoughEyeDetector, ImageToImageFilter);

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Display */
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

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
  itkSetMacro(CenterOfHeadMass, InputPointType);

  /** Set the interior radius of the shell-like RoI */
  itkSetMacro(R1, double);

  /** Set the exterior radius of the shell-like RoI */
  itkSetMacro(R2, double);

  /** Set the spread angle of the shell-like RoI */
  itkSetMacro(Theta, double);

  /** Get the left eye center coordinate */
  itkGetMacro(LE, InputPointType);

  /** Get the right eye center coordinate */
  itkGetMacro(RE, InputPointType);

  /** Set the debug output dir */
  itkSetMacro(ResultsDir, std::string);

  /** Set the write debug image level */
  itkSetMacro(WritedebuggingImagesLevel, unsigned int);

  /** Get the accumulator image */
  itkGetConstObjectMacro(AccumulatorImage, TInputImage);

  /** Get the RoI image */
  itkGetConstObjectMacro(RoIImage, TInputImage);

  /** Get the rotation angle of the alignment process */
  itkGetMacro(RotAngle, InputPointType);

  /** Get the adult interpupilary distance */
  itkGetConstMacro(Ipd, double);

  /** Get the maximum output pixel value */
  itkGetConstMacro(MaxInputPixelValue, OutputPixelType);

  /** Get the minimum output pixel value */
  itkGetConstMacro(MinInputPixelValue, OutputPixelType);

  /** Get the versor transform of the detector */
  itkGetModifiableObjectMacro(VersorTransform, VersorTransformType);
  itkGetModifiableObjectMacro(InvVersorTransform, VersorTransformType);

  /** Get/Set the failure report */
  itkGetConstMacro(Failure, bool);
  itkSetMacro(Failure, bool);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( IntConvertibleToOutputCheck,
                   ( Concept::Convertible<int, OutputPixelType> ) );
  itkConceptMacro( InputGreaterThanDoubleCheck,
                   ( Concept::GreaterThanComparable<InputPixelType, double> ) );
  itkConceptMacro( OutputPlusIntCheck,
                   ( Concept::AdditiveOperators<OutputPixelType, int> ) );
  itkConceptMacro( OutputDividedByIntCheck,
                   ( Concept::DivisionOperators<OutputPixelType, int> ) );
  /** End concept checking */
#endif
protected:

  BRAINSHoughEyeDetector();

  void GenerateData() ITK_OVERRIDE;

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
  OutputPointType    m_RotAngle;
  double             m_Ipd;     // adult interpupilary distance
  bool               m_Failure; // indicating wether the detector realizes the
                                // failure

  OutputPixelType m_MaxInputPixelValue;
  OutputPixelType m_MinInputPixelValue;

  VersorTransformType::Pointer m_VersorTransform;
  VersorTransformType::Pointer m_InvVersorTransform;
private:

  BRAINSHoughEyeDetector(const Self &)
  {
  }

  void operator=(const Self &)
  {
  }
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "BRAINSHoughEyeDetector.hxx"
#endif

#endif
