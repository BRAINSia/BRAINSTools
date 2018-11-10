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
 *  itkICCDeformableRegistrationFilter.h
 *  ICCDeformationTools
 *
 *  Created by Yongqiang Zhao on 4/21/09.
 *  Copyright 2009 The University of Iowa. All rights reserved.
 *
 */

#ifndef __itkICCDeformableRegistrationFilter_h
#define __itkICCDeformableRegistrationFilter_h

#include "itkPDEDeformableRegistrationFilter.h"
#include "itkICCDeformableFunction.h"

#include "itkMultiplyImageFilter.h"
#include "itkExponentialDisplacementFieldImageFilter.h"
#include "itkWarpVectorImageFilter.h"
#include "itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h"
#include "itkAddImageFilter.h"
#include <vector>
#include <complex>

#include "itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter.h"
#include "itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter.h"

#include "itkSpatialObject.h"
#include "itkWarpImageFilter.h"
#include "itkLandmarkSpatialObject.h"
#include "itkPointSet.h"

namespace itk
{
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
class ICCDeformableRegistrationFilter :
  public         PDEDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                                 TDisplacementField>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(ICCDeformableRegistrationFilter);

  /** Standard class type alias. */
  using Self = ICCDeformableRegistrationFilter;
  using Superclass = PDEDeformableRegistrationFilter<
      TFixedImage, TMovingImage, TDisplacementField>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ICCDeformableRegistrationFilter,
                PDEDeformableRegistrationFilter );

  /** FixedImage image type. */
  using FixedImageType = typename Superclass::FixedImageType;
  using FixedImagePointer = typename Superclass::FixedImagePointer;

  /** MovingImage image type. */
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingImagePointer = typename Superclass::MovingImagePointer;

  /** Deformation field type. */
  using DisplacementFieldType = typename Superclass::DisplacementFieldType;
  using DisplacementFieldPointer = typename Superclass::DisplacementFieldPointer;

  using OutputImageType = DisplacementFieldType;

  /** Inherit some enums from the superclass. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;
  using MaskPixelType = unsigned char;
  using MaskImageType = Image<MaskPixelType, Self::ImageDimension>;
  using MaskImagePointer = typename MaskImageType::Pointer;

  using MaskType = SpatialObject<Self::ImageDimension>;
  using MaskPointer = typename MaskType::Pointer;

  using DataObjectPointer = typename DataObject::Pointer;

  using LandmarkType = LandmarkSpatialObject<Self::ImageDimension>;
  using LandmarkPointer = typename LandmarkType::Pointer;

  using PointSetType = PointSet<typename MovingImageType::PixelType, Self::ImageDimension>;
  using PointSetPointer = typename PointSetType::Pointer;

  using ComplexPixelType = std::complex<float>;
  using DisplacementFieldFFTPixelType = Vector<ComplexPixelType, 3>;
  using ComplexImageType = Image<ComplexPixelType, 3>;
  using DisplacementFieldFFTType = Image<DisplacementFieldFFTPixelType, 3>;

  using FFTWComplexToRealImageType = VectorFFTWHalfHermitianToRealInverseFFTImageFilter<typename TDisplacementField::PixelType,
                                                             3>;
  using FFTWRealToComplexImageType = VectorFFTWRealToHalfHermitianForwardFFTImageFilter<typename TDisplacementField::PixelType,
                                                             3>;
  typedef typename FFTWComplexToRealImageType::Pointer
    FFTWComplexToRealImagePointer;
  typedef typename FFTWRealToComplexImageType::Pointer
    FFTWRealToComplexImagePointer;

  using DisplacementFieldFFTPointer = typename DisplacementFieldFFTType::Pointer;
  using ComplexImagePointer = typename ComplexImageType::Pointer;

  /** FiniteDifferenceFunction type. */
  typedef typename
    Superclass::FiniteDifferenceFunctionType               FiniteDifferenceFunctionType;

  /** Take timestep type from the FiniteDifferenceFunction. */
  typedef typename
    FiniteDifferenceFunctionType::TimeStepType             TimeStepType;

  using MaskWarperType = WarpImageFilter<MaskImageType,
                          MaskImageType, DisplacementFieldType>;
  using MaskWarperPointer = typename MaskWarperType::Pointer;

  /** DemonsRegistrationFilterFunction type. */
  using ICCDeformableFunctionType = ICCDeformableFunction<
      FixedImageType,
      MovingImageType, DisplacementFieldType>;
  typedef typename
    ICCDeformableFunctionType::GradientType        GradientType;

  /** Get the metric value. The metric value is the mean square difference
   * in intensity between the fixed image and transforming moving image
   * computed over the the overlapping region between the two images.
   * This value is calculated for the current iteration */
  virtual double GetForwardMetric() const;

  virtual double GetBackwardMetric() const;

  virtual const double & GetForwardRMSChange() const;

  virtual const double & GetBackwardRMSChange() const;

  itkSetMacro(BackgroundFilledValue, float);
  itkGetMacro(BackgroundFilledValue, float);
//  virtual void SetUseGradientType( GradientType gtype );
//  virtual GradientType GetUseGradientType() const;

  void CopyInputToOutput() override;

  void PostProcessOutput() override { }

  /** Use a first-order approximation of the exponential.
   *  This amounts to using an update rule of the type
   *  s <- s o (Id + u) instead of s <- s o exp(u) */
  itkSetMacro( UseFirstOrderExp, bool );
  itkGetMacro( UseFirstOrderExp, bool );
  itkBooleanMacro( UseFirstOrderExp );

  itkSetMacro( UseConsistentLandmark, bool);
  itkGetMacro( UseConsistentLandmark, bool);
  itkBooleanMacro( UseConsistentLandmark);

  itkSetMacro( UseConsistentIntensity, bool);
  itkGetMacro( UseConsistentIntensity, bool);
  itkBooleanMacro(UseConsistentIntensity);
  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void SetIntensityDifferenceThreshold(double);

  virtual double GetIntensityDifferenceThreshold() const;

  /** Set/Get the maximum length in terms of pixels of
   *  the vectors in the update buffer. */
  virtual void SetMaximumUpdateStepLength(double);

  virtual double GetMaximumUpdateStepLength() const;

  virtual std::vector<DisplacementFieldPointer> & GetUpdateBuffers()
  {
    return m_UpdateBuffers;
  }

//  virtual std::vector<DisplacementFieldPointer> & GetDisplacementFields()
//    { return m_DisplacementFields;}

  DisplacementFieldType * GetForwardDisplacementField()
  {
    return this->GetOutput(0);
  }

  DisplacementFieldType * GetBackwardDisplacementField()
  {
    return this->GetOutput(1);
  }

  itkSetMacro( RegularizationWeight, float );
  itkGetMacro( RegularizationWeight, float );

  itkSetMacro( InverseWeight, float );
  itkGetMacro( InverseWeight, float );

  itkSetMacro( SimilarityWeight, float );
  itkGetMacro( SimilarityWeight, float );

  itkSetMacro( LandmarkWeight, float );
  itkGetMacro( LandmarkWeight, float );

  itkSetMacro( HarmonicPercent, float );
  itkGetMacro( HarmonicPercent, float );

  itkSetMacro( MinJac, float );
  itkGetMacro( MinJac, float );

  using Superclass::MakeOutput;
  ProcessObject::DataObjectPointer MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx) override;

  virtual void SetInitialForwardDisplacementField( DisplacementFieldType * ptr )
  {
    this->SetInput( 3, ptr );
  }

  virtual void SetInitialBackwardDisplacementField( DisplacementFieldType * ptr )
  {
    this->SetInput( 4, ptr );
  }

  void GenerateInputRequestedRegion() override;

  virtual void SetMovingImageMask( MaskType *mask);

  virtual const MaskType * GetMovingImageMask() const;

  virtual void SetFixedImageMask(MaskType *mask);

  virtual const MaskType * GetFixedImageMask() const;

  virtual void SetFixedLandmark(PointSetType * landmark);

  virtual void SetMovingLandmark(PointSetType * landmark);

protected:
  ICCDeformableRegistrationFilter();
  ~ICCDeformableRegistrationFilter() override
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const override;

  /** Initialize the state of filter and equation before each iteration. */
  void InitializeIteration() override;

  /** This method allocates storage in m_UpdateBuffer.  It is called from
   * FiniteDifferenceFilter::GenerateData(). */
  void AllocateUpdateBuffer() override;

  /** Apply update. */
  void ApplyUpdate(const TimeStepType& dt) override;

  TimeStepType CalculateChange() override;

  void Initialize() override;

  void ComputeLinearElastic(DisplacementFieldFFTPointer &, float normalizer);

  void ComputeInverseConsistency(DisplacementFieldFFTPointer &, DisplacementFieldFFTPointer &, float normalizer);

  double ComputeMinJac(DisplacementFieldPointer & );

  /** This method returns a pointer to a FiniteDifferenceFunction object that
 * will be used by the filter to calculate updates at image pixels.
 * \returns A FiniteDifferenceObject pointer. */
  itkGetConstReferenceObjectMacro(BackwardDifferenceFunction,
                                  FiniteDifferenceFunctionType );

  /** This method sets the pointer to a FiniteDifferenceFunction object that
   * will be used by the filter to calculate updates at image pixels.
   * \returns A FiniteDifferenceObject pointer. */
  itkSetObjectMacro(BackwardDifferenceFunction, FiniteDifferenceFunctionType);

  using ThreadRegionType = typename DisplacementFieldFFTType::RegionType;
  using OutputImageRegionType = ThreadRegionType;
  virtual void ThreadedComputeLinearElastic(DisplacementFieldFFTPointer &, float normalizer,
                                            const ThreadRegionType & regionToProcess, int threadId);

  virtual void ThreadedComputeInverseConsistency(DisplacementFieldFFTPointer& inv0, DisplacementFieldFFTPointer& inv1,
                                                 float normalizer, const ThreadRegionType & regionToProcess, int);

  struct ThreadStruct
    {
    ICCDeformableRegistrationFilter *Filter;
    DisplacementFieldFFTPointer coeff;
    DisplacementFieldFFTPointer inverse0;
    DisplacementFieldFFTPointer inverse1;
    float normalizer_InverseConsistency;
    float normalizer_Regularization;
    };
private:
  /** Downcast the DifferenceFunction using a dynamic_cast to ensure that it is of the correct type.
    * this method will throw an exception if the function is not of the expected type. */
  ICCDeformableFunctionType *  GetForwardRegistrationFunctionType();

  const ICCDeformableFunctionType *  GetForwardRegistrationFunctionType() const;

  ICCDeformableFunctionType *  GetBackwardRegistrationFunctionType();

  const ICCDeformableFunctionType *  GetBackwardRegistrationFunctionType() const;

  static ITK_THREAD_RETURN_TYPE ComputeLinearElasticThreaderCallback(void * arg);

  static ITK_THREAD_RETURN_TYPE ComputeInverseConsistencyThreaderCallback(void * arg);

  unsigned int SplitRequestedRegion(unsigned int i, unsigned int num, OutputImageRegionType& splitRegion) override;

  /** Exp and composition type alias */
  using MultiplyByConstantType = MultiplyImageFilter<DisplacementFieldType,
                              Image<TimeStepType, DisplacementFieldType::ImageDimension> ,
                              DisplacementFieldType>;

  using FieldExponentiatorType = ExponentialDisplacementFieldImageFilter<
      DisplacementFieldType, DisplacementFieldType>;

  using VectorWarperType = WarpVectorImageFilter<
      DisplacementFieldType,
      DisplacementFieldType, DisplacementFieldType>;

  using FieldInterpolatorType = VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<
      DisplacementFieldType, double>;

  using AdderType = AddImageFilter<
      DisplacementFieldType,
      DisplacementFieldType, DisplacementFieldType>;

  using MultiplyByConstantPointer = typename MultiplyByConstantType::Pointer;
  using FieldExponentiatorPointer = typename FieldExponentiatorType::Pointer;
  using VectorWarperPointer = typename VectorWarperType::Pointer;
  using FieldInterpolatorPointer = typename FieldInterpolatorType::Pointer;
  using FieldInterpolatorOutputType = typename FieldInterpolatorType::OutputType;
  using AdderPointer = typename AdderType::Pointer;
  typename ICCDeformableFunctionType::Pointer           m_DifferenceFunction;
  typename FiniteDifferenceFunctionType::Pointer        m_BackwardDifferenceFunction;

//  MultiplyByConstantPointer m_Multiplier;
  FieldExponentiatorPointer m_Exponentiator;
  VectorWarperPointer       m_Warper;
//  AdderPointer              m_Adder;
  bool                                       m_UseFirstOrderExp;
  bool                                       m_UseConsistentIntensity;
  bool                                       m_UseConsistentLandmark;
  std::vector<DisplacementFieldPointer>      m_UpdateBuffers;
  std::vector<DisplacementFieldPointer>      m_InverseUpdateBuffers;
  std::vector<DisplacementFieldFFTPointer>   m_Coefficients;
  std::vector<AdderPointer>                  m_Adders;
  std::vector<MultiplyByConstantPointer>     m_Multipliers;
  std::vector<FFTWComplexToRealImagePointer> m_FFTc2rs;
  ComplexImagePointer                        m_sqr11, m_sqr12, m_sqr13, m_sqr22, m_sqr23, m_sqr33;
  ComplexImagePointer                        m_SmoothFilter;

  float m_Alpha, m_Gamma, m_Beta;
  float m_RegularizationWeight;
  float m_InverseWeight;
  float m_SimilarityWeight;
  float m_LandmarkWeight;
  float m_MaximumUpdateStepLength;
  float m_HarmonicPercent;
  float m_MinJac;
  float m_BackgroundFilledValue;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkICCDeformableRegistrationFilter.hxx"
#endif

#endif
