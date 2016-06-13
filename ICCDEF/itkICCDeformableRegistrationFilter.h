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
template <class TFixedImage, class TMovingImage, class TDisplacementField>
class ICCDeformableRegistrationFilter :
  public         PDEDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                                 TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef ICCDeformableRegistrationFilter Self;
  typedef PDEDeformableRegistrationFilter<
      TFixedImage, TMovingImage, TDisplacementField>      Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ICCDeformableRegistrationFilter,
                PDEDeformableRegistrationFilter );

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType    FixedImageType;
  typedef typename Superclass::FixedImagePointer FixedImagePointer;

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;

  /** Deformation field type. */
  typedef typename Superclass::DisplacementFieldType    DisplacementFieldType;
  typedef typename Superclass::DisplacementFieldPointer DisplacementFieldPointer;

  typedef DisplacementFieldType OutputImageType;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);
  typedef unsigned char                                                MaskPixelType;
  typedef Image<MaskPixelType, itkGetStaticConstMacro(ImageDimension)> MaskImageType;
  typedef typename MaskImageType::Pointer                              MaskImagePointer;

  typedef SpatialObject<itkGetStaticConstMacro(ImageDimension)> MaskType;
  typedef typename MaskType::Pointer                            MaskPointer;

  typedef typename DataObject::Pointer DataObjectPointer;

  typedef LandmarkSpatialObject<itkGetStaticConstMacro(ImageDimension)> LandmarkType;
  typedef typename LandmarkType::Pointer                                LandmarkPointer;

  typedef PointSet<typename MovingImageType::PixelType, itkGetStaticConstMacro(ImageDimension)> PointSetType;
  typedef typename PointSetType::Pointer                                                        PointSetPointer;

  typedef std::complex<float>                     ComplexPixelType;
  typedef Vector<ComplexPixelType, 3>             DisplacementFieldFFTPixelType;
  typedef Image<ComplexPixelType, 3>              ComplexImageType;
  typedef Image<DisplacementFieldFFTPixelType, 3> DisplacementFieldFFTType;

  typedef VectorFFTWHalfHermitianToRealInverseFFTImageFilter<typename TDisplacementField::PixelType,
                                                             3> FFTWComplexToRealImageType;
  typedef VectorFFTWRealToHalfHermitianForwardFFTImageFilter<typename TDisplacementField::PixelType,
                                                             3> FFTWRealToComplexImageType;
  typedef typename FFTWComplexToRealImageType::Pointer
    FFTWComplexToRealImagePointer;
  typedef typename FFTWRealToComplexImageType::Pointer
    FFTWRealToComplexImagePointer;

  typedef typename DisplacementFieldFFTType::Pointer DisplacementFieldFFTPointer;
  typedef typename ComplexImageType::Pointer         ComplexImagePointer;

  /** FiniteDifferenceFunction type. */
  typedef typename
    Superclass::FiniteDifferenceFunctionType               FiniteDifferenceFunctionType;

  /** Take timestep type from the FiniteDifferenceFunction. */
  typedef typename
    FiniteDifferenceFunctionType::TimeStepType             TimeStepType;

  typedef WarpImageFilter<MaskImageType,
                          MaskImageType, DisplacementFieldType>            MaskWarperType;
  typedef typename MaskWarperType::Pointer MaskWarperPointer;

  /** DemonsRegistrationFilterFunction type. */
  typedef ICCDeformableFunction<
      FixedImageType,
      MovingImageType, DisplacementFieldType>                 ICCDeformableFunctionType;
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

  virtual void CopyInputToOutput() ITK_OVERRIDE;

  virtual void PostProcessOutput() ITK_OVERRIDE { }

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
  virtual ProcessObject::DataObjectPointer MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

  virtual void SetInitialForwardDisplacementField( DisplacementFieldType * ptr )
  {
    this->SetInput( 3, ptr );
  }

  virtual void SetInitialBackwardDisplacementField( DisplacementFieldType * ptr )
  {
    this->SetInput( 4, ptr );
  }

  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  virtual void SetMovingImageMask( MaskType *mask);

  virtual const MaskType * GetMovingImageMask() const;

  virtual void SetFixedImageMask(MaskType *mask);

  virtual const MaskType * GetFixedImageMask() const;

  virtual void SetFixedLandmark(PointSetType * landmark);

  virtual void SetMovingLandmark(PointSetType * landmark);

protected:
  ICCDeformableRegistrationFilter();
  ~ICCDeformableRegistrationFilter()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  /** Initialize the state of filter and equation before each iteration. */
  virtual void InitializeIteration() ITK_OVERRIDE;

  /** This method allocates storage in m_UpdateBuffer.  It is called from
   * FiniteDifferenceFilter::GenerateData(). */
  virtual void AllocateUpdateBuffer() ITK_OVERRIDE;

  /** Apply update. */
  virtual void ApplyUpdate(const TimeStepType& dt) ITK_OVERRIDE;

  virtual TimeStepType CalculateChange() ITK_OVERRIDE;

  virtual void Initialize() ITK_OVERRIDE;

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

  typedef typename DisplacementFieldFFTType::RegionType ThreadRegionType;
  typedef ThreadRegionType                              OutputImageRegionType;
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
  ITK_DISALLOW_COPY_AND_ASSIGN(ICCDeformableRegistrationFilter);

  /** Downcast the DifferenceFunction using a dynamic_cast to ensure that it is of the correct type.
    * this method will throw an exception if the function is not of the expected type. */
  ICCDeformableFunctionType *  GetForwardRegistrationFunctionType();

  const ICCDeformableFunctionType *  GetForwardRegistrationFunctionType() const;

  ICCDeformableFunctionType *  GetBackwardRegistrationFunctionType();

  const ICCDeformableFunctionType *  GetBackwardRegistrationFunctionType() const;

  static ITK_THREAD_RETURN_TYPE ComputeLinearElasticThreaderCallback(void * arg);

  static ITK_THREAD_RETURN_TYPE ComputeInverseConsistencyThreaderCallback(void * arg);

  virtual unsigned int SplitRequestedRegion(unsigned int i, unsigned int num, OutputImageRegionType& splitRegion) ITK_OVERRIDE;

  /** Exp and composition typedefs */
  typedef MultiplyImageFilter<DisplacementFieldType,
                              Image<TimeStepType, DisplacementFieldType::ImageDimension> ,
                              DisplacementFieldType> MultiplyByConstantType;

  typedef ExponentialDisplacementFieldImageFilter<
      DisplacementFieldType, DisplacementFieldType>        FieldExponentiatorType;

  typedef WarpVectorImageFilter<
      DisplacementFieldType,
      DisplacementFieldType, DisplacementFieldType>        VectorWarperType;

  typedef VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<
      DisplacementFieldType, double>                      FieldInterpolatorType;

  typedef AddImageFilter<
      DisplacementFieldType,
      DisplacementFieldType, DisplacementFieldType>        AdderType;

  typedef typename MultiplyByConstantType::Pointer   MultiplyByConstantPointer;
  typedef typename FieldExponentiatorType::Pointer   FieldExponentiatorPointer;
  typedef typename VectorWarperType::Pointer         VectorWarperPointer;
  typedef typename FieldInterpolatorType::Pointer    FieldInterpolatorPointer;
  typedef typename FieldInterpolatorType::OutputType FieldInterpolatorOutputType;
  typedef typename AdderType::Pointer                AdderPointer;
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
