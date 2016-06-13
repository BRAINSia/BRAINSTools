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
 *  itkICCDeformableFunction.h
 *
 *  Created by Yong Qiang Zhao on 4/21/09.
 *  Copyright 2009 __UI__. All rights reserved.
 *
 */

#ifndef __itkICCDeformableFunction_h
#define __itkICCDeformableFunction_h

#include "itkPDEDeformableRegistrationFunction.h"
#include "itkVectorFFTWRealToHalfHermitianForwardFFTImageFilter.h"
#include "itkVectorFFTWHalfHermitianToRealInverseFFTImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkSpatialObject.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkLandmarkSpatialObject.h"
#include "itkPointSet.h"

namespace itk
{
template <class TFixedImage, class TMovingImage, class TDisplacementField>
class ICCDeformableFunction :
  public         PDEDeformableRegistrationFunction<TFixedImage,
                                                   TMovingImage, TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef ICCDeformableFunction Self;
  typedef PDEDeformableRegistrationFunction<
      TFixedImage, TMovingImage, TDisplacementField>    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ICCDeformableFunction,
                PDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;
  typedef typename MovingImageType::PixelType     MovingPixelType;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType    FixedImageType;
  typedef typename Superclass::FixedImagePointer FixedImagePointer;
  typedef typename FixedImageType::IndexType     IndexType;
  typedef typename FixedImageType::SizeType      SizeType;
  typedef typename FixedImageType::SpacingType   SpacingType;
  typedef typename FixedImageType::DirectionType DirectionType;
  typedef typename FixedImageType::PixelType     FixedPixelType;

  /** Deformation field type. */
  typedef typename Superclass::DisplacementFieldType DisplacementFieldType;
  typedef typename Superclass::DisplacementFieldTypePointer
    DisplacementFieldTypePointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType     TimeStepType;

  /** Interpolator type. */
  typedef double CoordRepType;
  typedef InterpolateImageFunction<
      MovingImageType, CoordRepType>                   InterpolatorMovingType;
  typedef typename InterpolatorMovingType::Pointer   InterpolatorMovingPointer;
  typedef typename InterpolatorMovingType::PointType MovingPointType;

  typedef LinearInterpolateImageFunction<
      MovingImageType, CoordRepType>                   DefaultMovingInterpolatorType;

  typedef LinearInterpolateImageFunction<
      FixedImageType, CoordRepType>                   DefaultFixedInterpolatorType;

  /** Warper type */
  typedef WarpImageFilter<
      MovingImageType,
      MovingImageType, DisplacementFieldType>           WarperMovingType;
  typedef typename WarperMovingType::Pointer WarperMovingPointer;

  typedef unsigned char                                                MaskPixelType;
  typedef Image<MaskPixelType, itkGetStaticConstMacro(ImageDimension)> MaskImageType;
  typedef typename MaskImageType::Pointer                              MaskImagePointer;

  typedef WarpImageFilter<MaskImageType,
                          MaskImageType, DisplacementFieldType>           MaskWarperType;
  typedef typename MaskWarperType::Pointer MaskWarperPointer;

  typedef std::complex<MovingPixelType> ComplexPixelType;
  typedef itk::Vector<ComplexPixelType,
                      itkGetStaticConstMacro(ImageDimension)>             DisplacementFieldFFTPixelType;
  typedef itk::Image<ComplexPixelType, itkGetStaticConstMacro(ImageDimension)>              ComplexImageType;
  typedef itk::Image<DisplacementFieldFFTPixelType, itkGetStaticConstMacro(ImageDimension)> DisplacementFieldFFTType;

  typedef typename DisplacementFieldFFTType::Pointer DisplacementFieldFFTPointer;
  typedef typename ComplexImageType::Pointer         ComplexImagePointer;

  typedef VectorFFTWRealToHalfHermitianForwardFFTImageFilter<typename TDisplacementField::PixelType,
                                                             itkGetStaticConstMacro(ImageDimension)>
    FFTWRealToComplexType;
  typedef VectorFFTWHalfHermitianToRealInverseFFTImageFilter<typename TDisplacementField::PixelType,
                                                             itkGetStaticConstMacro(ImageDimension)>
    FFTWComplexToRealType;

  typedef SubtractImageFilter<FixedImageType, FixedImageType, FixedImageType> SubtractImageType;
  typedef DerivativeImageFilter<MovingImageType, MovingImageType>             DerivativeType;
  typedef ImageRegionIterator<DisplacementFieldType>                          IterationDisplacementFieldType;
  typedef ImageRegionIterator<FixedImageType>                                 IterationImageType;
  typedef ImageRegionConstIterator<FixedImageType>                            ConstIterationImageType;
//  typedef SubtractImageFilter<DisplacementFieldType,  DisplacementFieldType, DisplacementFieldType>
// SubtractDisplacementFieldType;
//  typedef SubtractImageFilter<DisplacementFieldFFTType, DisplacementFieldFFTType, DisplacementFieldFFTType>
// SubtractDisplacementFieldFFTType;
//  typedef MultiplyByConstantImageFilter<DisplacementFieldType, float, DisplacementFieldType>
// MultiplyByConstantImageType;
//  typedef MultiplyByConstantImageFilter<DisplacementFieldFFTType, float, DisplacementFieldFFTType>
// MultiplyByConstantFFTType;
  typedef AddImageFilter<DisplacementFieldType, DisplacementFieldType, DisplacementFieldType> AddImageType;
  typedef Image<float, itkGetStaticConstMacro(ImageDimension)>                                FloatImageType;

  // typename TDisplacementField::PixelType
  /** Covariant vector type. */
  typedef CovariantVector<double, itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

  typedef SpatialObject<itkGetStaticConstMacro(ImageDimension)> MaskType;
  typedef typename MaskType::Pointer                            MaskPointer;

  typedef LandmarkSpatialObject<itkGetStaticConstMacro(ImageDimension)> LandmarkType;
  typedef typename LandmarkType::Pointer                                LandmarkPointer;

  typedef PointSet<MovingPixelType, itkGetStaticConstMacro(ImageDimension)> PointSetType;
  typedef typename PointSetType::Pointer                                    PointSetPointer;
  /** Set the moving image interpolator. */
  void SetMovingImageInterpolator( InterpolatorMovingType * ptr )
  {
    m_MovingImageInterpolator = ptr; m_MovingImageWarper->SetInterpolator( ptr );
  }

  /** Get the moving image interpolator. */
  InterpolatorMovingType * GetMovingImageInterpolator(void)
  {
    return m_MovingImageInterpolator;
  }

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void * itkNotUsed(GlobalData) ) const ITK_OVERRIDE
  {
    return m_TimeStep;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void * GetGlobalDataPointer() const ITK_OVERRIDE
  {
    GlobalDataStruct *global = new GlobalDataStruct();

    global->m_SumOfSquaredDifference  = 0.0;
    global->m_NumberOfPixelsProcessed = 0L;
    global->m_SumOfSquaredChange      = 0;
    return global;
  }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const ITK_OVERRIDE;

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration() ITK_OVERRIDE;

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */

  virtual void ComputeMetric( void *globalData);

  virtual PixelType  ComputeUpdate(const NeighborhoodType & neighborhood, void *globalData,
                                   const FloatOffsetType & offset = FloatOffsetType( 0.0) ) ITK_OVERRIDE;

  /** Get the metric value. The metric value is the mean square difference
   * in intensity between the fixed image and transforming moving image
   * computed over the the overlapping region between the two images. */
  virtual double GetMetric() const
  {
    return m_Metric;
  }

  /** Get the rms change in deformation field. */
  virtual const double & GetRMSChange() const
  {
    return m_RMSChange;
  }

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void SetIntensityDifferenceThreshold(double);

  virtual double GetIntensityDifferenceThreshold() const;

  /** Set/Get whether using consistent internsity registration**/
  virtual void SetUseConsistentIntensity(bool cl)
  {
    m_UseConsistentIntensity = cl;
  }

  virtual bool GetUseConsistentIntensity(void)
  {
    return m_UseConsistentIntensity;
  }

  virtual void SetUseConsistentLandmark(bool cl)
  {
    m_UseConsistentLandmark = cl;
  }

  virtual bool GetUseConsistentLandmark(void)
  {
    return m_UseConsistentLandmark;
  }

  /** Set/Get the maximum update step length. In Thirion this is 0.5.
   *  Setting it to 0 implies no restriction (beware of numerical
   *  instability in this case. */
  virtual void SetMaximumUpdateStepLength(double sm)
  {
    this->m_MaximumUpdateStepLength = sm;
  }

  virtual double GetMaximumUpdateStepLength() const
  {
    return this->m_MaximumUpdateStepLength;
  }

  virtual void SetBackgroundFilledValue(float bf)
  {
    this->m_BackgroundFilledValue = bf;
  }

  virtual float GetBackgroundFilledValue() const
  {
    return this->m_BackgroundFilledValue;
  }

  /** Type of available image forces */
  enum GradientType
    {
    Symmetric = 0,
    Fixed = 1,
    WarpedMoving = 2,
    MappedMoving = 3
    };

  /** Set/Get the type of used image forces */
  virtual void SetUseGradientType( GradientType gtype )
  {
    m_UseGradientType = gtype;
  }

  virtual GradientType GetUseGradientType() const
  {
    return m_UseGradientType;
  }

  DisplacementFieldTypePointer GetInverseUpdateBuffer()
  {
    return m_InverseUpdateBuffer;
  }

  DisplacementFieldTypePointer GetUpdateBuffer()
  {
    return m_UpdateBuffer;
  }

  DisplacementFieldFFTPointer GetCoefficient()
  {
    return m_Coefficient;
  }

  void SetCoefficient(DisplacementFieldFFTPointer& coeff)
  {
    m_Coefficient = coeff;
  }

  virtual void SetSimilarityWeight(double sm)
  {
    this->m_SimilarityWeight = sm;
  }

  virtual double GetSimilarityWeight() const
  {
    return this->m_SimilarityWeight;
  }

  virtual void SetLandmarkWeight(double sm)
  {
    this->m_LandmarkWeight = sm;
  }

  virtual double GetLandmarkWeight() const
  {
    return this->m_LandmarkWeight;
  }

  virtual void SetSmoothFilter(ComplexImagePointer& filter)
  {
    m_SmoothFilter = filter;
  }

  virtual void SetMovingImageMask(MaskType *mask)
  {
    m_MovingMask = mask;
  }

  virtual void SetFixedImageMask(MaskType *mask)
  {
    m_FixedMask = mask;
  }

  virtual const MaskType * GetFixedImageMask(void) const
  {
    return m_FixedMask;
  }

  virtual const MaskType * GetMovingImageMask(void) const
  {
    return m_MovingMask;
  }

  virtual MaskImageType * GetMovingWarpedMask(void)
  {
    return m_WarpedMaskImage;
  }

  virtual void SetFixedLandmark(PointSetType *lk)
  {
    m_FixedLandmark = lk;
  }

  virtual void SetMovingLandmark(PointSetType *lk)
  {
    m_MovingLandmark = lk;
  }

protected:
  ICCDeformableFunction();
  ~ICCDeformableFunction()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
    {
    double m_SumOfSquaredDifference;
    unsigned long m_NumberOfPixelsProcessed;
    double m_SumOfSquaredChange;
    };
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(ICCDeformableFunction);

  /** Cache fixed image information. */
  MovingPointType m_FixedImageOrigin;
  SpacingType     m_FixedImageSpacing;
  DirectionType   m_FixedImageDirection;

  double m_Normalizer;

  /** Function to compute derivatives of the fixed image. */
//  GradientCalculatorPointer       m_FixedImageGradientCalculator;

  /** Function to compute derivatives of the moving image (unwarped). */
//  MovingImageGradientCalculatorPointer m_MappedMovingImageGradientCalculator;

  GradientType m_UseGradientType;

  /** Function to interpolate the moving image. */
  InterpolatorMovingPointer m_MovingImageInterpolator;

  /** Filter to warp moving image for fast gradient computation. */
  WarperMovingPointer m_MovingImageWarper;
  MaskWarperPointer   m_MovingMaskImageWarper;

  /** The global timestep. */
  TimeStepType m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double m_IntensityDifferenceThreshold;

  /** Maximum update step length in pixels (default is 0.5 as in Thirion). */
  double m_MaximumUpdateStepLength;

  /** The metric value is the mean square difference in intensity between
   * the fixed image and transforming moving image computed over the
   * the overlapping region between the two images. */
  mutable double               m_Metric;
  mutable double               m_SumOfSquaredDifference;
  mutable unsigned long        m_NumberOfPixelsProcessed;
  mutable double               m_RMSChange;
  mutable double               m_SumOfSquaredChange;
  DisplacementFieldTypePointer m_InverseUpdateBuffer;
  DisplacementFieldTypePointer m_UpdateBuffer;
  DisplacementFieldFFTPointer  m_Coefficient;

  float               m_SimilarityWeight;
  float               m_LandmarkWeight;
  ComplexImagePointer m_SmoothFilter;

  /** Mutex lock to protect modification to metric. */
  mutable SimpleFastMutexLock m_MetricCalculationLock;
  MaskPointer                 m_MovingMask;
  MaskPointer                 m_FixedMask;
  PointSetPointer             m_FixedLandmark;
  PointSetPointer             m_MovingLandmark;
  MaskImagePointer            m_WarpedMaskImage;
  float                       m_FixedImageBackground;
  float                       m_MovingImageBackground;
  float                       m_BackgroundFilledValue;
  bool                        m_UseConsistentIntensity;
  bool                        m_UseConsistentLandmark;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkICCDeformableFunction.hxx"
#endif

#endif
