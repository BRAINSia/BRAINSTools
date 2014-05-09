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
#ifndef  __BRAINSFitHelper_h
#define  __BRAINSFitHelper_h

/**
  * \author Hans J. Johnson
  *
  * The intension of the BRIANSFitHelper is to provide a simple non-templated
  * class that can be used in other programs in a way that is very similar to
  * the command line version of the program from the SlicerExecutionModel
  * version of the BRAINSFitPrimary program.
  *
  * Almost all the command line options are available in this version, but
  * there is no need to read or write files to disk in order to use this class.
  *
  */
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>

#include "BRAINSCommonLibWin32Header.h"

// TODO:  This needs to be moved to the top, and header files moved to this
// header where needed.
#include "BRAINSFitHelperTemplate.h"
#include "BRAINSFitUtils.h"

#include "GenericTransformImage.h"
#include "ReadMask.h"
#include "BRAINSMacro.h"

#include "itkIO.h"
#include "itkFindCenterOfBrainFilter.h"
#include "itkImageRandomNonRepeatingConstIteratorWithIndex.h"

namespace itk
{
/** Method for verifying that the ordering of the transformTypes is consistent
  * with converting routines. */
extern void ValidateTransformRankOrdering(const std::vector<std::string> & transformType);
}

namespace itk
{
class BRAINSFitHelper : public Object
{
public:
  /** Standard class typedefs. */
  typedef BRAINSFitHelper          Self;
  typedef ProcessObject            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef float                         PixelType;
  typedef itk::Image<PixelType, 3>      FixedImageType;
  typedef FixedImageType::ConstPointer  FixedImageConstPointer;
  typedef FixedImageType::Pointer       FixedImagePointer;

  typedef itk::Image<PixelType, 3>      MovingImageType;
  typedef MovingImageType::ConstPointer MovingImageConstPointer;
  typedef MovingImageType::Pointer      MovingImagePointer;

  /** Constants for the image dimensions */
  itkStaticConstMacro(FixedImageDimension, unsigned int, FixedImageType::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int, MovingImageType::ImageDimension);

  typedef itk::CompositeTransform<double, MovingImageDimension>       CompositeTransformType;

  typedef SpatialObject<itkGetStaticConstMacro(FixedImageDimension)>  FixedBinaryVolumeType;
  typedef SpatialObject<itkGetStaticConstMacro(MovingImageDimension)> MovingBinaryVolumeType;
  typedef FixedBinaryVolumeType::Pointer                              FixedBinaryVolumePointer;
  typedef MovingBinaryVolumeType::Pointer                             MovingBinaryVolumePointer;

  typedef itk::BRAINSFitHelperTemplate<FixedImageType, MovingImageType> HelperType;

  typedef HelperType::MetricType                             GenericMetricType;

  typedef itk::AffineTransform<double, 3>                                                         AffineTransformType;
  typedef itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, AffineTransformType>  AffineRegistrationType;
  typedef AffineRegistrationType::MetricSamplingStrategyType                                      SamplingStrategyType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BRAINSFitHelper, ProcessObject);

  /** Set/Get the Fixed image. */
  itkSetObjectMacro(FixedVolume, FixedImageType);
  itkGetConstObjectMacro(FixedVolume, FixedImageType);

  /** Set/Get the Moving image. */
  itkSetObjectMacro(MovingVolume, MovingImageType)
  itkGetConstObjectMacro(MovingVolume, MovingImageType);

  /** The preprocessedMoving volume SHOULD NOT BE SET, you can get it out of the
    *  algorithm.*/
  itkGetConstObjectMacro(PreprocessedMovingVolume, MovingImageType);

  itkSetObjectMacro(FixedBinaryVolume, FixedBinaryVolumeType);
  itkGetModifiableObjectMacro(FixedBinaryVolume, FixedBinaryVolumeType);
  itkSetObjectMacro(MovingBinaryVolume, MovingBinaryVolumeType);
  itkGetModifiableObjectMacro(MovingBinaryVolume, MovingBinaryVolumeType);

  itkSetMacro(OutputFixedVolumeROI,  std::string);
  itkGetConstMacro(OutputFixedVolumeROI,  std::string);
  itkSetMacro(OutputMovingVolumeROI, std::string);
  itkGetConstMacro(OutputMovingVolumeROI, std::string);

  // TODO:  This should be converted to use the
  //       interpolation mechanisms from GenericTransform
  typedef enum
    {
    LINEAR_INTERP = 0,
    WINDOWSINC_INTERP = 1
    } InterpolationType;

  itkSetMacro(NumberOfSamples,                   unsigned int);
  itkGetConstMacro(NumberOfSamples,              unsigned int);
  itkSetMacro(SamplingPercentage,                double);
  itkGetConstMacro(SamplingPercentage,           double);
  itkSetMacro(NumberOfHistogramBins,             unsigned int);
  itkGetConstMacro(NumberOfHistogramBins,        unsigned int);
  itkSetMacro(NumberOfMatchPoints,               unsigned int);
  itkGetConstMacro(NumberOfMatchPoints,          unsigned int);
  VECTORitkSetMacro(NumberOfIterations,   std::vector<int> /**/);
  VECTORitkSetMacro(MinimumStepLength,    std::vector<double> ); // It will not be used by ITKv4 registration, but it is kept for backward compatibility.
  itkSetMacro(MaximumStepLength,             double);
  itkGetConstMacro(MaximumStepLength,             double);
  itkSetMacro(RelaxationFactor,              double);
  itkGetConstMacro(RelaxationFactor,              double);
  itkSetMacro(TranslationScale,              double);
  itkGetConstMacro(TranslationScale,              double);
  itkSetMacro(ReproportionScale,             double);
  itkGetConstMacro(ReproportionScale,             double);
  itkSetMacro(SkewScale,                     double);
  itkGetConstMacro(SkewScale,                     double);
  itkSetMacro(UseCachingOfBSplineWeightsMode, std::string);
  itkGetConstMacro(UseCachingOfBSplineWeightsMode, std::string);
  itkSetMacro(CostFunctionConvergenceFactor, double);
  itkGetConstMacro(CostFunctionConvergenceFactor, double);
  itkSetMacro(ProjectedGradientTolerance,    double);
  itkGetConstMacro(ProjectedGradientTolerance,    double);
  itkSetMacro(MaxBSplineDisplacement,        double);
  itkGetConstMacro(MaxBSplineDisplacement,        double);
  itkSetMacro(BackgroundFillValue,           double);
  itkGetConstMacro(BackgroundFillValue,           double);
  VECTORitkSetMacro(TransformType, std::vector<std::string> );
  itkSetMacro(InitializeTransformMode, std::string);
  itkGetConstMacro(InitializeTransformMode, std::string);
  itkSetMacro(MaskInferiorCutOffFromCenter, double);
  itkGetConstMacro(MaskInferiorCutOffFromCenter, double);
  itkSetMacro(CurrentGenericTransform,  CompositeTransformType::Pointer);
  itkGetConstMacro(CurrentGenericTransform,  CompositeTransformType::Pointer);
  VECTORitkSetMacro(SplineGridSize, std::vector<int>       );

  itkGetConstMacro(ActualNumberOfIterations,      unsigned int);
  itkGetConstMacro(PermittedNumberOfIterations,   unsigned int);

  itkGetConstMacro(FinalMetricValue,         double);
  /** Set/Get the Debugging level for filter verboseness */
  itkSetMacro(DebugLevel, unsigned int);
  itkGetConstMacro(DebugLevel, unsigned int);
  itkSetMacro(DisplayDeformedImage, bool);
  itkGetConstMacro(DisplayDeformedImage, bool);
  itkSetMacro(PromptUserAfterDisplay, bool);
  itkGetConstMacro(PromptUserAfterDisplay, bool);
  itkSetMacro(ObserveIterations,        bool);
  itkGetConstMacro(ObserveIterations,        bool);
  itkSetMacro(UseROIBSpline, bool);
  itkGetConstMacro(UseROIBSpline, bool);
  itkSetMacro(MetricSeed, int);
  itkGetConstMacro(MetricSeed, int);

  void SetSamplingStrategy(std::string strategy)
  {
  if( strategy == "Random" )
    {
    m_SamplingStrategy = AffineRegistrationType::RANDOM;
    }
  else if( strategy == "Regular" )
    {
    m_SamplingStrategy = AffineRegistrationType::REGULAR;
    }
  else if( (strategy == "None") || (strategy == "") )
    {
    m_SamplingStrategy = AffineRegistrationType::NONE;
    }
  else
    {
    std::cout << "ERROR: samplingStrategy is incorrectly specified" << std::endl;
    exit( -1 );
    }
  }

  /** Method to set the Permission to vary by level  */
  void SetPermitParameterVariation(std::vector<int> perms)
  {
    m_PermitParameterVariation.resize( perms.size() );
    for( unsigned int i = 0; i < perms.size(); ++i )
      {
      m_PermitParameterVariation[i] = perms[i];
      }
  }

  itkSetMacro(HistogramMatch, bool);
  itkGetConstMacro(HistogramMatch, bool);

  itkSetMacro(RemoveIntensityOutliers, float);
  itkGetConstMacro(RemoveIntensityOutliers, float);

  itkSetMacro(CostMetric, std::string);
  itkGetConstMacro(CostMetric, std::string);

  itkSetMacro(ForceMINumberOfThreads, int);
  itkGetConstMacro(ForceMINumberOfThreads, int);

  itkSetMacro(NormalizeInputImages, bool);

  /** Method that initiates the registration. */
  void Update(void);

  void PrintCommandLine(const bool dumpTempVolumes, const std::string & suffix) const;

protected:
  BRAINSFitHelper();
  virtual ~BRAINSFitHelper()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
    * the registration. */
  void  GenerateData();

private:

  BRAINSFitHelper(const Self &); // purposely not implemented
  void operator=(const Self &);  // purposely not implemented

  template <class TLocalCostMetric>
  void SetupRegistration(GenericMetricType *localCostMetric);

  template <class TLocalCostMetric>
  void RunRegistration();

  FixedImagePointer  m_FixedVolume;
  MovingImagePointer m_MovingVolume;
  MovingImagePointer m_PreprocessedMovingVolume;

  FixedBinaryVolumePointer  m_FixedBinaryVolume;
  MovingBinaryVolumePointer m_MovingBinaryVolume;
  std::string               m_OutputFixedVolumeROI;
  std::string               m_OutputMovingVolumeROI;
  std::vector<int>          m_PermitParameterVariation;

  unsigned int m_NumberOfSamples;
  double       m_SamplingPercentage;
  unsigned int m_NumberOfHistogramBins;
  bool         m_HistogramMatch;
  float        m_RemoveIntensityOutliers;
  unsigned int m_NumberOfMatchPoints;
  // TODO:  Would be better to have unsigned int
  std::vector<int>         m_NumberOfIterations;
  double                   m_MaximumStepLength;
  std::vector<double>      m_MinimumStepLength;
  double                   m_RelaxationFactor;
  double                   m_TranslationScale;
  double                   m_ReproportionScale;
  double                   m_SkewScale;
  std::string              m_UseCachingOfBSplineWeightsMode;
  double                   m_BackgroundFillValue;
  std::vector<std::string> m_TransformType;
  std::string              m_InitializeTransformMode;
  double                   m_MaskInferiorCutOffFromCenter;
  std::vector<int>         m_SplineGridSize;
  double                   m_CostFunctionConvergenceFactor;
  double                   m_ProjectedGradientTolerance;
  double                   m_MaxBSplineDisplacement;
  unsigned int             m_ActualNumberOfIterations;
  unsigned int             m_PermittedNumberOfIterations;
  unsigned int                               m_DebugLevel;
  CompositeTransformType::Pointer            m_CurrentGenericTransform;
  bool                                       m_DisplayDeformedImage;
  bool                                       m_PromptUserAfterDisplay;
  double                                     m_FinalMetricValue;
  bool                                       m_ObserveIterations;
  std::string                                m_CostMetric;
  bool                                       m_UseROIBSpline;
  int                                        m_MetricSeed;
  itk::Object::Pointer                       m_Helper;
  SamplingStrategyType                       m_SamplingStrategy;
  bool                                       m_NormalizeInputImages;
  // DEBUG OPTION:
  int m_ForceMINumberOfThreads;
};  // end BRAINSFitHelper class

template <class TLocalCostMetric>
void
BRAINSFitHelper::SetupRegistration(GenericMetricType *localCostMetric)
{
  typedef typename TLocalCostMetric::FixedSampledPointSetType                    MetricSamplePointSetType;

  localCostMetric->SetVirtualDomainFromImage(this->m_FixedVolume);

  localCostMetric->SetFixedImage(this->m_FixedVolume);
  localCostMetric->SetMovingImage(this->m_PreprocessedMovingVolume);

  // The preference is just set SamplingPercentage and pass that to the registrationMethodv4 directly.
  // However, the NumberOfSamples option is kept for backward compatibility,
  // and if it is set, it overwrites the SamplingPercentage option.
  if(this->m_NumberOfSamples == 0 || this->m_SamplingPercentage == 0)
    {
    this->m_SamplingPercentage = 1;
    }

  if(this->m_NumberOfSamples > 0)
    {
    const unsigned long numberOfAllSamples = this->m_FixedVolume->GetBufferedRegion().GetNumberOfPixels();
    this->m_SamplingPercentage = static_cast<double>( this->m_NumberOfSamples )/numberOfAllSamples;
    if( this->m_SamplingStrategy == AffineRegistrationType::NONE )
      {
      this->m_SamplingStrategy = AffineRegistrationType::REGULAR;
      }
    }

  if( this->m_MovingBinaryVolume.IsNotNull() )
    {
    localCostMetric->SetMovingImageMask(this->m_MovingBinaryVolume);
    }
  if( this->m_FixedBinaryVolume.IsNotNull() )
    {
    if( this->m_SamplingPercentage != 1 )
      {
      // In this case the registration framework does not do sampling inside the mask area.
      // It picks samples from the whole image, and those samples that are inside the mask are selected by metric.
      // However, we want to pick all of our intented samples from the mask area, so we do sampling here.

      // First overwrite the sampling strategy to be none
      this->m_SamplingStrategy = AffineRegistrationType::NONE;

      // then pick the samples inside the mask and pass them to metric
      typename MetricSamplePointSetType::Pointer samplePointSet = MetricSamplePointSetType::New();
      samplePointSet->Initialize();

      typedef typename MetricSamplePointSetType::PointType SamplePointType;

      itk::ImageRandomNonRepeatingConstIteratorWithIndex<FixedImageType> NRit( this->m_FixedVolume,
                                                                              this->m_FixedVolume->GetBufferedRegion() );
      const unsigned long numberOfAllSamples = this->m_FixedVolume->GetBufferedRegion().GetNumberOfPixels();
      NRit.SetNumberOfSamples( numberOfAllSamples );
      NRit.GoToBegin();

      const unsigned long sampleCount = static_cast<unsigned long>(vcl_ceil( numberOfAllSamples * this->m_SamplingPercentage ) );
      unsigned long index = 0;
      size_t currentCount = 0;

      while( ( !NRit.IsAtEnd() ) && ( currentCount < sampleCount ) )
        {
        SamplePointType testPoint;
        this->m_FixedVolume->TransformIndexToPhysicalPoint(NRit.GetIndex(), testPoint);
        if( this->m_FixedBinaryVolume->IsInside(testPoint) )
          {
          samplePointSet->SetPoint( index, testPoint );
          ++index;
          ++currentCount;
          }
        ++NRit;
        }

      if( samplePointSet.IsNull() )
        {
        itkGenericExceptionMacro("samplePointSet is empty.");
        }

      localCostMetric->SetUseFixedSampledPointSet( true );
      localCostMetric->SetFixedSampledPointSet( samplePointSet );
      }
    else
      {
      // The samplingPercentage will be passed to registration filter as it is (i.e. 1),
      // No sampling happens, but metric just accepts the points that are inside the mask area.
      localCostMetric->SetFixedImageMask(this->m_FixedBinaryVolume);
      }
    }

  typename HelperType::Pointer
  myHelper = BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::New();
  myHelper->SetTransformType(this->m_TransformType);
  myHelper->SetFixedVolume(this->m_FixedVolume);
  myHelper->SetMovingVolume(this->m_PreprocessedMovingVolume);
  myHelper->SetHistogramMatch(this->m_HistogramMatch);
  myHelper->SetRemoveIntensityOutliers(this->m_RemoveIntensityOutliers);
  myHelper->SetNumberOfMatchPoints(this->m_NumberOfMatchPoints);
  myHelper->SetFixedBinaryVolume(this->m_FixedBinaryVolume);
  myHelper->SetMovingBinaryVolume(this->m_MovingBinaryVolume);
  myHelper->SetOutputFixedVolumeROI(this->m_OutputFixedVolumeROI);
  myHelper->SetOutputMovingVolumeROI(this->m_OutputMovingVolumeROI);
  myHelper->SetPermitParameterVariation(this->m_PermitParameterVariation);
  myHelper->SetSamplingPercentage(this->m_SamplingPercentage);
  myHelper->SetNumberOfHistogramBins(this->m_NumberOfHistogramBins);
  myHelper->SetNumberOfIterations(this->m_NumberOfIterations);
  myHelper->SetMaximumStepLength(this->m_MaximumStepLength);
  myHelper->SetRelaxationFactor(this->m_RelaxationFactor);
  myHelper->SetTranslationScale(this->m_TranslationScale);
  myHelper->SetReproportionScale(this->m_ReproportionScale);
  myHelper->SetSkewScale(this->m_SkewScale);
  myHelper->SetBackgroundFillValue(this->m_BackgroundFillValue);
  myHelper->SetInitializeTransformMode(this->m_InitializeTransformMode);
  myHelper->SetMaskInferiorCutOffFromCenter(this->m_MaskInferiorCutOffFromCenter);
  myHelper->SetCurrentGenericTransform(this->m_CurrentGenericTransform);
  myHelper->SetSplineGridSize(this->m_SplineGridSize);
  myHelper->SetCostFunctionConvergenceFactor(this->m_CostFunctionConvergenceFactor);
  myHelper->SetProjectedGradientTolerance(this->m_ProjectedGradientTolerance);
  myHelper->SetMaxBSplineDisplacement(this->m_MaxBSplineDisplacement);
  myHelper->SetDisplayDeformedImage(this->m_DisplayDeformedImage);
  myHelper->SetPromptUserAfterDisplay(this->m_PromptUserAfterDisplay);
  myHelper->SetDebugLevel(this->m_DebugLevel);
  myHelper->SetCostMetricObject(localCostMetric);
  myHelper->SetForceMINumberOfThreads(this->m_ForceMINumberOfThreads);
  myHelper->SetUseROIBSpline(this->m_UseROIBSpline);
  myHelper->SetSamplingStrategy(this->m_SamplingStrategy);
  if( this->m_DebugLevel > 7 )
    {
    this->PrintCommandLine(true, "BF");
    }
  this->m_Helper = static_cast<itk::Object *>(myHelper.GetPointer() );
}

template <class TLocalCostMetric>
void
BRAINSFitHelper::RunRegistration()
{
  const typename HelperType::Pointer myHelper( dynamic_cast<HelperType *>(this->m_Helper.GetPointer() ) );
  if( myHelper.IsNull() )
    {
    std::cout << "ERROR:  Invalid BRAINSFitHelper conversion" << __FILE__ << " " << __LINE__ << std::endl;
    }
  myHelper->Update();
  this->m_CurrentGenericTransform = myHelper->GetCurrentGenericTransform();
  this->m_ActualNumberOfIterations = myHelper->GetActualNumberOfIterations();
  this->m_PermittedNumberOfIterations = myHelper->GetPermittedNumberOfIterations();
  //this->m_GenericTransformList.resize(myHelper->GetGenericTransformListPtr()->size() );
  //std::copy(myHelper->GetGenericTransformListPtr()->begin(),
  //          myHelper->GetGenericTransformListPtr()->end(), this->m_GenericTransformList.begin() );

/*
  // Find the final metric value based on the used metric type and returned output transform
  // typedef typename HelperType::MetricType    GenericMetricType;
  const typename GenericMetricType::Pointer returnMetric = myHelper->GetModifiableCostMetricObject();
  typename TLocalCostMetric::Pointer finalMetric = dynamic_cast<TLocalCostMetric *>(returnMetric.GetPointer() );
  if( finalMetric.IsNull() )
    {
    std::cout << "ERROR:  Invalid CostMetric conversion" << __FILE__ << " " << __LINE__ << std::endl;
    }
  finalMetric->SetTransform(this->m_CurrentGenericTransform);
  finalMetric->Initialize();
  this->m_FinalMetricValue = finalMetric->GetValue();
*/
  this->m_FinalMetricValue = -123.456789;
}
}   // end namespace itk

#endif  // __BRAINSFITHELPER__
