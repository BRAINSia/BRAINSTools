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
#ifndef __BRAINSFitHelper_h
#define __BRAINSFitHelper_h

/**
 * \author Hans J. Johnson
 *
 * The intention of the BRIANSFitHelper is to provide a simple non-templated
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

// INFO:  This needs to be moved to the top, and header files moved to this
// header where needed.
#include "BRAINSFitHelperTemplate.h"
#include "BRAINSFitUtils.h"

#include "GenericTransformImage.h"
#include "ReadMask.h"
#include "BRAINSMacro.h"

#include "itkIO.h"
#include "itkFindCenterOfBrainFilter.h"
#include "itkImageRandomNonRepeatingConstIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

namespace itk
{
class BRAINSFitHelper : public Object
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(BRAINSFitHelper);

  /** Standard class type alias. */
  using Self = BRAINSFitHelper;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using RealType = double;
  using PixelType = float;
  using FixedImageType = itk::Image<PixelType, 3>;
  using FixedImageConstPointer = FixedImageType::ConstPointer;
  using FixedImagePointer = FixedImageType::Pointer;

  using MovingImageType = itk::Image<PixelType, 3>;
  using MovingImageConstPointer = MovingImageType::ConstPointer;
  using MovingImagePointer = MovingImageType::Pointer;

  /** Constants for the image dimensions */
  static constexpr unsigned int FixedImageDimension = FixedImageType::ImageDimension;
  static constexpr unsigned int MovingImageDimension = MovingImageType::ImageDimension;

  using CompositeTransformType = itk::CompositeTransform<RealType, MovingImageDimension>;

  using FixedBinaryVolumeType = SpatialObject<Self::FixedImageDimension>;
  using MovingBinaryVolumeType = SpatialObject<Self::MovingImageDimension>;
  using FixedBinaryVolumePointer = FixedBinaryVolumeType::Pointer;
  using MovingBinaryVolumePointer = MovingBinaryVolumeType::Pointer;

  using HelperType = itk::BRAINSFitHelperTemplate<FixedImageType, MovingImageType>;

  using GenericMetricType = HelperType::MetricType;
  using MultiMetricType = HelperType::MultiMetricType;

  using AffineRegistrationType = itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType>;
  using AffineTransformType = itk::AffineTransform<RealType, 3>;
  using SamplingStrategyType = AffineRegistrationType::MetricSamplingStrategyEnum;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(BRAINSFitHelper);

  /** Set/Get the Fixed image. */
  itkSetObjectMacro(FixedVolume, FixedImageType);
  itkGetConstObjectMacro(FixedVolume, FixedImageType);

  /** Set/Get the Fixed image. */
  itkSetObjectMacro(FixedVolume2, FixedImageType);
  itkGetConstObjectMacro(FixedVolume2, FixedImageType);

  /** Set/Get the Moving image. */
  itkSetObjectMacro(MovingVolume, MovingImageType);
  itkGetConstObjectMacro(MovingVolume, MovingImageType);

  /** Set/Get the Moving image. */
  itkSetObjectMacro(MovingVolume2, MovingImageType);
  itkGetConstObjectMacro(MovingVolume2, MovingImageType);

  /** The preprocessedMoving volume SHOULD NOT BE SET, you can get it out of the
   *  algorithm.*/
  itkGetConstObjectMacro(PreprocessedMovingVolume, MovingImageType);

  /** The preprocessedMoving2 volume SHOULD NOT BE SET, you can get it out of the
   *  algorithm.*/
  itkGetConstObjectMacro(PreprocessedMovingVolume2, MovingImageType);

  itkSetObjectMacro(FixedBinaryVolume, FixedBinaryVolumeType);
  itkGetModifiableObjectMacro(FixedBinaryVolume, FixedBinaryVolumeType);
  itkSetObjectMacro(MovingBinaryVolume, MovingBinaryVolumeType);
  itkGetModifiableObjectMacro(MovingBinaryVolume, MovingBinaryVolumeType);

  itkSetMacro(OutputFixedVolumeROI, std::string);
  itkGetConstMacro(OutputFixedVolumeROI, std::string);
  itkSetMacro(OutputMovingVolumeROI, std::string);
  itkGetConstMacro(OutputMovingVolumeROI, std::string);

  // INFO:  This should be converted to use the
  //       interpolation mechanisms from GenericTransform
  using InterpolationType = enum

  {

    LINEAR_INTERP = 0,

    WINDOWSINC_INTERP = 1

  };

  itkSetMacro(SamplingPercentage, RealType);
  itkGetConstMacro(SamplingPercentage, RealType);
  itkSetMacro(NumberOfHistogramBins, unsigned int);
  itkGetConstMacro(NumberOfHistogramBins, unsigned int);
  itkSetMacro(NumberOfMatchPoints, unsigned int);
  itkGetConstMacro(NumberOfMatchPoints, unsigned int);
  VECTORitkSetMacro(NumberOfIterations, std::vector<int> /**/);
  VECTORitkSetMacro(MinimumStepLength, std::vector<double>);
  itkSetMacro(MaximumStepLength, double);
  itkGetConstMacro(MaximumStepLength, double);
  itkSetMacro(RelaxationFactor, double);
  itkGetConstMacro(RelaxationFactor, double);
  itkSetMacro(TranslationScale, double);
  itkGetConstMacro(TranslationScale, double);
  itkSetMacro(ReproportionScale, double);
  itkGetConstMacro(ReproportionScale, double);
  itkSetMacro(SkewScale, double);
  itkGetConstMacro(SkewScale, double);
  itkSetMacro(CostFunctionConvergenceFactor, double);
  itkGetConstMacro(CostFunctionConvergenceFactor, double);
  itkSetMacro(ProjectedGradientTolerance, double);
  itkGetConstMacro(ProjectedGradientTolerance, double);
  itkSetMacro(MaxBSplineDisplacement, double);
  itkGetConstMacro(MaxBSplineDisplacement, double);
  itkSetMacro(BackgroundFillValue, double);
  itkGetConstMacro(BackgroundFillValue, double);
  VECTORitkSetMacro(TransformType, std::vector<std::string>);
  itkSetMacro(InitializeTransformMode, std::string);
  itkGetConstMacro(InitializeTransformMode, std::string);
  itkSetMacro(MaskInferiorCutOffFromCenter, double);
  itkGetConstMacro(MaskInferiorCutOffFromCenter, double);
  itkSetMacro(MaximumNumberOfEvaluations, int);
  itkGetConstMacro(MaximumNumberOfEvaluations, int);
  itkSetMacro(MaximumNumberOfCorrections, int);
  itkGetConstMacro(MaximumNumberOfCorrections, int);
  itkSetMacro(CurrentGenericTransform, CompositeTransformType::Pointer);
  itkGetConstMacro(CurrentGenericTransform, CompositeTransformType::Pointer);
  itkSetMacro(RestoreState, CompositeTransformType::Pointer);
  itkGetConstMacro(RestoreState, CompositeTransformType::Pointer);
  VECTORitkSetMacro(SplineGridSize, std::vector<int>);

  itkGetConstMacro(ActualNumberOfIterations, unsigned int);
  itkGetConstMacro(PermittedNumberOfIterations, unsigned int);

  itkGetConstMacro(FinalMetricValue, double);
  /** Set/Get the Debugging level for filter verboseness */
  itkSetMacro(DebugLevel, unsigned int);
  itkGetConstMacro(DebugLevel, unsigned int);
  itkSetMacro(DisplayDeformedImage, bool);
  itkGetConstMacro(DisplayDeformedImage, bool);
  itkSetMacro(PromptUserAfterDisplay, bool);
  itkGetConstMacro(PromptUserAfterDisplay, bool);
  itkSetMacro(ObserveIterations, bool);
  itkGetConstMacro(ObserveIterations, bool);
  itkSetMacro(UseROIBSpline, bool);
  itkGetConstMacro(UseROIBSpline, bool);
  itkSetMacro(WriteOutputTransformInFloat, bool);
  itkGetConstMacro(WriteOutputTransformInFloat, bool);

  void
  SetSamplingStrategy(const std::string & strategy)
  {
    if (strategy == "Random")
    {
      m_SamplingStrategy = AffineRegistrationType::MetricSamplingStrategyEnum::RANDOM;
    }
    else
    {
      std::cout << "ERROR: samplingStrategy is incorrectly specified, only Random supported" << std::endl;
      exit(-1);
    }
  }

  itkSetMacro(HistogramMatch, bool);
  itkGetConstMacro(HistogramMatch, bool);

  itkSetMacro(RemoveIntensityOutliers, float);
  itkGetConstMacro(RemoveIntensityOutliers, float);

  itkSetMacro(CostMetricName, std::string);
  itkGetConstMacro(CostMetricName, std::string);

  itkSetMacro(SaveState, std::string);
  itkGetConstMacro(SaveState, std::string);

  itkSetMacro(NormalizeInputImages, bool);
  itkSetMacro(InitializeRegistrationByCurrentGenericTransform, bool);

  /** BRAINSFit uses ANTs to run only 3 levels of SyN registration unless SyNFull flag is false!
   * In that case only a single level SyN registration is run.
   */
  itkSetMacro(SyNFull, bool);
  itkGetConstMacro(SyNFull, bool);

  /** Method that initiates the registration. */
  void
  Update();

  void
  PrintCommandLine(bool dumpTempVolumes, const std::string & suffix) const;

protected:
  BRAINSFitHelper();
  ~BRAINSFitHelper() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  void
  GenerateData();

private:
  template <typename TLocalCostMetric>
  void
  SetupRegistration(GenericMetricType * costMetric);

  template <typename TLocalCostMetric>
  void
  RunRegistration();

  FixedImagePointer  m_FixedVolume;
  FixedImagePointer  m_FixedVolume2; // For multi-modal SyN
  MovingImagePointer m_MovingVolume;
  MovingImagePointer m_MovingVolume2; // For multi-modal SyN
  MovingImagePointer m_PreprocessedMovingVolume;
  MovingImagePointer m_PreprocessedMovingVolume2; // For multi-modal SyN

  FixedBinaryVolumePointer  m_FixedBinaryVolume;
  FixedBinaryVolumePointer  m_FixedBinaryVolume2; // For multi-modal SyN
  MovingBinaryVolumePointer m_MovingBinaryVolume;
  MovingBinaryVolumePointer m_MovingBinaryVolume2; // For multi-modal SyN
  std::string               m_OutputFixedVolumeROI;
  std::string               m_OutputMovingVolumeROI;

  RealType     m_SamplingPercentage{ 1.0 };
  unsigned int m_NumberOfHistogramBins{ 50 };
  bool         m_HistogramMatch{ false };
  float        m_RemoveIntensityOutliers{ 0.00 };
  unsigned int m_NumberOfMatchPoints{ 10 };
  // INFO:  Would be better to have unsigned int
  std::vector<int>                m_NumberOfIterations;
  double                          m_MaximumStepLength{ 0.2 };
  std::vector<double>             m_MinimumStepLength;
  double                          m_RelaxationFactor{ 0.5 };
  double                          m_TranslationScale{ 1000.0 };
  double                          m_ReproportionScale{ 1.0 };
  double                          m_SkewScale{ 1.0 };
  double                          m_BackgroundFillValue{ 0.0 };
  std::vector<std::string>        m_TransformType;
  std::string                     m_InitializeTransformMode;
  double                          m_MaskInferiorCutOffFromCenter{ 1000 };
  std::vector<int>                m_SplineGridSize;
  double                          m_CostFunctionConvergenceFactor{ 1e+9 };
  double                          m_ProjectedGradientTolerance{ 1e-5 };
  double                          m_MaxBSplineDisplacement{ 0.0 };
  unsigned int                    m_ActualNumberOfIterations{ 0 };
  unsigned int                    m_PermittedNumberOfIterations{ 0 };
  unsigned int                    m_DebugLevel{ 0 };
  CompositeTransformType::Pointer m_CurrentGenericTransform;
  CompositeTransformType::Pointer m_RestoreState;
  bool                            m_DisplayDeformedImage{ false };
  bool                            m_PromptUserAfterDisplay{ false };
  double                          m_FinalMetricValue{ 0.0 };
  bool                            m_ObserveIterations{ true };
  std::string                     m_CostMetricName;
  std::string                     m_SaveState;
  bool                            m_UseROIBSpline{ false };
  itk::Object::Pointer            m_Helper;
  SamplingStrategyType            m_SamplingStrategy;
  bool                            m_NormalizeInputImages{ false };
  bool                            m_InitializeRegistrationByCurrentGenericTransform{ true };
  int                             m_MaximumNumberOfEvaluations{ 900 };
  int                             m_MaximumNumberOfCorrections{ 12 };
  bool                            m_SyNFull{ true };
  bool                            m_WriteOutputTransformInFloat{ false };
}; // end BRAINSFitHelper class

template <typename TLocalCostMetric>
void
BRAINSFitHelper::SetupRegistration(GenericMetricType * costMetric)
{
  using MetricSamplePointSetType = typename TLocalCostMetric::FixedSampledPointSetType;

  typename TLocalCostMetric::Pointer localCostMetric = dynamic_cast<TLocalCostMetric *>(costMetric);
  if (localCostMetric.IsNull())
  {
    itkGenericExceptionMacro("ERROR in metric type conversion");
  }

  localCostMetric->SetVirtualDomainFromImage(this->m_FixedVolume);

  localCostMetric->SetFixedImage(this->m_FixedVolume);
  localCostMetric->SetMovingImage(this->m_PreprocessedMovingVolume);

  if (this->m_SamplingPercentage <= 0.0)
  {
    itkGenericExceptionMacro("ERROR: SamplingPercentage can not be less than or equal to zero");
  }

  if (this->m_MovingBinaryVolume.IsNotNull())
  {
    localCostMetric->SetMovingImageMask(this->m_MovingBinaryVolume);
  }
  if (this->m_FixedBinaryVolume.IsNotNull())
  {
    // In this case the registration framework does not do sampling inside the mask area.
    // It picks samples from the whole image, and those samples that are inside the mask are selected by metric.
    // However, we want to pick all of our intended samples from the mask area, so we do sampling here.

    // First overwrite the sampling strategy to be none
    this->m_SamplingStrategy = AffineRegistrationType::MetricSamplingStrategyEnum::NONE;

    // then pick the samples inside the mask and pass them to metric
    typename MetricSamplePointSetType::Pointer samplePointSet = MetricSamplePointSetType::New();
    samplePointSet->Initialize();

    using SamplePointType = typename MetricSamplePointSetType::PointType;
    const unsigned long numberOfAllSamples = this->m_FixedVolume->GetBufferedRegion().GetNumberOfPixels();

    const auto sampleCount = static_cast<unsigned long>(std::ceil(numberOfAllSamples * this->m_SamplingPercentage));

    using RandomizerType = Statistics::MersenneTwisterRandomVariateGenerator;
    typename RandomizerType::Pointer randomizer = RandomizerType::New();
    randomizer->SetSeed(1234);

    itk::ImageRandomNonRepeatingConstIteratorWithIndex<FixedImageType> NRit(this->m_FixedVolume,
                                                                            this->m_FixedVolume->GetBufferedRegion());
    NRit.ReinitializeSeed(121212);

    const typename FixedImageType::SpacingType oneThirdVirtualSpacing = this->m_FixedVolume->GetSpacing() / 3.0;
    NRit.SetNumberOfSamples(numberOfAllSamples); // Take random samples from entire image.
    NRit.GoToBegin();
    unsigned long samplesInsideMask = 0;
    while (!NRit.IsAtEnd() && (samplesInsideMask < sampleCount))
    {
      SamplePointType testPoint;
      this->m_FixedVolume->TransformIndexToPhysicalPoint(NRit.GetIndex(), testPoint);
      if (this->m_FixedBinaryVolume->IsInsideInWorldSpace(testPoint))
      {
        // randomly perturb the point within a voxel (approximately)
        for (unsigned int d = 0; d < FixedImageDimension; d++)
        {
          testPoint[d] += randomizer->GetNormalVariate() * oneThirdVirtualSpacing[d];
        }
        samplePointSet->SetPoint(samplesInsideMask, testPoint);
        ++samplesInsideMask;
      }
      ++NRit;
    }

    if (samplePointSet.IsNull())
    {
      itkGenericExceptionMacro("samplePointSet is empty.");
    }

    localCostMetric->SetUseSampledPointSet(true);
    localCostMetric->SetFixedSampledPointSet(samplePointSet);
  }

  unsigned int numberOfinputImageSets = 1;
  if (m_FixedVolume2.IsNotNull() && m_MovingVolume2.IsNotNull())
  {
    numberOfinputImageSets = 2;
    std::cout << "Multi-modal registration is run! Number of image pairs: " << numberOfinputImageSets << std::endl;
    std::cout << "In BRAINSFit the same metric is used for both modalities: " << this->m_CostMetricName << std::endl;
  }
  typename MultiMetricType::Pointer multiMetric = MultiMetricType::New();
  for (unsigned int i = 0; i < numberOfinputImageSets; i++)
  {
    multiMetric->AddMetric(localCostMetric); // In the case of multi-modality,
                                             // The same metric is used for both modalities.
  }

  typename HelperType::Pointer myHelper = BRAINSFitHelperTemplate<FixedImageType, MovingImageType>::New();
  myHelper->SetTransformType(this->m_TransformType);
  myHelper->SetFixedVolume(this->m_FixedVolume);
  myHelper->SetFixedVolume2(this->m_FixedVolume2);
  myHelper->SetMovingVolume(this->m_PreprocessedMovingVolume);
  myHelper->SetMovingVolume2(this->m_PreprocessedMovingVolume2);
  myHelper->SetHistogramMatch(this->m_HistogramMatch);
  myHelper->SetRemoveIntensityOutliers(this->m_RemoveIntensityOutliers);
  myHelper->SetNumberOfMatchPoints(this->m_NumberOfMatchPoints);
  myHelper->SetFixedBinaryVolume(this->m_FixedBinaryVolume);
  myHelper->SetMovingBinaryVolume(this->m_MovingBinaryVolume);
  myHelper->SetOutputFixedVolumeROI(this->m_OutputFixedVolumeROI);
  myHelper->SetOutputMovingVolumeROI(this->m_OutputMovingVolumeROI);
  myHelper->SetSamplingPercentage(this->m_SamplingPercentage);
  myHelper->SetSamplingStrategy(this->m_SamplingStrategy);
  myHelper->SetNumberOfHistogramBins(this->m_NumberOfHistogramBins);
  myHelper->SetNumberOfIterations(this->m_NumberOfIterations);
  myHelper->SetMaximumStepLength(this->m_MaximumStepLength);
  myHelper->SetMinimumStepLength(this->m_MinimumStepLength);
  myHelper->SetRelaxationFactor(this->m_RelaxationFactor);
  myHelper->SetTranslationScale(this->m_TranslationScale);
  myHelper->SetReproportionScale(this->m_ReproportionScale);
  myHelper->SetSkewScale(this->m_SkewScale);
  myHelper->SetBackgroundFillValue(this->m_BackgroundFillValue);
  myHelper->SetInitializeTransformMode(this->m_InitializeTransformMode);
  myHelper->SetMaskInferiorCutOffFromCenter(this->m_MaskInferiorCutOffFromCenter);
  myHelper->SetCurrentGenericTransform(this->m_CurrentGenericTransform);
  myHelper->SetRestoreState(this->m_RestoreState);
  myHelper->SetSplineGridSize(this->m_SplineGridSize);
  myHelper->SetCostFunctionConvergenceFactor(this->m_CostFunctionConvergenceFactor);
  myHelper->SetProjectedGradientTolerance(this->m_ProjectedGradientTolerance);
  myHelper->SetMaxBSplineDisplacement(this->m_MaxBSplineDisplacement);
  myHelper->SetDisplayDeformedImage(this->m_DisplayDeformedImage);
  myHelper->SetPromptUserAfterDisplay(this->m_PromptUserAfterDisplay);
  myHelper->SetDebugLevel(this->m_DebugLevel);
  myHelper->SetCostMetricObject(multiMetric);
  myHelper->SetUseROIBSpline(this->m_UseROIBSpline);
  myHelper->SetInitializeRegistrationByCurrentGenericTransform(this->m_InitializeRegistrationByCurrentGenericTransform);
  myHelper->SetMaximumNumberOfEvaluations(this->m_MaximumNumberOfEvaluations);
  myHelper->SetMaximumNumberOfCorrections(this->m_MaximumNumberOfCorrections);
  myHelper->SetSyNMetricType(this->m_CostMetricName);
  myHelper->SetSaveState(this->m_SaveState);
  myHelper->SetSyNFull(this->m_SyNFull);
  if (this->m_DebugLevel > 7)
  {
    this->PrintCommandLine(true, "BF");
  }
  this->m_Helper = static_cast<itk::Object *>(myHelper.GetPointer());
}

template <typename TLocalCostMetric>
void
BRAINSFitHelper::RunRegistration()
{
  const typename HelperType::Pointer myHelper(dynamic_cast<HelperType *>(this->m_Helper.GetPointer()));
  if (myHelper.IsNull())
  {
    std::cout << "ERROR:  Invalid BRAINSFitHelper conversion" << __FILE__ << " " << __LINE__ << std::endl;
  }
  myHelper->Update();
  this->m_CurrentGenericTransform = myHelper->GetCurrentGenericTransform();
  this->m_ActualNumberOfIterations = myHelper->GetActualNumberOfIterations();
  this->m_PermittedNumberOfIterations = myHelper->GetPermittedNumberOfIterations();
  /*
    // Find the final metric value based on the used metric type and returned output transform
    // using GenericMetricType = typename HelperType::MetricType;
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
} // end namespace itk

#endif // __BRAINSFITHELPER__
