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
#ifndef __BRAINSFitHelperTemplate_h
#define __BRAINSFitHelperTemplate_h

/**
 * \author Hans J. Johnson
 *
 * The intension of the BRIANSFitHelperTemplate is to provide a
 * class that can be used in other programs in a way that is very similar to
 * the command line version of the program from the SlicerExecutionModel
 * version of the BRAINSFitPrimary program.
 *
 * Almost all the command line options are available in this version, but
 * there is no need to read or write files to disk in order to use this class.
 */
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>

#include "itkIO.h"
#include "itkVector.h"
#include "itkMedianImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"

#include "itkCenteredVersorTransformInitializer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkDisplacementFieldTransform.h"

#include "itkFindCenterOfBrainFilter.h"

#include "itkBSplineTransform.h"
#include "itkBSplineTransformParametersAdaptor.h"
#include "itkBSplineTransformInitializer.h"

#include "itkScalableAffineTransform.h"

#ifdef USE_ANTS
#  include "BRAINSFitSyN.h"
#endif
#include "BRAINSFitUtils.h"

#include "BRAINSTypes.h"

#include "ConvertToRigidAffine.h"
#include "genericRegistrationHelper.h"
#include "ReadMask.h"
#include "BRAINSMacro.h"
#include "GenericTransformImage.h"
#include "BRAINSCommonLibWin32Header.h"


namespace itk
{
/** Method for verifying that the ordering of the transformTypes is consistent
 * with converting routines. */
extern void
ValidateTransformRankOrdering( const std::vector< std::string > & transformType );
} // namespace itk

namespace itk
{
template < typename FixedImageType, typename MovingImageType >
class BRAINSFitHelperTemplate : public Object
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( BRAINSFitHelperTemplate );

  /** Standard class type alias. */
  using Self = BRAINSFitHelperTemplate;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using RealType = double;

  using FixedImageConstPointer = typename FixedImageType::ConstPointer;
  using FixedImagePointer = typename FixedImageType::Pointer;

  using MovingImageConstPointer = typename MovingImageType::ConstPointer;
  using MovingImagePointer = typename MovingImageType::Pointer;

  /** Constants for the image dimensions */
  static constexpr unsigned int FixedImageDimension = FixedImageType::ImageDimension;
  static constexpr unsigned int MovingImageDimension = MovingImageType::ImageDimension;

  using MetricType = itk::ObjectToObjectMetricBaseTemplate< RealType >;
  using MultiMetricType =
    typename itk::ObjectToObjectMultiMetricv4< FixedImageDimension, MovingImageDimension, FixedImageType, RealType >;
  using ImageMetricType =
    typename itk::ImageToImageMetricv4< FixedImageType, MovingImageType, FixedImageType, RealType >;

  using CompositeTransformType = itk::CompositeTransform< RealType, MovingImageDimension >;
  using CompositeTransformPointer = typename CompositeTransformType::Pointer;
  using IdentityTransformType = IdentityTransform< RealType, MovingImageDimension >;

  using FixedBinaryVolumeType = SpatialObject< Self::FixedImageDimension >;
  using MovingBinaryVolumeType = SpatialObject< Self::MovingImageDimension >;
  using FixedBinaryVolumePointer = typename FixedBinaryVolumeType::Pointer;
  using MovingBinaryVolumePointer = typename MovingBinaryVolumeType::Pointer;

  using AffineRegistrationType = itk::ImageRegistrationMethodv4< FixedImageType, MovingImageType >;
  using TranslationTransformType = itk::TranslationTransform< RealType, MovingImageDimension >;
  using AffineTransformType = itk::AffineTransform< RealType, MovingImageDimension >;
  using ScalableAffineTransformType = itk::ScalableAffineTransform< RealType, MovingImageDimension >;
  using SamplingStrategyType = typename AffineRegistrationType::MetricSamplingStrategyType;

  using MatrixOffsetTransformBaseType = typename AffineTransformType::Superclass;
  using MatrixOffsetTransformBasePointer = typename MatrixOffsetTransformBaseType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( BRAINSFitHelperTemplate, ProcessObject );

  /** Set/Get the Fixed image. */
  itkSetObjectMacro( FixedVolume, FixedImageType );
  itkGetConstObjectMacro( FixedVolume, FixedImageType );

  /** Set/Get the second Fixed image for multi-modal SyN. */
  itkSetObjectMacro( FixedVolume2, FixedImageType );
  itkGetConstObjectMacro( FixedVolume2, FixedImageType );

  /** Set/Get the Moving image. */
  itkSetObjectMacro( MovingVolume, MovingImageType ) itkGetConstObjectMacro( MovingVolume, MovingImageType );

  /** Set/Get the second Moving image for multi-modal SyN. */
  itkSetObjectMacro( MovingVolume2, MovingImageType ) itkGetConstObjectMacro( MovingVolume2, MovingImageType );

  itkSetObjectMacro( FixedBinaryVolume, FixedBinaryVolumeType );
  itkGetModifiableObjectMacro( FixedBinaryVolume, FixedBinaryVolumeType );
  itkSetObjectMacro( MovingBinaryVolume, MovingBinaryVolumeType );
  itkGetModifiableObjectMacro( MovingBinaryVolume, MovingBinaryVolumeType );

  itkSetMacro( OutputFixedVolumeROI, std::string );
  itkGetConstMacro( OutputFixedVolumeROI, std::string );
  itkSetMacro( OutputMovingVolumeROI, std::string );
  itkGetConstMacro( OutputMovingVolumeROI, std::string );

  itkSetObjectMacro( CostMetricObject, MetricType );
  itkGetModifiableObjectMacro( CostMetricObject, MetricType );

  // TODO:  This should be converted to use the
  //       interpolation mechanisms from GenericTransform
  typedef enum
  {
    LINEAR_INTERP = 0,
    WINDOWSINC_INTERP = 1
  } InterpolationType;

  itkSetMacro( SamplingPercentage, RealType );
  itkGetConstMacro( SamplingPercentage, RealType );
  itkSetMacro( NumberOfHistogramBins, unsigned int );
  itkGetConstMacro( NumberOfHistogramBins, unsigned int );
  itkSetMacro( NumberOfMatchPoints, unsigned int );
  itkGetConstMacro( NumberOfMatchPoints, unsigned int );
  VECTORitkSetMacro( NumberOfIterations, std::vector< int > /**/ );
  VECTORitkSetMacro( MinimumStepLength, std::vector< double > );
  itkSetMacro( MaximumStepLength, double );
  itkGetConstMacro( MaximumStepLength, double );
  itkSetMacro( RelaxationFactor, double );
  itkGetConstMacro( RelaxationFactor, double );
  itkSetMacro( TranslationScale, double );
  itkGetConstMacro( TranslationScale, double );
  itkSetMacro( ReproportionScale, double );
  itkGetConstMacro( ReproportionScale, double );
  itkSetMacro( SkewScale, double );
  itkGetConstMacro( SkewScale, double );
  itkSetMacro( CostFunctionConvergenceFactor, double );
  itkGetConstMacro( CostFunctionConvergenceFactor, double );
  itkSetMacro( ProjectedGradientTolerance, double );
  itkGetConstMacro( ProjectedGradientTolerance, double );
  itkSetMacro( MaxBSplineDisplacement, double );
  itkGetConstMacro( MaxBSplineDisplacement, double );
  itkSetMacro( BackgroundFillValue, double );
  itkGetConstMacro( BackgroundFillValue, double );
  itkSetMacro( InitializeTransformMode, std::string );
  itkGetConstMacro( InitializeTransformMode, std::string );
  itkSetMacro( MaskInferiorCutOffFromCenter, double );
  itkGetConstMacro( MaskInferiorCutOffFromCenter, double );
  itkSetMacro( CurrentGenericTransform, CompositeTransformPointer );
  itkGetConstMacro( CurrentGenericTransform, CompositeTransformPointer );
  itkSetMacro( RestoreState, CompositeTransformPointer );
  itkGetConstMacro( RestoreState, CompositeTransformPointer );
  itkSetMacro( MaximumNumberOfEvaluations, int );
  itkGetConstMacro( MaximumNumberOfEvaluations, int );
  itkSetMacro( MaximumNumberOfCorrections, int );
  itkGetConstMacro( MaximumNumberOfCorrections, int );

  // cppcheck-suppress unusedFunction
  VECTORitkSetMacro( TransformType, std::vector< std::string > );
  // cppcheck-suppress unusedFunction
  VECTORitkSetMacro( SplineGridSize, std::vector< int > );

  itkGetConstMacro( ActualNumberOfIterations, unsigned int );
  itkGetConstMacro( PermittedNumberOfIterations, unsigned int );

  itkGetConstMacro( FinalMetricValue, double );
  /** Set/Get the Debugging level for filter verboseness */
  itkSetMacro( DebugLevel, unsigned int );
  itkGetConstMacro( DebugLevel, unsigned int );
  itkSetMacro( DisplayDeformedImage, bool );
  itkGetConstMacro( DisplayDeformedImage, bool );
  itkSetMacro( PromptUserAfterDisplay, bool );
  itkGetConstMacro( PromptUserAfterDisplay, bool );
  itkSetMacro( ObserveIterations, bool );
  itkGetConstMacro( ObserveIterations, bool );
  itkSetMacro( UseROIBSpline, bool );
  itkGetConstMacro( UseROIBSpline, bool );

  itkSetMacro( HistogramMatch, bool );
  itkGetConstMacro( HistogramMatch, bool );

  itkSetMacro( RemoveIntensityOutliers, bool );
  itkGetConstMacro( RemoveIntensityOutliers, bool );

  /** BRAINSFit uses ANTs to run only 3 levels of SyN registration unless SyNFull flag is false!
   * In that case only a single level SyN registration is run.
   */
  itkSetMacro( SyNFull, bool );
  itkGetConstMacro( SyNFull, bool );

  /** Method that initiates the registration. */
  void
  Update( void );

  itkSetMacro( SamplingStrategy, SamplingStrategyType );
  itkGetConstMacro( SamplingStrategy, SamplingStrategyType );

  itkSetMacro( InitializeRegistrationByCurrentGenericTransform, bool );

  itkSetMacro( SyNMetricType, std::string );
  itkGetConstMacro( SyNMetricType, std::string );

  itkSetMacro( SaveState, std::string );
  itkGetConstMacro( SaveState, std::string );

protected:
  BRAINSFitHelperTemplate();
  ~BRAINSFitHelperTemplate() override {}

  void
  PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  void
  GenerateData();

  template < typename TransformType >
  typename TransformType::Pointer
  CollapseLinearTransforms( const CompositeTransformType * compositeTransform );

  /** instantiate and call the Registration Helper */
  template < typename TransformType, typename OptimizerType, typename FitCommonCodeMetricType >
  void
  FitCommonCode( int numberOfIterations, double minimumStepLength,
                 typename CompositeTransformType::Pointer & initialITKTransform );

private:
  FixedImagePointer m_FixedVolume;
  FixedImagePointer m_FixedVolume2; // For multi-modal SyN

  MovingImagePointer m_MovingVolume;
  MovingImagePointer m_MovingVolume2; // For multi-modal SyN

  FixedBinaryVolumePointer  m_FixedBinaryVolume;
  MovingBinaryVolumePointer m_MovingBinaryVolume;
  std::string               m_OutputFixedVolumeROI;
  std::string               m_OutputMovingVolumeROI;

  RealType     m_SamplingPercentage;
  unsigned int m_NumberOfHistogramBins;
  bool         m_HistogramMatch;
  float        m_RemoveIntensityOutliers;
  unsigned int m_NumberOfMatchPoints;

  // TODO:  Would be better to have unsigned int
  std::vector< int >           m_NumberOfIterations;
  double                       m_MaximumStepLength;
  std::vector< double >        m_MinimumStepLength;
  double                       m_RelaxationFactor;
  double                       m_TranslationScale;
  double                       m_ReproportionScale;
  double                       m_SkewScale;
  double                       m_BackgroundFillValue;
  std::vector< std::string >   m_TransformType;
  std::string                  m_InitializeTransformMode;
  double                       m_MaskInferiorCutOffFromCenter;
  std::vector< int >           m_SplineGridSize;
  double                       m_CostFunctionConvergenceFactor;
  double                       m_ProjectedGradientTolerance;
  double                       m_MaxBSplineDisplacement;
  unsigned int                 m_ActualNumberOfIterations;
  unsigned int                 m_PermittedNumberOfIterations;
  unsigned int                 m_DebugLevel;
  CompositeTransformPointer    m_CurrentGenericTransform;
  CompositeTransformPointer    m_RestoreState;
  bool                         m_DisplayDeformedImage;
  bool                         m_PromptUserAfterDisplay;
  double                       m_FinalMetricValue;
  bool                         m_ObserveIterations;
  typename MetricType::Pointer m_CostMetricObject;
  bool                         m_UseROIBSpline;
  SamplingStrategyType         m_SamplingStrategy;
  bool                         m_InitializeRegistrationByCurrentGenericTransform;
  int                          m_MaximumNumberOfEvaluations;
  int                          m_MaximumNumberOfCorrections;
  std::string                  m_SyNMetricType;
  std::string                  m_SaveState;
  bool                         m_SyNFull;
  // DEBUG OPTION:
  int m_ForceMINumberOfThreads;
}; // end BRAINSFitHelperTemplate class
} // end namespace itk

#include "BRAINSFitHelperTemplate.hxx"

#endif // __BRAINSFITHELPER__
