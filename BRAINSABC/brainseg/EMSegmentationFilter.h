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
//
//
// //////////////////////////////////////////////////////////////////////////////
//
// Atlas based segmentation using the Expectation Maximization algorithm
//
// Designed for 3D MRI originally inspired by prastawa@cs.unc.edu 3/2004
//
// Van Leemput K, Maes F, Vandermeulen D, Suetens P. Automated model based
// tissue classification of MR images of the brain. IEEE TMI 1999; 18:897-908.
//
//
//
// //////////////////////////////////////////////////////////////////////////////
#ifndef __EMSegmentationFilter_h
#define __EMSegmentationFilter_h

#include "GeneratePurePlugMask.h"
#include <map>
#include <list>
class AtlasDefinition;

/**
 * \class EMSegmentationFilter
 */
template <typename TInputImage, typename TProbabilityImage>
class EMSegmentationFilter : public itk::ProcessObject
{
public:
  // Standard class type alias
  using Self = EMSegmentationFilter;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  // Method for creation through the object factory
  itkNewMacro(Self);

  // The dimension of the image we're working with
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  using CoordinateRepType = double;

  using RangeDBType = orderedmap<std::string, orderedmap<std::string, AtlasDefinition::BoundsType>>;
  // Image types
  using InputImageType = TInputImage;
  using InputImagePointer = typename TInputImage::Pointer;
  using InputImageIndexType = typename TInputImage::IndexType;
  using InputImageOffsetType = typename TInputImage::OffsetType;
  using InputImagePixelType = typename TInputImage::PixelType;
  using InputImageRegionType = typename TInputImage::RegionType;
  using InputImageSizeType = typename TInputImage::SizeType;
  using InputImageSpacingType = typename TInputImage::SpacingType;

  using BackgroundValueVector = std::vector<InputImagePixelType>;

  using InputImageVector = std::vector<InputImagePointer>;
  using MapOfInputImageVectors = orderedmap<std::string, InputImageVector>;

  using ByteImagePointer = typename ByteImageType::Pointer;
  using ByteImageIndexType = typename ByteImageType::IndexType;
  using ByteImageOffsetType = typename ByteImageType::OffsetType;
  using ByteImagePixelType = typename ByteImageType::PixelType;
  using ByteImageRegionType = typename ByteImageType::RegionType;
  using ByteImageSizeType = typename ByteImageType::SizeType;
  using ByteImageVectorType = std::vector<ByteImagePointer>;

  using ShortImageType = itk::Image<short, Self::ImageDimension>;
  using ShortImagePointer = typename ShortImageType::Pointer;
  using ShortImageIndexType = typename ShortImageType::IndexType;
  using ShortImageOffsetType = typename ShortImageType::OffsetType;
  using ShortImagePixelType = typename ShortImageType::PixelType;
  using ShortImageRegionType = typename ShortImageType::RegionType;
  using ShortImageSizeType = typename ShortImageType::SizeType;

  using ProbabilityImageType = TProbabilityImage;
  using ProbabilityImagePointer = typename ProbabilityImageType::Pointer;
  using ProbabilityImageIndexType = typename ProbabilityImageType::IndexType;
  using ProbabilityImageOffsetType = typename ProbabilityImageType::OffsetType;
  using ProbabilityImagePixelType = typename ProbabilityImageType::PixelType;
  using ProbabilityImageRegionType = typename ProbabilityImageType::RegionType;
  using ProbabilityImageSizeType = typename ProbabilityImageType::SizeType;
  using ProbabilityImageSpacingType = typename ProbabilityImageType::SpacingType;
  using ProbabilityImageVectorType = std::vector<ProbabilityImagePointer>;

  using VectorType = vnl_vector<FloatingPrecision>;
  using IntVectorType = vnl_vector<unsigned int>;
  using BoolVectorType = std::vector<bool>;
  using MatrixType = vnl_matrix<FloatingPrecision>;
  using MatrixInverseType = vnl_matrix_inverse<FloatingPrecision>;

  using BSplineTransformType = itk::BSplineTransform<CoordinateRepType, 3, 3>;
  using BSplineTransformPointer = BSplineTransformType::Pointer;

  using GenericTransformType = itk::Transform<double, 3, 3>;

  using MeasurementVectorType = itk::Array<FloatingPrecision>;
  using SampleType = itk::Statistics::ListSample<MeasurementVectorType>;

  using InputImageNNInterpolationType = itk::NearestNeighborInterpolateImageFunction<InputImageType, double>;
  using MaskNNInterpolationType = itk::NearestNeighborInterpolateImageFunction<ByteImageType, double>;

  using InputImageInterpolatorVector = std::vector<typename InputImageNNInterpolationType::Pointer>;
  using MapOfInputImageInterpolatorVectors = orderedmap<std::string, InputImageInterpolatorVector>;

  itkSetMacro(UseKNN, bool);
  itkGetMacro(UseKNN, bool);

  itkSetMacro(UsePurePlugs, bool);
  itkGetMacro(UsePurePlugs, bool);

  itkSetMacro(PurePlugsThreshold, float);
  itkGetMacro(PurePlugsThreshold, float);

  void
  SetNumberOfSubSamplesInEachPlugArea(unsigned int nx, unsigned int ny, unsigned int nz)
  {
    m_NumberOfSubSamplesInEachPlugArea[0] = nx;
    m_NumberOfSubSamplesInEachPlugArea[1] = ny;
    m_NumberOfSubSamplesInEachPlugArea[2] = nz;
    this->Modified();
  }

  // Set/Get the maximum polynomial degree of the bias field estimate
  itkSetMacro(MaxBiasDegree, unsigned int);
  itkGetMacro(MaxBiasDegree, unsigned int);
  //
  // Set/Get the Debugging level for filter verboseness
  itkSetMacro(DebugLevel, unsigned int);
  itkGetMacro(DebugLevel, unsigned int);

  itkSetMacro(BiasLikelihoodTolerance, FloatingPrecision);
  itkGetMacro(BiasLikelihoodTolerance, FloatingPrecision);

  itkSetMacro(OutputDebugDir, std::string);
  itkGetMacro(OutputDebugDir, std::string);

  itkSetMacro(LikelihoodTolerance, FloatingPrecision);
  itkGetMacro(LikelihoodTolerance, FloatingPrecision);

  itkSetMacro(MaximumIterations, unsigned int);
  itkGetMacro(MaximumIterations, unsigned int);

  itkSetMacro(SampleSpacing, FloatingPrecision);
  itkGetMacro(SampleSpacing, FloatingPrecision);

  void
  SetInputImages(const MapOfInputImageVectors newInputImages);

  void
  SetRawInputImages(const MapOfInputImageVectors newInputImages);

  void
  SetOriginalAtlasImages(const MapOfInputImageVectors newAtlasImages);

  //
  //  itkGetMacro(WarpedAtlasImages,std::vector<InputImagePointer>);
  MapOfInputImageVectors
  GenerateWarpedAtlasImages();

  itkSetMacro(TemplateBrainMask, ByteImagePointer);
  itkGetMacro(TemplateBrainMask, ByteImagePointer);

  void
  SetPriors(ProbabilityImageVectorType priors);

  void
  SetPriorWeights(VectorType w);

  IntVectorType
  GetPriorLabelCodeVector() const
  {
    return m_PriorLabelCodeVector;
  }

  void
  SetPriorLabelCodeVector(IntVectorType ng);

  BoolVectorType
  GetPriorUseForBiasVector() const
  {
    return m_PriorUseForBiasVector;
  }

  void
  SetPriorUseForBiasVector(const BoolVectorType & ng);

  BoolVectorType
  GetPriorIsForegroundPriorVector() const
  {
    return m_PriorIsForegroundPriorVector;
  }

  void
  SetPriorIsForegroundPriorVector(const BoolVectorType & ng);

  void
  SetPriorNames(const std::vector<std::string> & newPriorNames)
  {
    this->m_PriorNames = newPriorNames;
  }

  std::vector<std::string>
  GetPriorNames() const
  {
    return this->m_PriorNames;
  }

  ByteImagePointer
  GetOutput();

  ByteImagePointer
  GetCleanedOutput();

  ByteImagePointer
  GetThresholdedOutput();

  ProbabilityImageVectorType
  GetPosteriors();

  MapOfInputImageVectors
  GetCorrected();

  MapOfInputImageVectors
  GetRawCorrected();

  void
  Update() override;

  itkGetMacro(AtlasTransformType, std::string);
  itkSetMacro(AtlasTransformType, std::string);

  // Standard ITK style get/set macros for DoWarp
  itkGetMacro(UpdateTransformation, bool);
  itkSetMacro(UpdateTransformation, bool);
  itkBooleanMacro(UpdateTransformation);

  /*  For backwards compatibility */
  /*
  void WarpingOn() { this->SetDoWarp(true); }
  void WarpingOff() { this->SetDoWarp(false); }
  */

  itkGetMacro(TemplateGenericTransform, GenericTransformType::Pointer);
  itkSetMacro(TemplateGenericTransform, GenericTransformType::Pointer);

  /**
   * Set the index of the background region
   */
  itkGetMacro(AirIndex, LOOPITERTYPE);
  itkSetMacro(AirIndex, LOOPITERTYPE);

  void
  SetWarpGrid(unsigned int gx, unsigned int gy, unsigned int gz)
  {
    m_WarpGrid[0] = gx;
    m_WarpGrid[1] = gy;
    m_WarpGrid[2] = gz;
    this->Modified();
  }

  void
  SetTissueTypeThresholdMapsRange(const RangeDBType & newRangeDB)
  {
    this->m_TissueTypeThresholdMapsRange = newRangeDB;
    this->Modified();
  }

protected:
  EMSegmentationFilter();
  ~EMSegmentationFilter() override;

  void
  CheckInput();

  FloatingPrecision
  ComputeLogLikelihood() const;

  void
  EMLoop();

  void
  UpdateTransformation(const unsigned int CurrentEMIteration);

  ByteImageVectorType
  UpdateIntensityBasedClippingOfPriors(const unsigned int                 CurrentEMIteration,
                                       const MapOfInputImageVectors &     intensityList,
                                       const ProbabilityImageVectorType & WarpedPriorsList,
                                       ByteImagePointer &                 ForegroundBrainRegion);

  ByteImageVectorType
  ForceToOne(ProbabilityImageVectorType & WarpedPriorsList);

private:
  void
  WritePartitionTable(const unsigned int CurrentEMIteration) const;

  void
  WriteDebugLabels(const unsigned int CurrentEMIteration) const;

  void
  WriteDebugHeadRegion(const unsigned int CurrentEMIteration) const;

  void
  WriteDebugPosteriors(const unsigned int                 ComputeIterationID,
                       const std::string                  ClassifierID,
                       const ProbabilityImageVectorType & Posteriors) const;

  void
  WriteDebugBlendClippedPriors(const unsigned int CurrentEMIteration) const;

  void
  WriteDebugWarpedAtlasPriors(const unsigned int CurrentEMIteration) const;

  void
  WriteDebugWarpedAtlasImages(const unsigned int CurrentEMIteration) const;

  void
  WriteDebugForegroundMask(const ByteImageType::Pointer & currForgroundBrainMask,
                           const unsigned int             CurrentEMIteration) const;

  void
  WriteDebugCorrectedImages(const MapOfInputImageVectors & correctImageList,
                            const unsigned int             CurrentEMIteration) const;

  unsigned int
  ComputePriorLookupTable();

  void
  InitializePosteriors();

  void
  kNNCore(SampleType *                          trainSampleSet,
          const vnl_vector<FloatingPrecision> & labelVector,
          const vnl_matrix<FloatingPrecision> & testMatrix,
          vnl_matrix<FloatingPrecision> &       liklihoodMatrix,
          unsigned int                          K);

  typename TProbabilityImage::Pointer
  assignVectorToImage(const typename TProbabilityImage::Pointer prior, const vnl_vector<FloatingPrecision> & vector);

  std::vector<typename TProbabilityImage::Pointer>
  ComputekNNPosteriors(const ProbabilityImageVectorType & Priors,
                       const MapOfInputImageVectors &     intensityImages,
                       ByteImagePointer &                 labelsImage,
                       const IntVectorType &              labelClasses,
                       const std::vector<bool> &          priorIsForegroundPriorVector);

  typename TProbabilityImage::Pointer
  ComputeOnePosterior(const FloatingPrecision                   priorScale,
                      const typename TProbabilityImage::Pointer prior,
                      const vnl_matrix<FloatingPrecision>       currCovariance,
                      typename RegionStats::MeanMapType &       currMeans,
                      const MapOfInputImageVectors &            intensityImages);

  std::vector<typename TProbabilityImage::Pointer>
  ComputeEMPosteriors(const std::vector<typename TProbabilityImage::Pointer> & Priors,
                      const vnl_vector<FloatingPrecision> &                    PriorWeights,
                      const MapOfInputImageVectors &                           IntensityImages,
                      std::vector<RegionStats> &                               ListOfClassStatistics);

  std::vector<typename TProbabilityImage::Pointer>
  ComputePosteriors(const std::vector<typename TProbabilityImage::Pointer> & Priors,
                    const vnl_vector<FloatingPrecision> &                    PriorWeights,
                    const MapOfInputImageVectors &                           IntensityImages,
                    std::vector<RegionStats> &                               ListOfClassStatistics,
                    const IntVectorType &                                    priorLabelCodeVector,
                    std::vector<bool> &                                      priorIsForegroundPriorVector,
                    typename ByteImageType::Pointer &                        nonAirRegion,
                    const unsigned int                                       IterationID);

  std::vector<RegionStats>
  ComputeDistributions(const ByteImageVectorType &        SubjectCandidateRegions,
                       const ProbabilityImageVectorType & probAllDistributions);

  void
  BlendPosteriorsAndPriors(const double                       blendPosteriorPercentage,
                           const ProbabilityImageVectorType & ProbList1,
                           const ProbabilityImageVectorType & ProbList2,
                           ProbabilityImageVectorType &       ReturnBlendedProbList);

  void
  CheckLoopAgainstFilterOutput(ByteImagePointer & loopImg, ByteImagePointer & filterImg);

  ProbabilityImageVectorType
  WarpImageList(ProbabilityImageVectorType &        originalList,
                const typename TInputImage::Pointer referenceOutput,
                const BackgroundValueVector &       backgroundValues,
                const GenericTransformType::Pointer warpTransform);
  MapOfInputImageVectors
  WarpImageList(MapOfInputImageVectors &            originalList,
                const InputImagePointer             referenceOutput,
                const GenericTransformType::Pointer warpTransform);
  // External Templates to improve compilation times.
  MapOfInputImageVectors
  CorrectBias(const unsigned int                 degree,
              const unsigned int                 CurrentEMIteration,
              const ByteImageVectorType &        CandidateRegions,
              MapOfInputImageVectors &           inputImages,
              const ByteImageType::Pointer       currentBrainMask,
              const ByteImageType::Pointer       currentForegroundMask,
              const ProbabilityImageVectorType & probImages,
              const BoolVectorType &             probUseForBias,
              const int                          DebugLevel,
              const std::string &                OutputDebugDir);

  InputImagePointer
  GetFirstInputImage()
  {
    return GetMapVectorFirstElement(this->m_InputImages);
  }
  InputImagePointer
  GetFirstOriginalAtlasImage()
  {
    return GetMapVectorFirstElement(this->m_OriginalAtlasImages);
  }
  template <typename TMap>
  unsigned long
  MapofListsSize(const TMap & theMap)
  {
    unsigned long count = 0;
    for (typename TMap::iterator mapIt = theMap.begin(); mapIt != theMap.end(); ++mapIt)
    {
      for (typename TMap::mapped_type::iterator it = mapIt->second.begin(); it != mapIt->second.end(); ++it)
      {
        ++count;
      }
    }
    return count;
  }

  ProbabilityImageVectorType m_WarpedPriors;
  ProbabilityImageVectorType m_OriginalSpacePriors;
  ProbabilityImageVectorType m_Posteriors;

  std::string m_AtlasTransformType;

  // Variable set if the inputs are modified
  MapOfInputImageVectors m_InputImages;
  MapOfInputImageVectors m_RawInputImages;
  MapOfInputImageVectors m_CorrectedImages;
  MapOfInputImageVectors m_RawCorrectedImages;
  MapOfInputImageVectors m_OriginalAtlasImages;
  MapOfInputImageVectors m_WarpedAtlasImages;

  ByteImagePointer m_TemplateBrainMask;

  // final output
  ByteImagePointer m_DirtyLabels;
  ByteImagePointer m_CleanedLabels;
  ByteImagePointer m_ThresholdedLabels;
  ByteImagePointer m_DirtyThresholdedLabels;

  // exclude region from outside image space created by warping an all ones
  // image with zero default value, anded across images
  ByteImagePointer m_NonAirRegion;

  FloatingPrecision m_SampleSpacing;

  unsigned int      m_MaxBiasDegree;
  FloatingPrecision m_BiasLikelihoodTolerance;
  FloatingPrecision m_LikelihoodTolerance;
  unsigned int      m_MaximumIterations;

  VectorType m_PriorWeights;
  bool       m_PriorWeightsSet;

  IntVectorType m_PriorLabelCodeVector;
  bool          m_PriorLabelCodeVectorSet;

  BoolVectorType m_PriorUseForBiasVector;
  bool           m_PriorUseForBiasVectorSet;

  BoolVectorType m_PriorIsForegroundPriorVector;
  bool           m_PriorIsForegroundPriorVectorSet;

  BackgroundValueVector m_PriorsBackgroundValues;

  std::string m_OutputDebugDir;

  std::vector<RegionStats> m_ListOfClassStatistics;

  bool m_UseKNN;

  bool             m_UsePurePlugs;
  float            m_PurePlugsThreshold;
  unsigned int     m_NumberOfSubSamplesInEachPlugArea[3];
  ByteImagePointer m_PurePlugsMask;

  bool         m_UpdateTransformation;
  unsigned int m_DebugLevel;

  GenericTransformType::Pointer m_TemplateGenericTransform;
  unsigned int                  m_WarpGrid[3];

  FloatingPrecision m_WarpLikelihoodTolerance;
  bool              m_UpdateRequired;

  std::vector<std::string> m_PriorNames;
  std::vector<size_t>      m_ClassToPriorMapping;
  RangeDBType              m_TissueTypeThresholdMapsRange;
  LOOPITERTYPE             m_AirIndex;
};

#ifndef MU_MANUAL_INSTANTIATION
#  include "EMSegmentationFilter.hxx"
#endif

#endif
