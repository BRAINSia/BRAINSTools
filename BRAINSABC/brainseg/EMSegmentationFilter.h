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

#include "BRAINSABCUtilities.h"
#include <map>
#include <list>
class AtlasDefinition;

/**
 * \class EMSegmentationFilter
 */
template <class TInputImage, class TProbabilityImage>
class EMSegmentationFilter : public itk::ProcessObject
{
public:

  // Standard class typedefs
  typedef EMSegmentationFilter          Self;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  // Method for creation through the object factory
  itkNewMacro(Self);

  // The dimension of the image we're working with
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  typedef double CoordinateRepType;

  typedef typename std::map<std::string,  std::map<std::string, AtlasDefinition::BoundsType> > RangeDBType;
  // Image types
  typedef TInputImage                       InputImageType;
  typedef typename TInputImage::Pointer     InputImagePointer;
  typedef typename TInputImage::IndexType   InputImageIndexType;
  typedef typename TInputImage::OffsetType  InputImageOffsetType;
  typedef typename TInputImage::PixelType   InputImagePixelType;
  typedef typename TInputImage::RegionType  InputImageRegionType;
  typedef typename TInputImage::SizeType    InputImageSizeType;
  typedef typename TInputImage::SpacingType InputImageSpacingType;

  typedef std::vector<InputImagePixelType> BackgroundValueVector;

  typedef std::vector<InputImagePointer> InputImageVector;
  typedef std::map<std::string, InputImageVector> MapOfInputImageVectors;

  typedef typename ByteImageType::Pointer    ByteImagePointer;
  typedef typename ByteImageType::IndexType  ByteImageIndexType;
  typedef typename ByteImageType::OffsetType ByteImageOffsetType;
  typedef typename ByteImageType::PixelType  ByteImagePixelType;
  typedef typename ByteImageType::RegionType ByteImageRegionType;
  typedef typename ByteImageType::SizeType   ByteImageSizeType;
  typedef std::vector<ByteImagePointer>      ByteImageVectorType;

  typedef itk::Image<short, itkGetStaticConstMacro(ImageDimension)> ShortImageType;
  typedef typename ShortImageType::Pointer                          ShortImagePointer;
  typedef typename ShortImageType::IndexType                        ShortImageIndexType;
  typedef typename ShortImageType::OffsetType                       ShortImageOffsetType;
  typedef typename ShortImageType::PixelType                        ShortImagePixelType;
  typedef typename ShortImageType::RegionType                       ShortImageRegionType;
  typedef typename ShortImageType::SizeType                         ShortImageSizeType;

  typedef TProbabilityImage                          ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer     ProbabilityImagePointer;
  typedef typename ProbabilityImageType::IndexType   ProbabilityImageIndexType;
  typedef typename ProbabilityImageType::OffsetType  ProbabilityImageOffsetType;
  typedef typename ProbabilityImageType::PixelType   ProbabilityImagePixelType;
  typedef typename ProbabilityImageType::RegionType  ProbabilityImageRegionType;
  typedef typename ProbabilityImageType::SizeType    ProbabilityImageSizeType;
  typedef typename ProbabilityImageType::SpacingType ProbabilityImageSpacingType;
  typedef std::vector<ProbabilityImagePointer>       ProbabilityImageVectorType;

  typedef vnl_vector<FloatingPrecision>         VectorType;
  typedef vnl_vector<unsigned int>              IntVectorType;
  typedef std::vector<bool>                     BoolVectorType;
  typedef vnl_matrix<FloatingPrecision>         MatrixType;
  typedef vnl_matrix_inverse<FloatingPrecision> MatrixInverseType;

  typedef itk::BSplineTransform<CoordinateRepType, 3, 3 > BSplineTransformType;
  typedef BSplineTransformType::Pointer                   BSplineTransformPointer;

  typedef itk::Transform<double, 3, 3>  GenericTransformType;

  typedef std::vector< FloatingPrecision >                      MeasurementVectorType;
  typedef itk::Statistics::ListSample< MeasurementVectorType >  SampleType;

  itkSetMacro(UseKNN, bool);
  itkGetMacro(UseKNN, bool);

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

  void SetInputImages(const MapOfInputImageVectors newInputImages);

  void SetRawInputImages(const MapOfInputImageVectors newInputImages);

  void SetOriginalAtlasImages(const MapOfInputImageVectors newTemplateImages);

  //
  //  itkGetMacro(WarpedAtlasImages,std::vector<InputImagePointer>);
  MapOfInputImageVectors GenerateWarpedAtlasImages(void);

  itkSetMacro(TemplateBrainMask, ByteImagePointer);
  itkGetMacro(TemplateBrainMask, ByteImagePointer);

  void SetPriors(ProbabilityImageVectorType probs);

  void SetPriorWeights(VectorType w);

  IntVectorType GetPriorLabelCodeVector( void ) const
    {
      return m_PriorLabelCodeVector;
    }

  void SetPriorLabelCodeVector(IntVectorType n);

  BoolVectorType GetPriorUseForBiasVector( void ) const
    {
      return m_PriorUseForBiasVector;
    }

  void SetPriorUseForBiasVector(const BoolVectorType& n);

  BoolVectorType GetPriorIsForegroundPriorVector( void ) const
    {
      return m_PriorIsForegroundPriorVector;
    }

  void SetPriorIsForegroundPriorVector(const BoolVectorType& n);

  void SetPriorNames(const std::vector<std::string> & newPriorNames)
    {
      this->m_PriorNames = newPriorNames;
    }

  std::vector<std::string> GetPriorNames(void) const
    {
      return this->m_PriorNames;
    }

  ByteImagePointer GetOutput(void);

  ByteImagePointer GetCleanedOutput(void);

  ByteImagePointer GetThresholdedOutput(void);

  ProbabilityImageVectorType GetPosteriors();

  MapOfInputImageVectors GetCorrected();

  MapOfInputImageVectors  GetRawCorrected();

  void Update() ITK_OVERRIDE;

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

  void SetWarpGrid(unsigned int gx, unsigned int gy, unsigned int gz)
    {
      m_WarpGrid[0] = gx;
      m_WarpGrid[1] = gy;
      m_WarpGrid[2] = gz;
      this->Modified();
    }

  void SetTissueTypeThresholdMapsRange(const RangeDBType & newRangeDB)
    {
      this->m_TissueTypeThresholdMapsRange = newRangeDB;
      this->Modified();
    }

protected:

  EMSegmentationFilter(void);
  ~EMSegmentationFilter(void);

  void CheckInput(void);

  FloatingPrecision ComputeLogLikelihood(void) const;

  void EMLoop(void);

  void UpdateTransformation(const unsigned int CurrentEMIteration);

  ByteImageVectorType
  UpdateIntensityBasedClippingOfPriors(const unsigned int CurrentEMIteration,
                                       const MapOfInputImageVectors  &intensityList,
                                       const ProbabilityImageVectorType &WarpedPriorsList,
                                       ByteImagePointer &NonAirRegion);

  ByteImageVectorType ForceToOne(const unsigned int CurrentEMIteration,
                                 // const MapOfInputImageVectors &intensityList,
                                 ProbabilityImageVectorType &WarpedPriorsList,
                                 typename ByteImageType::Pointer NonAirRegion);
private:

  void WritePartitionTable(const unsigned int CurrentEMIteration) const;

  void WriteDebugLabels(const unsigned int CurrentEMIteration) const;

  void WriteDebugHeadRegion(const unsigned int CurrentEMIteration) const;

  void WriteDebugPosteriors(const unsigned int CurrentEMIteration,
                            const std::string ClassifierID,
                            const ProbabilityImageVectorType & Posteriors) const;

  void WriteDebugBlendClippedPriors(const unsigned int CurrentEMIteration) const;

  void WriteDebugWarpedAtlasPriors(const unsigned int CurrentEMIteration) const;

  void WriteDebugWarpedAtlasImages(const unsigned int CurrentEMIteration) const;

  void WriteDebugForegroundMask(const ByteImageType::Pointer & currForegroundMask,
                                const unsigned int CurrentEMIteration) const;

  void WriteDebugCorrectedImages(const MapOfInputImageVectors & correctImageList,
                                 const unsigned int CurrentEMIteration ) const;

  unsigned int ComputePriorLookupTable(void);

  void InitializePosteriors(void);

  typename TInputImage::Pointer
  NormalizeInputIntensityImage(const typename TInputImage::Pointer inputImage);

  void
  kNNCore( SampleType * trainMatrix,
           const vnl_vector<FloatingPrecision> & labelVector,
           const vnl_matrix<FloatingPrecision> & testMatrix,
           vnl_matrix<FloatingPrecision> & liklihoodMatrix,
           unsigned int K );

  typename TProbabilityImage::Pointer
  assignVectorToImage(const typename TProbabilityImage::Pointer prior,
                      const vnl_vector<FloatingPrecision> & vector);

  std::vector<typename TProbabilityImage::Pointer>
  ComputekNNPosteriors(const ProbabilityImageVectorType & Priors,
                        const MapOfInputImageVectors & IntensityImages,
                        ByteImagePointer & CleanedLabels,
                        const IntVectorType & labelClasses,
                        const std::vector<bool> & priorIsForegroundPriorVector);

  typename TProbabilityImage::Pointer
  ComputeOnePosterior(const FloatingPrecision priorScale,
                      const typename TProbabilityImage::Pointer prior,
                      const vnl_matrix<FloatingPrecision> currCovariance,
                      typename RegionStats::MeanMapType &currMeans,
                      const MapOfInputImageVectors & intensityImages);

  std::vector<typename TProbabilityImage::Pointer>
  ComputeEMPosteriors(const std::vector<typename TProbabilityImage::Pointer> & Priors,
                      const vnl_vector<FloatingPrecision> & PriorWeights,
                      const MapOfInputImageVectors & IntensityImages,
                      std::vector<RegionStats> & ListOfClassStatistics);

  std::vector<typename TProbabilityImage::Pointer>
  ComputePosteriors(const std::vector<typename TProbabilityImage::Pointer> & Priors,
                    const vnl_vector<FloatingPrecision> & PriorWeights,
                    const MapOfInputImageVectors & IntensityImages,
                    std::vector<RegionStats> & ListOfClassStatistics,
                    const IntVectorType & priorLabelCodeVector,
                    std::vector<bool> & priorIsForegroundPriorVector,
                    typename ByteImageType::Pointer & nonAirRegion,
                    const unsigned int IterationID);

  std::vector<RegionStats> ComputeDistributions(const ByteImageVectorType &SubjectCandidateRegions,
                                                const ProbabilityImageVectorType &probAllDistributions);

  void BlendPosteriorsAndPriors(const double blendPosteriorPercentage,
                                const ProbabilityImageVectorType & ProbList1,
                                const ProbabilityImageVectorType & ProbList2,
                                ProbabilityImageVectorType & ReturnBlendedProbList);

  void CheckLoopAgainstFilterOutput(ByteImagePointer &loopImg, ByteImagePointer & filterImg);

  ProbabilityImageVectorType  WarpImageList(ProbabilityImageVectorType& originalList,
                                            const typename TInputImage::Pointer referenceOutput,
                                            const BackgroundValueVector & backgroundValues,
                                            const GenericTransformType::Pointer warpTransform);
  MapOfInputImageVectors WarpImageList(MapOfInputImageVectors &originalList,
                                     const InputImagePointer referenceOutput,
                                     const GenericTransformType::Pointer warpTransform);
   // External Templates to improve compilation times.
  MapOfInputImageVectors
  CorrectBias(const unsigned int degree,
              const unsigned int CurrentEMIteration,
              const ByteImageVectorType & CandidateRegions,
              MapOfInputImageVectors & inputImages,
              const ByteImageType::Pointer currentBrainMask,
              const ByteImageType::Pointer currentForegroundMask,
              const ProbabilityImageVectorType & probImages,
              const BoolVectorType & probUseForBias,
              const FloatingPrecision sampleSpacing,
              const int DebugLevel,
              const std::string& OutputDebugDir);

  InputImagePointer GetFirstInputImage()
    {
      return GetMapVectorFirstElement(this->m_InputImages);
    }
  InputImagePointer GetFirstOriginalAtlasImage()
    {
      return GetMapVectorFirstElement(this->m_OriginalAtlasImages);
    }
  template <typename TMap>
  unsigned long MapofListsSize(const TMap &theMap)
    {
      unsigned long count = 0;
      for(typename TMap::iterator mapIt = theMap.begin(); mapIt != theMap.end(); ++mapIt)
        {
        for(typename TMap::mapped_type::iterator it = mapIt->second.begin();
            it != mapIt->second.end(); ++it)
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

  ByteImagePointer           m_TemplateBrainMask;

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

  bool         m_UseKNN;

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
#include "EMSegmentationFilter.hxx"
#endif

#endif
