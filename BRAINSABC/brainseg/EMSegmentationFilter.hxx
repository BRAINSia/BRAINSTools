#ifndef __EMSegmentationFilter_hxx
#define __EMSegmentationFilter_hxx

#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <cmath>
#include <cstdlib>

#include "itkAddImageFilter.h"
#include "itkBSplineDownsampleImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBlendImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiModeHistogramThresholdBinaryImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkNumericTraits.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkResampleImageFilter.h"
// #include "itkMersenneTwisterRandomVariateGenerator.h"

#include "vnl/algo/vnl_determinant.h"
#include "vnl/vnl_math.h"

#include "itkBRAINSROIAutoImageFilter.h"
#include "BRAINSFitBSpline.h"
#include "BRAINSFitUtils.h"
#include "BRAINSFitSyN.h"

// #include "QHullMSTClusteringProcess.h"
#include "AtlasDefinition.h"
#include "EMSegmentationFilter.h"
#include "ExtractSingleLargestRegion.h"
#include "PrettyPrintTable.h"
#include "ComputeDistributions.h"

template <class TInputImage, class TProbabilityImage>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::EMSegmentationFilter()
{
  m_DirtyLabels = NULL;
  m_CleanedLabels = NULL;
  m_ThresholdedLabels = NULL;
  m_DirtyThresholdedLabels = NULL;

  m_SampleSpacing = 2.0;

  // Bias
  m_MaxBiasDegree = 4;
  m_BiasLikelihoodTolerance = 1e-2;
  // NOTE: warp tol needs to be <= bias tol
  m_WarpLikelihoodTolerance = 1e-3;

  // EM convergence parameters
  m_LikelihoodTolerance = 1e-5;
  m_MaximumIterations = 40;

  m_PriorWeights = VectorType(0);
  m_PriorWeightsSet = false;

  // m_PriorGaussianClusterCountVector = IntVectorType(0);
  // m_PriorGaussianClusterCountVectorSet=false;

  m_PriorLabelCodeVector = IntVectorType(0);
  m_PriorLabelCodeVectorSet = false;

  m_PriorUseForBiasVector = BoolVectorType(0);
  m_PriorUseForBiasVectorSet = false;

  m_PriorIsForegroundPriorVector = BoolVectorType(false);
  m_PriorIsForegroundPriorVectorSet = false;

  m_PriorsBackgroundValues.resize(0);

  m_WarpedPriors.clear();
  m_OriginalSpacePriors.clear();
  m_Posteriors.clear();

  m_InputImages.clear();
  m_RawInputImages.clear();
  m_CorrectedImages.clear();
  m_RawCorrectedImages.clear();
  m_InputVolumeTypes.clear();

  m_TemplateBrainMask = NULL;
  m_OriginalAtlasImages.clear();
  m_WarpedAtlasImages.clear();

  m_OutputDebugDir = "";
  // m_PriorLookupTable = IntVectorType(0);

  m_NonAirRegion = 0;

  m_AtlasTransformType = "invalid_TransformationTypeNotSet";

  m_UpdateTransformation = false;

  m_DebugLevel = 0;

  m_TemplateGenericTransform = NULL;

  m_WarpGrid[0] = 5;
  m_WarpGrid[1] = 5;
  m_WarpGrid[2] = 5;

  m_UpdateRequired = true;
  m_AirIndex = 1;
  this->m_PriorNames.clear();
  // m_ClassToPriorMapping.clear();
}

template <class TInputImage, class TProbabilityImage>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::~EMSegmentationFilter()
{
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::CheckInput()
{
  if( m_WarpedPriors.size() < 1 )
    {
    itkExceptionMacro(<< "Must have one or more class probabilities" << std::endl );
    }

  if( m_PriorWeightsSet == false )
    {
    itkExceptionMacro(<< "The PriorWeights were not set." << std::endl );
    }
  if( m_PriorLabelCodeVectorSet == false )
    {
    itkExceptionMacro(<< "The PriorLabelCodeVector was not set." << std::endl );
    }
  if( m_PriorUseForBiasVectorSet == false )
    {
    itkExceptionMacro(<< "The PriorUseForBiasVector was not set." << std::endl );
    }
  if( m_PriorIsForegroundPriorVectorSet == false )
    {
    itkExceptionMacro(<< "The PriorIsForegroundPriorVector was not set." << std::endl );
    }

  if( m_WarpedPriors.size() != m_PriorWeights.size() )
    {
    itkExceptionMacro(<< "The PriorWeights vector size must match the number of priors listed." << std::endl );
    }
  if( m_WarpedPriors.size() != m_PriorLabelCodeVector.size() )
    {
    itkExceptionMacro(<< "The PriorLabelCodeVector vector size must match the number of priors listed." << std::endl );
    }
  if( m_WarpedPriors.size() != m_PriorUseForBiasVector.size() )
    {
    itkExceptionMacro(<< "The PriorUseForBiasVector vector size must match the number of priors listed." << std::endl );
    }
  if( m_WarpedPriors.size() != m_PriorIsForegroundPriorVector.size() )
    {
    itkExceptionMacro(
      << "The PriorIsForegroundPriorVector vector size must match the number of priors listed." << std::endl );
    }

  if( m_MaximumIterations == 0 )
    {
    itkWarningMacro(<< "Maximum iterations set to zero" << std::endl );
    }

  if( m_InputImages.empty() )
    {
    itkExceptionMacro(<< "No input images" << std::endl );
    }

    {
    const InputImageSizeType size = m_InputImages[0]->GetLargestPossibleRegion().GetSize();
    for( unsigned i = 1; i < m_InputImages.size(); i++ )
      {
      if( m_InputImages[i]->GetImageDimension() != 3 )
        {
        itkExceptionMacro(<< "InputImage [" << i << "] has invalid dimension: only supports 3D images" << std::endl );
        }
      const InputImageSizeType isize = m_InputImages[i]->GetLargestPossibleRegion().GetSize();
      if( size != isize )
        {
        itkExceptionMacro(
          << "Image data [" << i << "] 3D size mismatch " << size << " != " << isize << "." << std::endl );
        }
      }
    for( unsigned i = 0; i < m_WarpedPriors.size(); i++ )
      {
      if( m_WarpedPriors[i]->GetImageDimension() != 3 )
        {
        itkExceptionMacro(<< "Warped Prior [" << i << "] has invalid dimension: only supports 3D images" << std::endl );
        }
      const ProbabilityImageSizeType psize = m_WarpedPriors[i]->GetLargestPossibleRegion().GetSize();
      if( size != psize )
        {
        itkExceptionMacro(
          << "Warped prior [" << i << "] and atlas data 3D size mismatch" << size << " != " << psize << "."
          << std::endl );
        }
      }
    }

    {
    const InputImageSizeType atlasSize = m_OriginalAtlasImages[0]->GetLargestPossibleRegion().GetSize();
    for( unsigned i = 1; i < m_OriginalAtlasImages.size(); i++ )
      {
      if( m_OriginalAtlasImages[i]->GetImageDimension() != 3 )
        {
        itkExceptionMacro(<< "Atlas Image [" << i << "] has invalid dimension: only supports 3D images" << std::endl );
        }
      const InputImageSizeType asize = m_OriginalAtlasImages[i]->GetLargestPossibleRegion().GetSize();
      if( atlasSize != asize )
        {
        itkExceptionMacro(
          << "Image data [" << i << "] 3D size mismatch " << atlasSize << " != " << asize << "." << std::endl );
        }
      }
    for( unsigned i = 0; i < m_OriginalSpacePriors.size(); i++ )
      {
      if( m_OriginalSpacePriors[i]->GetImageDimension() != 3 )
        {
        itkExceptionMacro(<< "Prior [" << i << "] has invalid dimension: only supports 3D images" << std::endl );
        }
      const ProbabilityImageSizeType psize = m_OriginalSpacePriors[i]->GetLargestPossibleRegion().GetSize();
      if( atlasSize != psize )
        {
        itkExceptionMacro(
          << "Normalized prior [" << i << "] and atlas 3D size mismatch" << atlasSize << " != " << psize << "."
          << std::endl );
        }
      }
    }
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetInputImages(const std::vector<InputImagePointer> & newInputImages)
{
  muLogMacro(<< "SetInputImages" << std::endl);

  if( newInputImages.size() == 0 )
    {
    itkExceptionMacro(<< "No input images" << std::endl );
    }

  m_InputImages = newInputImages;

  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetRawInputImages(const std::vector<InputImagePointer> & newInputImages)
{
  muLogMacro(<< "SetRawInputImages" << std::endl);

  if( newInputImages.size() == 0 )
    {
    itkExceptionMacro(<< "No input images" << std::endl );
    }

  m_RawInputImages = newInputImages;

  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetOriginalAtlasImages(const std::vector<InputImagePointer> & newAtlasImages)
{
  muLogMacro(<< "SetAtlasImages" << std::endl);

  if( newAtlasImages.size() == 0 )
    {
    itkExceptionMacro(<< "No template images" << std::endl );
    }
  m_OriginalAtlasImages = newAtlasImages;

  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
void EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugPosteriors(const unsigned int ComputeIterationID) const
{
  if( this->m_DebugLevel > 9 )
    {
    // write out posteriors
    const unsigned int numPosteriors = this->m_Posteriors.size();
    const unsigned int write_posteriors_level = ComputeIterationID; // DEBUG:
                                                                    //  This
                                                                    // code is
                                                                    // for
                                                                    // debugging
                                                                    // purposes
                                                                    // only;
    std::stringstream write_posteriors_level_stream("");
    write_posteriors_level_stream << write_posteriors_level;
    for( unsigned int iprob = 0; iprob < numPosteriors; iprob++ )
      {
      typedef itk::ImageFileWriter<TProbabilityImage> ProbabilityImageWriterType;
      typename ProbabilityImageWriterType::Pointer writer = ProbabilityImageWriterType::New();

      std::stringstream template_index_stream("");
      template_index_stream << iprob;
      const std::string fn = this->m_OutputDebugDir + "/POSTERIOR_INDEX_" + template_index_stream.str() + "_"
        + this->m_PriorNames[iprob] + "_LEVEL_" + write_posteriors_level_stream.str() + ".nii.gz";

      muLogMacro(<< "Writing posterior images... " << fn <<  std::endl);
      writer->SetInput(m_Posteriors[iprob]);
      writer->SetFileName(fn);
      writer->UseCompressionOn();
      writer->Update();
      }
    }
  return;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetPriors(std::vector<ProbabilityImagePointer> priors)
{
  muLogMacro(<< "Set and Normalize for segmentation." << std::endl);
  // Need to normalize priors before getting started.
  this->m_OriginalSpacePriors = priors;
  ZeroNegativeValuesInPlace<TProbabilityImage>(this->m_OriginalSpacePriors);
  std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
  NormalizeProbListInPlace<TProbabilityImage>(this->m_OriginalSpacePriors);
  std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
  this->m_OriginalSpacePriors = priors;
  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetPriorWeights(VectorType w)
{
  muLogMacro(<< "SetPriorWeights" << std::endl);

  if( w.size() != m_OriginalSpacePriors.size() )
    {
    itkExceptionMacro(<< "Number of prior weights invalid" << w.size() << " != " << m_OriginalSpacePriors.size() );
    }
  for( unsigned i = 0; i < w.size(); i++ )
    {
    if( w[i] == 0.0 )
      {
      itkExceptionMacro(<< "Prior weight " << i << " is zero" << std::endl );
      }
    }

  m_PriorWeights = w;
  m_PriorWeightsSet = true;
  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetPriorLabelCodeVector(IntVectorType ng)
{
  muLogMacro(<< "SetPriorLabelCodeVector" << std::endl );
  if( ng.size() == 0 )
    {
    itkExceptionMacro(<< "Number of clusters info invalid" << std::endl );
    }
  const unsigned int numPriors = m_WarpedPriors.size();
  for( unsigned int i = 0; i < numPriors; i++ )
    {
    if( ng[i] == 0 )
      {
      itkExceptionMacro(<< "PriorLabelCode" << i << " is zero" << std::endl );
      }
    }
  m_PriorLabelCodeVector = ng;
  m_PriorLabelCodeVectorSet = true;
  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetPriorUseForBiasVector(const BoolVectorType& ng)
{
  muLogMacro(<< "SetPriorUseForBiasVector" << std::endl );
  if( ng.size() == 0 )
    {
    itkExceptionMacro(<< "Vector size for PriorUseForBiasVector info invalid" << std::endl );
    }
  const unsigned int numPriors = m_WarpedPriors.size();
  for( unsigned int i = 0; i < numPriors; i++ )
    {
    if( ng[i] != 0 && ng[i] != 1 )
      {
      itkExceptionMacro(<< "PriorUseForBiasVector" << i << " can only be 0 or 1" << std::endl );
      }
    }
  m_PriorUseForBiasVector = ng;
  m_PriorUseForBiasVectorSet = true;
  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SetPriorIsForegroundPriorVector(const BoolVectorType& ng)
{
  muLogMacro(<< "SetPriorIsForegroundPriorVector" << std::endl );
  if( ng.size() == 0 )
    {
    itkExceptionMacro(<< "Vector size for PriorIsForegroundPriorVector info invalid" << std::endl );
    }
  const unsigned int numPriors = m_WarpedPriors.size();
  for( unsigned int i = 0; i < numPriors; i++ )
    {
    if( ng[i] != 0 && ng[i] != 1 )
      {
      itkExceptionMacro(<< "PriorIsForegroundPriorVector" << i << " can only be 0 or 1" << std::endl );
      }
    }
  m_PriorIsForegroundPriorVector = ng;
  m_PriorIsForegroundPriorVectorSet = true;
  this->Modified();
  m_UpdateRequired = true;
}

template <class TInputImage, class TProbabilityImage>
typename EMSegmentationFilter<TInputImage, TProbabilityImage>::
ByteImagePointer
EMSegmentationFilter<TInputImage, TProbabilityImage>
::GetThresholdedOutput(void)
{
  // TODO:  This assumes that GetOutput was already called.  This should be made
  // more intelligent
  return m_DirtyThresholdedLabels;
}

template <class TInputImage, class TProbabilityImage>
typename EMSegmentationFilter<TInputImage, TProbabilityImage>::
ByteImagePointer
EMSegmentationFilter<TInputImage, TProbabilityImage>
::GetCleanedOutput(void)
{
  // TODO:  This assumes that GetOutput was already called.  This should be made
  // more intelligent
  return m_CleanedLabels;
}

template <class TInputImage, class TProbabilityImage>
typename EMSegmentationFilter<TInputImage, TProbabilityImage>::
ByteImagePointer
EMSegmentationFilter<TInputImage, TProbabilityImage>
::GetOutput(void)
{
  this->Update();
  return m_DirtyLabels;
}

template <class TInputImage, class TProbabilityImage>
std::vector<
  typename EMSegmentationFilter<TInputImage, TProbabilityImage>::
  ProbabilityImagePointer>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::GetPosteriors()
{
  return m_Posteriors;
}

template <class TInputImage, class TProbabilityImage>
std::vector<typename EMSegmentationFilter<TInputImage, TProbabilityImage>::InputImagePointer>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::GetCorrected()
{
  return m_CorrectedImages;
}

template <class TInputImage, class TProbabilityImage>
std::vector<typename EMSegmentationFilter<TInputImage, TProbabilityImage>::InputImagePointer>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::GetRawCorrected()
{
  return m_RawCorrectedImages;
}

template <class TInputImage, class TProbabilityImage>
void
CheckLoopAgainstFilterOutput
  (typename EMSegmentationFilter<TInputImage, TProbabilityImage>::ByteImageType::Pointer & loopImg,
  typename EMSegmentationFilter<TInputImage, TProbabilityImage>::ByteImageType::Pointer & filterImg)
{
  typedef typename EMSegmentationFilter<TInputImage, TProbabilityImage>::ByteImageType
    ByteImageType;

  typedef typename itk::ImageRegionConstIterator<ByteImageType> IterType;

  IterType     maskIter(loopImg, loopImg->GetLargestPossibleRegion() );
  IterType     dilIter(filterImg, filterImg->GetLargestPossibleRegion() );
  unsigned int count = 0;
  for( maskIter.Begin(), dilIter.Begin();
       !maskIter.IsAtEnd() && !dilIter.IsAtEnd(); ++maskIter, ++dilIter )
    {
    if( maskIter.Value() != dilIter.Value() )
      {
      std::cerr << "mask = " << static_cast<float>(maskIter.Value() )
                << " dilated = " << static_cast<float>(dilIter.Value() )
                << " at vIndex " << maskIter.GetIndex()
                << std::endl;
      count++;
      }
    }
  if( count == 0 )
    {
    muLogMacro( << "DEBUG:  Filter output same as after loop output!" << std::endl);
    }
}

namespace
{
template <class ImageType>
typename ImageType::Pointer
CopyImage(const typename ImageType::Pointer & input )
{
  typedef itk::ImageDuplicator<ImageType> ImageDupeType;
  typename ImageDupeType::Pointer MyDuplicator = ImageDupeType::New();
  MyDuplicator->SetInputImage(input);
  MyDuplicator->Update();
  return MyDuplicator->GetOutput();
}
}

template <class TInputImage, class TProbabilityImage>
std::vector<RegionStats> ComputeDistributions(
  const std::vector<typename ByteImageType::Pointer> & SubjectCandidateRegions,
  const std::vector<typename TProbabilityImage::Pointer> & probAllDistributions,
  const std::vector<typename TInputImage::Pointer> & intensityImages,
  const unsigned int DebugLevel)
{
  muLogMacro(<< "Computing Distributions" << std::endl );
  const std::vector<typename TProbabilityImage::Pointer> & probabilityMaps = probAllDistributions;

  std::vector<RegionStats> outputStats;
  typedef vnl_matrix<FloatingPrecision> MatrixType;

  CombinedComputeDistributions<TInputImage, TProbabilityImage, MatrixType>(SubjectCandidateRegions, intensityImages,
                                                                           probabilityMaps,
                                                                           outputStats,
                                                                           DebugLevel,
                                                                           false);

  return outputStats;
}

static double ComputeCovarianceDeterminant( const vnl_matrix<FloatingPrecision> & currCovariance)
{
  const FloatingPrecision detcov = vnl_determinant(currCovariance);

  if( detcov <= 0.0 )
    {
    itkGenericExceptionMacro(
      << "Determinant of covariance "
      << " is <= 0.0 (" << detcov << "), covariance matrix:" << std::endl
      << currCovariance
      <<
      "\n\n\n This is indicative of providing two images that are related only through a linear depenancy\n"
      <<
      "at least two images are so close in their ratio of values that a degenerate covariance matrix\n"
      << "would result, thus making an unstable calculation\n\n\n");
    }
  return detcov;
}

template <class TInputImage, class TProbabilityImage>
typename TProbabilityImage::Pointer
ComputeOnePosterior(
  const FloatingPrecision priorScale,
  const typename TProbabilityImage::Pointer prior,
  const vnl_matrix<FloatingPrecision> currCovariance,
  const vnl_vector<FloatingPrecision> currMeans,
  const std::vector<typename TInputImage::Pointer> & intensityImages
  )
{
  typedef vnl_matrix<FloatingPrecision>         MatrixType;
  typedef vnl_matrix_inverse<FloatingPrecision> MatrixInverseType;

  const unsigned          numChannels = currMeans.size();
  const FloatingPrecision detcov = ComputeCovarianceDeterminant(currCovariance);

  // Normalizing constant for the Gaussian
  const FloatingPrecision denom =
    vcl_pow(2 * vnl_math::pi, numChannels / 2.0) * vcl_sqrt(detcov) + vnl_math::eps;
  const FloatingPrecision invdenom = 1.0 / denom;
  CHECK_NAN(invdenom, __FILE__, __LINE__, "\n  denom:" << denom );
  const MatrixType invcov = MatrixInverseType(currCovariance);

  typename TProbabilityImage::Pointer post = TProbabilityImage::New();
  post->CopyInformation(prior);
  post->SetRegions(prior->GetLargestPossibleRegion() );
  post->Allocate();

  const typename TProbabilityImage::SizeType size = post->GetLargestPossibleRegion().GetSize();
    {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for
#endif
    for( LOOPITERTYPE kk = 0; kk < (LOOPITERTYPE)size[2]; kk++ )
      {
      for( LOOPITERTYPE jj = 0; jj < (LOOPITERTYPE)size[1]; jj++ )
        {
        for( LOOPITERTYPE ii = 0; ii < (LOOPITERTYPE)size[0]; ii++ )
          {
          const typename TProbabilityImage::IndexType currIndex = {{ii, jj, kk}};
          // At a minimum, every class has at least a 0.001% chance of being
          // true no matter what.
          // I realize that this small value makes the priors equal slightly
          // larger than 100%, but everything
          // is renormalized anyway, so it is not really that big of a deal as
          // long as the main priors for
          // the desired class is significantly higher than 1%.
          const typename TProbabilityImage::PixelType minPriorValue = 0.0;
          const typename TProbabilityImage::PixelType priorValue = (prior->GetPixel(currIndex) + minPriorValue);
            {
            MatrixType X(numChannels, 1);
            for( unsigned int ichan = 0; ichan < numChannels; ichan++ )
              {
              X(ichan, 0) =
                intensityImages[ichan]->GetPixel(currIndex) - currMeans[ichan];
              }

            const MatrixType  Y = invcov * X;
            FloatingPrecision mahalo = 0.0;
            for( unsigned int ichan = 0; ichan < numChannels; ichan++ )
              {
              const FloatingPrecision & currVal = X(ichan, 0) * Y(ichan, 0);
              CHECK_NAN(currVal, __FILE__, __LINE__, "\n  currIndex: " << currIndex
                                                                       << "\n  mahalo: " << mahalo
                                                                       << "\n  ichan: " << ichan
                                                                       << "\n  invcov: " << invcov
                                                                       << "\n  X:  " << X
                                                                       << "\n  Y:  " << Y );
              mahalo += currVal;
              }

            // Note:  This is the maximum likelyhood estimate as described in
            // formula at bottom of
            //       http://en.wikipedia.org/wiki/Maximum_likelihood_estimation
            const FloatingPrecision likelihood = vcl_exp(-0.5 * mahalo) * invdenom;

            const typename TProbabilityImage::PixelType currentPosterior =
              static_cast<typename TProbabilityImage::PixelType>( (priorScale * priorValue * likelihood) );
            CHECK_NAN(currentPosterior, __FILE__, __LINE__, "\n  currIndex: " << currIndex
                                                                              << "\n  priorScale: " << priorScale << "\n  priorValue: " << priorValue << "\n  likelihood: " << likelihood
                                                                              << "\n  mahalo: " << mahalo
                                                                              << "\n  invcov: " << invcov
                                                                              << "\n  X:  " << X
                                                                              << "\n  Y:  " << Y );
            post->SetPixel(currIndex, currentPosterior);
            }
          }
        }
      }
    }
  return post;
}

template <class TInputImage, class TProbabilityImage>
std::vector<typename TProbabilityImage::Pointer>
ComputePosteriors(const std::vector<typename TProbabilityImage::Pointer> & Priors,
                  const vnl_vector<FloatingPrecision> & PriorWeights,
                  const std::vector<typename TInputImage::Pointer> & IntensityImages,
                  std::vector<RegionStats> & ListOfClassStatistics)
{
  // Compute initial distribution parameters
  muLogMacro(<< "ComputePosteriors" << std::endl );
  itk::TimeProbe ComputePosteriorsTimer;
  ComputePosteriorsTimer.Start();

  const unsigned int numClasses = Priors.size();
  muLogMacro(<< "Computing posteriors at full resolution" << std::endl);

  std::vector<typename TProbabilityImage::Pointer> Posteriors;
  Posteriors.resize(numClasses);
  for( unsigned int iclass = 0; iclass < numClasses; iclass++ )
    {
    const FloatingPrecision priorScale = PriorWeights[iclass];
    CHECK_NAN(priorScale, __FILE__, __LINE__, "\n  iclass: " << iclass );

    Posteriors[iclass] = ComputeOnePosterior<TInputImage, TProbabilityImage>(
        priorScale,
        Priors[iclass],
        ListOfClassStatistics[iclass].m_Covariance,
        ListOfClassStatistics[iclass].m_Means,
        IntensityImages
        );
    } // end class loop

  ComputePosteriorsTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime =
    ComputePosteriorsTimer.GetTotal();
  muLogMacro(<< "Computing Posteriors took " << elapsedTime << " " << ComputePosteriorsTimer.GetUnit() << std::endl);
  return Posteriors;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugLabels(const unsigned int CurrentEMIteration) const
{
  if( this->m_DebugLevel > 6 )
    {
    // write out labels
    std::stringstream CurrentEMIteration_stream("");
    CurrentEMIteration_stream << CurrentEMIteration;
      {
      typedef itk::ImageFileWriter<ByteImageType> LabelImageWriterType;
      typename LabelImageWriterType::Pointer writer = LabelImageWriterType::New();

      const std::string fn = this->m_OutputDebugDir + "/LABELS_LEVEL_" + CurrentEMIteration_stream.str() + ".nii.gz";

      muLogMacro(<< "Writing label images... " << fn <<  std::endl);
      writer->SetInput(m_CleanedLabels);
      writer->SetFileName(fn);
      writer->UseCompressionOn();
      writer->Update();
      }
    }
  if( this->m_DebugLevel > 6 )
    {
    // write out labels
    std::stringstream CurrentEMIteration_stream("");
    CurrentEMIteration_stream << CurrentEMIteration;
      {
      typedef itk::ImageFileWriter<ByteImageType> LabelImageWriterType;
      typename LabelImageWriterType::Pointer writer = LabelImageWriterType::New();

      const std::string fn = this->m_OutputDebugDir + "/LABELSDIRTY_LEVEL_" + CurrentEMIteration_stream.str()
        + ".nii.gz";

      muLogMacro(<< "Writing label images... " << fn <<  std::endl);
      writer->SetInput(m_DirtyLabels);
      writer->SetFileName(fn);
      writer->UseCompressionOn();
      writer->Update();
      }
    }
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugCorrectedImages(const std::vector<typename TInputImage::Pointer> & correctImageList,
                            const unsigned int CurrentEMIteration ) const
{
  if( this->m_DebugLevel > 8 )
    { // DEBUG:  This code is for debugging purposes only;
    std::stringstream CurrentEMIteration_stream("");
    CurrentEMIteration_stream << CurrentEMIteration;
    for( unsigned int vIndex = 0; vIndex < correctImageList.size(); vIndex++ )
      {
        { // DEBUG:  This code is for debugging purposes only;
        typedef itk::ImageFileWriter<InputImageType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->UseCompressionOn();

        std::stringstream template_index_stream("");
        template_index_stream << vIndex;
        const std::string fn = this->m_OutputDebugDir + "/CORRECTED_INDEX_" + template_index_stream.str() + "_LEVEL_"
          + CurrentEMIteration_stream.str() + ".nii.gz";
        writer->SetInput(correctImageList[vIndex]);
        writer->SetFileName(fn.c_str() );
        writer->Update();
        muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
        }
      }
    }
}

template <class TInputImage, class TProbabilityImage>
FloatingPrecision
EMSegmentationFilter<TInputImage, TProbabilityImage>
::ComputeLogLikelihood() const
{
  const InputImageSizeType size = m_Posteriors[0]->GetLargestPossibleRegion().GetSize();
  const unsigned int       computeInitialNumClasses = m_Posteriors.size();
  FloatingPrecision        logLikelihood = 0.0;

    {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for reduction(+:logLikelihood)
#endif
    for( LOOPITERTYPE kk = 0; kk < (LOOPITERTYPE)size[2]; kk++ )
      {
      for( LOOPITERTYPE jj = 0; jj < (LOOPITERTYPE)size[1]; jj++ )
        {
        for( LOOPITERTYPE ii = 0; ii < (LOOPITERTYPE)size[0]; ii++ )
          {
          const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
          FloatingPrecision               tmp = 1e-20;
          for( unsigned int iclass = 0; iclass < computeInitialNumClasses; iclass++ )
            {
            if( this->m_PriorIsForegroundPriorVector[iclass] ) // We should
                                                               // probably only
                                                               // compute the
                                                               // foreground.
              {
              tmp += m_Posteriors[iclass]->GetPixel(currIndex);
              }
            }
          logLikelihood = logLikelihood + vcl_log(tmp);
          }
        }
      }
    }

  return logLikelihood;
}

/**
 * \param referenceImage is the image to be used for defining the tissueRegion of iterest.
 * \param safetyRegion is the amount to dilate so that there is not such a tight region.
 */
template <class TInputImage, class TByteImage>
typename TByteImage::Pointer
ComputeTissueRegion(const typename TInputImage::Pointer referenceImage, const unsigned int safetyRegion)
{
  typedef itk::BRAINSROIAutoImageFilter<TInputImage, TByteImage> ROIAutoType;
  typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
  ROIFilter->SetInput(referenceImage);
  ROIFilter->SetDilateSize(safetyRegion); // Create a very tight fitting tissue
                                          // region here.
  ROIFilter->Update();
  typename TByteImage::Pointer tissueRegion = ROIFilter->GetOutput();
  return tissueRegion;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugHeadRegion(const unsigned int CurrentEMIteration) const
{
  if( this->m_DebugLevel > 7 )
    {
    std::stringstream CurrentEMIteration_stream("");
    CurrentEMIteration_stream << CurrentEMIteration;
      { // DEBUG:  This code is for debugging purposes only;
      typedef itk::ImageFileWriter<ByteImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();

      const std::string fn = this->m_OutputDebugDir + "/HEAD_REGION_LEVEL_" + CurrentEMIteration_stream.str()
        + ".nii.gz";
      writer->SetInput( this->m_NonAirRegion );
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
      }
    }
}

template <class TInputImage>
std::vector<typename TInputImage::Pointer>
WarpImageList(const std::vector<typename TInputImage::Pointer> & originalList,
              const typename TInputImage::Pointer referenceOutput,
              const std::vector<typename TInputImage::PixelType> & backgroundValues,
              const GenericTransformType::Pointer warpTransform)
{
  if( originalList.size() != backgroundValues.size() )
    {
    itkGenericExceptionMacro(<< "ERROR:  originalList and backgroundValues arrays sizes do not match" << std::endl);
    }
  std::vector<typename TInputImage::Pointer> warpedList(originalList.size() );

  typedef itk::ResampleImageFilter<TInputImage, TInputImage> ResamplerType;
    {
    for( unsigned int vIndex = 0; vIndex < originalList.size(); vIndex++ )
      {
      typename ResamplerType::Pointer warper = ResamplerType::New();
      warper->SetInput(originalList[vIndex]);
      warper->SetTransform(warpTransform);

      // warper->SetInterpolator(linearInt); // Default is linear
      warper->SetOutputParametersFromImage(referenceOutput);
      warper->SetDefaultPixelValue(backgroundValues[vIndex]);
      warper->Update();
      warpedList[vIndex] = warper->GetOutput();
      }
    }
  return warpedList;
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugWarpedAtlasImages(const unsigned int CurrentEMIteration) const
{
  std::stringstream CurrentEMIteration_stream("");

  CurrentEMIteration_stream << CurrentEMIteration;
  if( this->m_DebugLevel > 9 )
    {
    for( unsigned int vIndex = 0; vIndex < this->m_WarpedAtlasImages.size(); vIndex++ )
      {
      typedef itk::ImageFileWriter<InputImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();

      std::stringstream template_index_stream("");
      template_index_stream << vIndex;
      const std::string fn = this->m_OutputDebugDir + "/WARPED_ATLAS_INDEX_" + template_index_stream.str()
        + "_LEVEL_" + CurrentEMIteration_stream.str() + ".nii.gz";
      writer->SetInput(m_WarpedAtlasImages[vIndex]);
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
      }
    }
}

template <class TInputImage, class TProbabilityImage>
std::vector<typename TInputImage::Pointer>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::GenerateWarpedAtlasImages(void)
{
  // TODO:  Need to make this cleaner.
  std::vector<typename TInputImage::PixelType> backgroundValues(m_OriginalAtlasImages.size() );
  std::fill(backgroundValues.begin(), backgroundValues.end(), 0);

  this->m_WarpedAtlasImages =
    WarpImageList<TProbabilityImage>(this->m_OriginalAtlasImages, this->m_InputImages[0], backgroundValues,
                                     this->m_TemplateGenericTransform);
  return m_WarpedAtlasImages;
}

template <class TInputImage, class TProbabilityImage>
std::vector<typename ByteImageType::Pointer>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::UpdateIntensityBasedClippingOfPriors(const unsigned int CurrentEMIteration,
                                       const std::vector<typename TInputImage::Pointer> & intensityList,
                                       std::vector<typename TProbabilityImage::Pointer> & WarpedPriorsList,
                                       typename ByteImageType::Pointer ForegroundBrainRegion)
{
  // #################################################################
  // #################################################################
  // #################################################################
  // #################################################################
  // #################################################################
  // #################################################################
  // For each intensityList, get it's type, and then create an "anded mask of
  // candidate regions"
  // using the table from BRAINSMultiModeHistogramThresholder.
  std::vector<typename ByteImageType::Pointer> subjectCandidateRegions;

  subjectCandidateRegions.resize(WarpedPriorsList.size() );
    { // StartValid Regions Section
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
    for( LOOPITERTYPE i = 0; i < (LOOPITERTYPE)WarpedPriorsList.size(); i++ )
      {
      typename ByteImageType::Pointer probThreshImage = NULL;
        {
        typedef itk::BinaryThresholdImageFilter<TProbabilityImage, ByteImageType> ProbThresholdType;
        typename ProbThresholdType::Pointer probThresh = ProbThresholdType::New();
        probThresh->SetInput(  WarpedPriorsList[i] );
        probThresh->SetInsideValue(1);
        probThresh->SetOutsideValue(0);
        probThresh->SetLowerThreshold(0.05);  // Hueristic: Need greater than 1 in 20
                                              // chance of being this structure
                                              // from the spatial probabilities
        probThresh->SetUpperThreshold( vcl_numeric_limits<typename TProbabilityImage::PixelType>::max() );
        // No upper limit needed, values
        // should be between 0 and 1
        probThresh->Update();
        probThreshImage = probThresh->GetOutput();
        if( this->m_DebugLevel > 9 )
          {
          std::stringstream CurrentEMIteration_stream("");
          CurrentEMIteration_stream << CurrentEMIteration;
          // Write the subject candidate regions

          std::ostringstream oss;
          oss << this->m_OutputDebugDir << "CANDIDIDATE_PROBTHRESH_" << this->m_PriorNames[i] << "_LEVEL_"
              << CurrentEMIteration_stream.str() << ".nii.gz" << std::ends;
          std::string fn = oss.str();
          muLogMacro( << "Writing Subject Candidate Region." << fn << std::endl );
          muLogMacro( << std::endl );

          typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
          typename ByteWriterType::Pointer writer = ByteWriterType::New();
          writer->SetInput(probThreshImage);
          writer->SetFileName(fn.c_str() );
          writer->UseCompressionOn();
          writer->Update();
          }
        }

      const unsigned int numberOfModes = intensityList.size();
      typedef typename itk::MultiModeHistogramThresholdBinaryImageFilter<InputImageType,
                                                                         ByteImageType> ThresholdRegionFinderType;
      typename ThresholdRegionFinderType::ThresholdArrayType QuantileLowerThreshold(numberOfModes);
      typename ThresholdRegionFinderType::ThresholdArrayType QuantileUpperThreshold(numberOfModes);
      typename ThresholdRegionFinderType::Pointer thresholdRegionFinder = ThresholdRegionFinderType::New();
      // TODO:  Need to define PortionMaskImage from deformed probspace
      thresholdRegionFinder->SetBinaryPortionImage(ForegroundBrainRegion);
      for( LOOPITERTYPE modeIndex = 0; modeIndex < (LOOPITERTYPE)numberOfModes; modeIndex++ )
        {
        thresholdRegionFinder->SetInput(modeIndex, intensityList[modeIndex]);
        const std::string imageType = this->m_InputVolumeTypes[modeIndex];
        const std::string priorType = this->m_PriorNames[i];
        if( m_TissueTypeThresholdMapsRange[priorType].find(imageType) ==
            m_TissueTypeThresholdMapsRange[priorType].end() )
          {
          muLogMacro(
            << "NOT FOUND:" << "[" << priorType << "," << imageType << "]: [" << 0.00 << "," << 1.00 << "]"
            <<  std::endl);
          QuantileLowerThreshold.SetElement(modeIndex, 0.00);
          QuantileUpperThreshold.SetElement(modeIndex, 1.00);
          }
        else
          {
          const float lower = m_TissueTypeThresholdMapsRange[priorType][imageType].GetLower();
          const float upper = m_TissueTypeThresholdMapsRange[priorType][imageType].GetUpper();
          muLogMacro( <<  "[" << priorType << "," << imageType << "]: [" << lower << "," << upper << "]" <<  std::endl);
          QuantileLowerThreshold.SetElement(modeIndex, lower);
          QuantileUpperThreshold.SetElement(modeIndex, upper);
          m_TissueTypeThresholdMapsRange[priorType][imageType].Print();
          }
        }
      // Assume upto (2*0.025)% of intensities are noise that corrupts the image
      // min/max values
      thresholdRegionFinder->SetLinearQuantileThreshold(0.025);
      thresholdRegionFinder->SetQuantileLowerThreshold(QuantileLowerThreshold);
      thresholdRegionFinder->SetQuantileUpperThreshold(QuantileUpperThreshold);
      // thresholdRegionFinder->SetInsideValue(1);
      // thresholdRegionFinder->SetOutsideValue(0);//Greatly reduce the value to
      // zero.
      thresholdRegionFinder->Update();
      if( this->m_DebugLevel > 8 )
        {
        std::stringstream CurrentEMIteration_stream("");
        CurrentEMIteration_stream << CurrentEMIteration;
        // Write the subject candidate regions

        std::ostringstream oss;
        oss << this->m_OutputDebugDir << "CANDIDIDATE_INTENSITY_REGION_" << this->m_PriorNames[i] << "_LEVEL_"
            << CurrentEMIteration_stream.str() << ".nii.gz" << std::ends;
        std::string fn = oss.str();
        muLogMacro( << "Writing Subject Candidate Region." << fn << std::endl );
        muLogMacro( << std::endl );

        typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
        typename ByteWriterType::Pointer writer = ByteWriterType::New();
        writer->SetInput(thresholdRegionFinder->GetOutput() );
        writer->SetFileName(fn.c_str() );
        writer->UseCompressionOn();
        writer->Update();
        }

      // Now multiply the warped priors by the subject candidate regions.
      typename itk::MultiplyImageFilter<ByteImageType, ByteImageType, ByteImageType>::Pointer multFilter =
        itk::MultiplyImageFilter<ByteImageType, ByteImageType, ByteImageType>::New();
      multFilter->SetInput1( probThreshImage );
      multFilter->SetInput2( thresholdRegionFinder->GetOutput() );
      multFilter->Update();
      subjectCandidateRegions[i] = multFilter->GetOutput();
      } // End loop over all warped priors

      { // Ensure that every candidate region has some value
      const unsigned int candiateVectorSize = subjectCandidateRegions.size();
      itk::ImageRegionIteratorWithIndex<ByteImageType>
             firstCandidateIter(subjectCandidateRegions[0], subjectCandidateRegions[0]->GetLargestPossibleRegion() );
      size_t AllZeroCounts = 0;
      while( !firstCandidateIter.IsAtEnd() )
        {
        const typename ByteImageType::IndexType myIndex = firstCandidateIter.GetIndex();
        bool AllPixelsAreZero = true;
//        unsigned int maxPriorRegionIndex = 0;
        typename TProbabilityImage::PixelType maxProbValue = WarpedPriorsList[0]->GetPixel(myIndex);
        for( unsigned int k = 0; ( k < candiateVectorSize ) && AllPixelsAreZero; k++ )
          {
          const typename ByteImageType::PixelType value = subjectCandidateRegions[k]->GetPixel(myIndex);
          if( value > 0 )
            {
            AllPixelsAreZero = false;
            }
          typename TProbabilityImage::PixelType currValue = WarpedPriorsList[k]->GetPixel(myIndex);
          if( currValue > maxProbValue )
            {
            maxProbValue = currValue;
//            maxPriorRegionIndex = k;
            }
          }
        if( AllPixelsAreZero ) // If all candidate regions are zero, then force
                               // to most likely background value.
          {
          AllZeroCounts++;
          for( unsigned int k = 0; k < candiateVectorSize; k++ )
            {
            if( this->m_PriorIsForegroundPriorVector[k] == false )
              {
              subjectCandidateRegions[k]->SetPixel(myIndex, 1);
              WarpedPriorsList[k]->SetPixel(myIndex, 0.05);
              }
            }
          }
        ++firstCandidateIter;
        }

      if( AllZeroCounts != 0 )
        {
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "Locations with no candidate regions specified!" << AllZeroCounts << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "^^^^^^^^^^^^^^^^" << std::endl;
        }
      }
    if( this->m_DebugLevel > 5 )
      {
      std::stringstream CurrentEMIteration_stream("");
      CurrentEMIteration_stream << CurrentEMIteration;
      for( unsigned int i = 0; i < subjectCandidateRegions.size(); i++ )
        {
        // Write the subject candidate regions

        std::ostringstream oss;
        oss << this->m_OutputDebugDir << "CANDIDIDATE_FINAL" << this->m_PriorNames[i] << "_LEVEL_"
            << CurrentEMIteration_stream.str() << ".nii.gz" << std::ends;
        std::string fn = oss.str();
        muLogMacro( << "Writing Subject Candidate Region." << fn << std::endl );
        muLogMacro( << std::endl );

        typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
        typename ByteWriterType::Pointer writer = ByteWriterType::New();
        writer->SetInput(subjectCandidateRegions[i]);
        writer->SetFileName(fn.c_str() );
        writer->UseCompressionOn();
        writer->Update();
        }
      }
    } // END Valid regions section

  return subjectCandidateRegions;
}

template <class TInputImage, class TProbabilityImage>
std::vector<typename ByteImageType::Pointer>
EMSegmentationFilter<TInputImage, TProbabilityImage>
::ForceToOne(const unsigned int, const std::vector<typename TInputImage::Pointer> &,
             std::vector<typename TProbabilityImage::Pointer> & WarpedPriorsList,
             typename ByteImageType::Pointer )
{
  // #################################################################
  // #################################################################
  // #################################################################
  // #################################################################
  // #################################################################
  // #################################################################
  // For each intensityList, get it's type, and then create an "anded mask of
  // candidate regions"
  // using the table from BRAINSMultiModeHistogramThresholder.
  std::vector<typename ByteImageType::Pointer> subjectCandidateRegions;

  subjectCandidateRegions.resize(WarpedPriorsList.size() );
    { // StartValid Regions Section
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
    for( LOOPITERTYPE i = 0; i < (LOOPITERTYPE)WarpedPriorsList.size(); i++ )
      {
      subjectCandidateRegions[i] = ByteImageType::New();
      subjectCandidateRegions[i]->CopyInformation(WarpedPriorsList[i].GetPointer() );
      subjectCandidateRegions[i]->SetRegions( WarpedPriorsList[i]->GetLargestPossibleRegion() );
      subjectCandidateRegions[i]->Allocate();
      subjectCandidateRegions[i]->FillBuffer(1);
      }
    } // END Valid regions section

  return subjectCandidateRegions;
}

// ReturnBlendedProbList can be the same as one of the inputs!
template <class TProbabilityImage>
void BlendPosteriorsAndPriors(const double blendPosteriorPercentage,
                              const std::vector<typename TProbabilityImage::Pointer> & ProbList1,
                              const std::vector<typename TProbabilityImage::Pointer> & ProbList2,
                              std::vector<typename TProbabilityImage::Pointer> & ReturnBlendedProbList)
{
  for( unsigned int k = 0; k < ProbList2.size(); k++ )
    {
    std::cout << "Start Blending Prior:" << k << std::endl;
    typename TProbabilityImage::Pointer multInputImage = ProbList2[k];
    // BLEND Posteriors and Priors Here:
    // It is important to keep the warped priors as at least a small component
    // of this part
    // of the algorithm, because otherwise single pixels that exactly match the
    // mean of the
    // NOT* regions will become part of those NOT* regions regardless of spatial
    // locations.
    // std::cout << "\n\nWarpedPriors[" << k << "] \n" << ProbList2[k] <<
    // std::endl;
    if(
      ( ProbList1.size() == ProbList2.size() )
      && ProbList1.size() > k
      && ProbList1[k].IsNotNull()
      && ( blendPosteriorPercentage > 0.01 ) // Need to blend at more than 1%,
                                             // else just skip it
      )
      {
      // Really we need to use a heirarchial approach to solving this problem.
      //  It is not sufficient to have these heuristics
      // break things apart artificially.
      std::cout << "\nBlending Priors with Posteriors with formula: " << (blendPosteriorPercentage)
                << "*Posterior + " << (1.0 - blendPosteriorPercentage) << "*Prior" << std::endl;
      typedef itk::BlendImageFilter<TProbabilityImage, TProbabilityImage> BlenderType;
      typename BlenderType::Pointer myBlender = BlenderType::New();
      myBlender->SetInput1(ProbList1[k]);
      myBlender->SetInput2(ProbList2[k]);
      myBlender->SetBlend1(blendPosteriorPercentage);
      myBlender->SetBlend2(1.0 - blendPosteriorPercentage);
      myBlender->Update();
      multInputImage = myBlender->GetOutput();
      }
    else
      {
      std::cout << "Not Blending Posteriors into Priors" << std::endl;
      multInputImage = ProbList2[k];
      }
#if 1
    ReturnBlendedProbList[k] = multInputImage;
#else
    // Now multiply the warped priors by the subject candidate regions.
    typename itk::MultiplyImageFilter<TProbabilityImage, ByteImageType, TProbabilityImage>::Pointer multFilter =
      itk::MultiplyImageFilter<TProbabilityImage, ByteImageType, TProbabilityImage>::New();
    multFilter->SetInput1(multInputImage);
    multFilter->SetInput2(candidateRegions[k]);
    multFilter->Update();
    ReturnBlendedProbList[k] = multFilter->GetOutput();
    std::cout << "Stop Blending Prior:" << k << std::endl;
#endif
    }
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugWarpedAtlasPriors(const unsigned int CurrentEMIteration) const
{
  std::stringstream CurrentEMIteration_stream("");

  CurrentEMIteration_stream << CurrentEMIteration;
  if( this->m_DebugLevel > 9 )
    {
    for( unsigned int vIndex = 0; vIndex < this->m_WarpedPriors.size(); vIndex++ )
      {
      typedef itk::ImageFileWriter<InputImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();

      std::stringstream template_index_stream("");
      template_index_stream << this->m_PriorNames[vIndex];
      const std::string fn = this->m_OutputDebugDir + "/WARPED_PRIOR_" + template_index_stream.str() + "_LEVEL_"
        + CurrentEMIteration_stream.str() + ".nii.gz";
      writer->SetInput(m_WarpedPriors[vIndex]);
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
      }
    }
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugBlendClippedPriors( const unsigned int CurrentEMIteration) const
{
  std::stringstream CurrentEMIteration_stream("");

  CurrentEMIteration_stream << CurrentEMIteration;
  if( this->m_DebugLevel > 9 )
    { // DEBUG:  This code is for debugging purposes only;
    for( unsigned int k = 0; k < m_WarpedPriors.size(); k++ )
      {
      typedef itk::ImageFileWriter<InputImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();

      std::stringstream prior_index_stream("");
      prior_index_stream << k;
      // const std::string fn = this->m_OutputDebugDir +
      //
      // "PRIOR_INDEX_"+prior_index_stream.str()+"_LEVEL_"+CurrentEMIteration_stream.str()+".nii.gz";
      const std::string fn = this->m_OutputDebugDir + "BLENDCLIPPED_PRIOR_INDEX_" + this->m_PriorNames[k] + "_LEVEL_"
        + CurrentEMIteration_stream.str() + ".nii.gz";
      writer->SetInput(m_WarpedPriors[k]);
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
      }
    }
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::UpdateTransformation(const unsigned int /*CurrentEMIteration*/)
{
  muLogMacro(<< "Updating Warping with transform type: " << m_AtlasTransformType  << std::endl );
  if( m_UpdateTransformation == false )
    {
    muLogMacro(
      << "WARNING: WARNING: WARNING: WARNING:  Doing warping even though it was turned off from the command line"
      << std::endl);
    }
  for( unsigned int vIndex = 0; vIndex < m_OriginalAtlasImages.size() && vIndex < m_CorrectedImages.size(); vIndex++ )
    {
    typedef itk::BRAINSFitHelper HelperType;
    HelperType::Pointer atlasToSubjectRegistrationHelper = HelperType::New();
    atlasToSubjectRegistrationHelper->SetNumberOfSamples(500000);
    atlasToSubjectRegistrationHelper->SetNumberOfHistogramBins(50);
    std::vector<int> numberOfIterations(1);
    numberOfIterations[0] = 1500;
    atlasToSubjectRegistrationHelper->SetNumberOfIterations(numberOfIterations);
    //  atlasToSubjectRegistrationHelper->SetMaximumStepLength(maximumStepSize);
    atlasToSubjectRegistrationHelper->SetTranslationScale(1000);
    atlasToSubjectRegistrationHelper->SetReproportionScale(1.0);
    atlasToSubjectRegistrationHelper->SetSkewScale(1.0);
    //
    //
    //
    // atlasToSubjectRegistrationHelper->SetBackgroundFillValue(backgroundFillValue);
    // NOT VALID When using initializeTransformMode
    //
    //
    // atlasToSubjectRegistrationHelper->SetCurrentGenericTransform(currentGenericTransform);
    //
    //
    //
    // atlasToSubjectRegistrationHelper->SetMaskInferiorCutOffFromCenter(maskInferiorCutOffFromCenter);
    //  atlasToSubjectRegistrationHelper->SetUseWindowedSinc(useWindowedSinc);
    // Register each intrasubject image mode to first image
    atlasToSubjectRegistrationHelper->SetFixedVolume(m_CorrectedImages[vIndex]);
    // Register all atlas images to first image
      {
      if( false /* m_AtlasTransformType == ID_TRANSFORM */ )
        {
        muLogMacro(<< "Registering (Identity) atlas to first image." << std::endl);
        // TODO: m_AtlasToSubjectTransform = MakeRigidIdentity();
        }
      else // continue;
        {
        std::string preprocessMovingString("");
        // const bool histogramMatch=true;//Setting histogram matching to true

        // Setting histogram matching to false because it appears to have been
        // causing problems for some images.
        const bool histogramMatch = false;
        if( histogramMatch )
          {
          typedef itk::HistogramMatchingImageFilter<InputImageType,
                                                    InputImageType> HistogramMatchingFilterType;
          typename HistogramMatchingFilterType::Pointer histogramfilter
            = HistogramMatchingFilterType::New();

          histogramfilter->SetInput( m_OriginalAtlasImages[vIndex]  );
          histogramfilter->SetReferenceImage( m_CorrectedImages[vIndex] );

          histogramfilter->SetNumberOfHistogramLevels( 50 );
          histogramfilter->SetNumberOfMatchPoints( 10 );
          histogramfilter->ThresholdAtMeanIntensityOn();
          histogramfilter->Update();
          typename InputImageType::Pointer equalizedMovingImage = histogramfilter->GetOutput();
          atlasToSubjectRegistrationHelper->SetMovingVolume(equalizedMovingImage);
          preprocessMovingString = "histogram equalized ";
          }
        else
          {
          atlasToSubjectRegistrationHelper->SetMovingVolume(m_OriginalAtlasImages[vIndex]);
          preprocessMovingString = "";
          }
          {
          typedef itk::BRAINSROIAutoImageFilter<InputImageType, itk::Image<unsigned char, 3> > ROIAutoType;
          typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
          ROIFilter->SetInput(m_OriginalAtlasImages[vIndex]);
          ROIFilter->SetDilateSize(10);
          ROIFilter->Update();
          atlasToSubjectRegistrationHelper->SetMovingBinaryVolume(ROIFilter->GetSpatialObjectROI() );
          }
          {
          typedef itk::BRAINSROIAutoImageFilter<InputImageType, itk::Image<unsigned char, 3> > ROIAutoType;
          typename ROIAutoType::Pointer  ROIFilter = ROIAutoType::New();
          ROIFilter->SetInput(m_CorrectedImages[vIndex]);
          ROIFilter->SetDilateSize(10);
          ROIFilter->Update();
          atlasToSubjectRegistrationHelper->SetFixedBinaryVolume(ROIFilter->GetSpatialObjectROI() );
          }
        // KENT TODO:  Need to expose the transform type to the command line
        // --AtlasTransformType [Rigid|Affine|BSpline|SyN], defaults to SyN.
        if( m_AtlasTransformType == "Rigid" )
          {
          muLogMacro(
            << "Registering (Rigid) " << preprocessMovingString << "atlas(" << vIndex << ") to template(" << vIndex
            << ") image." << std::endl);
          std::vector<double> minimumStepSize(1);
          minimumStepSize[0] = 0.005; // NOTE: 0.005 for between subject
                                      // registration is probably about the
                                      // limit.
          atlasToSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
          std::vector<std::string> transformType(1);
          transformType[0] = "Rigid";
          atlasToSubjectRegistrationHelper->SetTransformType(transformType);
          }
        else if( m_AtlasTransformType == "Affine" )
          {
          muLogMacro(
            << "Registering (Affine) " << preprocessMovingString << "atlas(" << vIndex << ") to template(" << vIndex
            << ") image." << std::endl);
          std::vector<double> minimumStepSize(1);
          minimumStepSize[0] = 0.0025;
          atlasToSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
          std::vector<std::string> transformType(1);
          transformType[0] = "Affine";
          atlasToSubjectRegistrationHelper->SetTransformType(transformType);
          }
        else if( m_AtlasTransformType == "BSpline" )
          {
          muLogMacro(
            << "Registering (BSpline) " << preprocessMovingString << "atlas(" << vIndex << ") to template(" << vIndex
            << ") image." << std::endl);
          std::vector<double> minimumStepSize(1);
          minimumStepSize[0] = 0.0025;
          atlasToSubjectRegistrationHelper->SetMinimumStepLength(minimumStepSize);
          std::vector<std::string> transformType(1);
          transformType[0] = "BSpline";
          atlasToSubjectRegistrationHelper->SetTransformType(transformType);
          std::vector<int> splineGridSize(3);
          splineGridSize[0] = m_WarpGrid[0];
          splineGridSize[1] = m_WarpGrid[1];
          splineGridSize[2] = m_WarpGrid[2];
          atlasToSubjectRegistrationHelper->SetSplineGridSize(splineGridSize);
          // Setting max displace
          atlasToSubjectRegistrationHelper->SetMaxBSplineDisplacement(6.0);
          }
        // else if( m_AtlastTransformType == "SyN" )
        //  {
        //   muLogMacro(
        //    << "Registering (SyN) " << preprocessMovingString << "atlas(" << vIndex << ") to template(" << vIndex
        //    << ") image." << std::endl);
        //  }
        atlasToSubjectRegistrationHelper->SetCurrentGenericTransform(m_TemplateGenericTransform);
        if( this->m_DebugLevel > 9 )
          {
          static unsigned int atlasToSubjectCounter = 0;
          std::stringstream   ss;
          ss << std::setw(3) << std::setfill('0') << atlasToSubjectCounter;
          atlasToSubjectRegistrationHelper->PrintCommandLine(true, std::string("AtlasToSubjectUpdate") + ss.str() );
          muLogMacro( << __FILE__ << " " << __LINE__ << " "  <<   std::endl );
          atlasToSubjectCounter++;
          }
        atlasToSubjectRegistrationHelper->Update();
        unsigned int actualIterations = atlasToSubjectRegistrationHelper->GetActualNumberOfIterations();
        muLogMacro( << "Registration tool " << actualIterations << " iterations." << std::endl );
#if 0   // ERROR:  This is not working correctly, because the proper number of
        // iterations is not reportd by the optimizer.
        while( actualIterations == 0 )
          {
          double newGradientTolerance = atlasToSubjectRegistrationHelper->GetProjectedGradientTolerance() * 0.1;
          atlasToSubjectRegistrationHelper->SetProjectedGradientTolerance(newGradientTolerance);
          muLogMacro( << "Reducing Gradient Tolerance to " << newGradientTolerance << std::endl );
          atlasToSubjectRegistrationHelper->Update();
          actualIterations = atlasToSubjectRegistrationHelper->GetActualNumberOfIterations();
          muLogMacro( << "Registration tool " << actualIterations << " iterations." << std::endl );
          }

#endif
        m_TemplateGenericTransform = atlasToSubjectRegistrationHelper->GetCurrentGenericTransform();
        }
      }
    }
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WritePartitionTable(const unsigned int CurrentEMIteration) const
{
  const unsigned int numChannels = this->m_CorrectedImages.size();
  const unsigned int numPriors = this->m_WarpedPriors.size();

  muLogMacro(<< "\n\nEM iteration " << CurrentEMIteration <<  std::endl);
  muLogMacro(<< "---------------------" << std::endl);
  PrettyPrintTable EMIterationTable;
  for( unsigned int ichan = 0; ichan < numChannels; ichan++ )
    {
    EMIterationTable.add(0, ichan + 2, this->m_InputVolumeTypes[ichan]);
    }
  for( unsigned int iclass = 0; iclass < numPriors; iclass++ )
    {
    EMIterationTable.add(iclass + 1, 0, std::string("Class ") + this->m_PriorNames[iclass] + " mean ");
    EMIterationTable.add(iclass + 1, 1, std::string(": ") );
    for( unsigned int ichan = 0; ichan < numChannels; ichan++ )
      {
      EMIterationTable.add(iclass + 1, ichan + 2, this->m_ListOfClassStatistics[iclass].m_Means[ichan], "%8.2f");
      }
    }
    {
    std::ostringstream oss;
    EMIterationTable.Print(oss);
    muLogMacro( << oss.str() );
    }
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::WriteDebugForegroundMask(const ByteImageType::Pointer & currForgroundBrainMask,
                           const unsigned int CurrentEMIteration) const
{
  std::stringstream CurrentEMIteration_stream("");

  CurrentEMIteration_stream << CurrentEMIteration;
  typedef itk::ImageFileWriter<ByteImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();

  const std::string fn = this->m_OutputDebugDir + "/MASK_LEVEL_" + CurrentEMIteration_stream.str() + ".nii.gz";
  writer->SetInput(currForgroundBrainMask);
  writer->SetFileName(fn.c_str() );
  writer->Update();
  muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl );
}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::Update()
{
  if( m_AtlasTransformType == "invalid_TransformationTypeNotSet" )
    {
    raise itk::Exception("The AtlasTransformType has NOT been set!")
    }
  if( m_UpdateRequired )
    {
    // TODO:  This should be filled out from the XML file eventually
    this->m_PriorsBackgroundValues.resize(this->m_OriginalSpacePriors.size() );
    std::fill(this->m_PriorsBackgroundValues.begin(), this->m_PriorsBackgroundValues.end(), 0);
      {
      // HACK:  In the XML file, each prior should also specify the default
      // background value
      //       to be used during the warping process.
      //       For AIR, it should be 1.0,  for all others is should be 0.0
      //
      m_PriorsBackgroundValues[this->GetAirIndex()] = 1;
      }

    this->EMLoop();
    m_UpdateRequired = false;
    }
}

/**
 * The EMLoop is the global algorithmic framework for
 * completing the iterative parts of the processing.
 */

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::EMLoop()
{
  if( this->m_TemplateGenericTransform.IsNull() )
    {
    itkExceptionMacro( << "ERROR:  Must suppply an intial transformation!" );
    }

  this->m_NonAirRegion = ComputeTissueRegion<TInputImage, ByteImageType>(this->m_InputImages[0], 3);
  if( this->m_DebugLevel > 9 )
    {
    this->WriteDebugHeadRegion(0);
    }

  this->m_WarpedPriors =
    WarpImageList<TProbabilityImage>(this->m_OriginalSpacePriors, this->m_InputImages[0],
                                     this->m_PriorsBackgroundValues,
                                     this->m_TemplateGenericTransform);
  if( this->m_DebugLevel > 9 )
    {
    this->WriteDebugWarpedAtlasPriors(0);
    std::vector<typename TInputImage::PixelType> backgroundValues(m_OriginalAtlasImages.size() );
    std::fill(backgroundValues.begin(), backgroundValues.end(), 0);
    this->m_WarpedAtlasImages =
      WarpImageList<TProbabilityImage>(this->m_OriginalAtlasImages, this->m_InputImages[0], backgroundValues,
                                       this->m_TemplateGenericTransform);
    this->WriteDebugWarpedAtlasImages(0);
    }
  typename ByteImageType::Pointer
  currForgroundBrainMask = ComputeForegroundProbMask<TProbabilityImage>(this->m_WarpedPriors,
                                                                        this->m_PriorIsForegroundPriorVector);
  if( this->m_DebugLevel > 9 )
    {
    WriteDebugForegroundMask(currForgroundBrainMask, 0);
    }

  std::vector<ByteImagePointer> SubjectCandidateRegions = this->UpdateIntensityBasedClippingOfPriors(
      0, this->m_InputImages, this->m_WarpedPriors, currForgroundBrainMask);
    {
    BlendPosteriorsAndPriors<TProbabilityImage>(0.0, this->m_WarpedPriors, this->m_WarpedPriors, this->m_WarpedPriors);
    std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
    NormalizeProbListInPlace<TProbabilityImage>(this->m_WarpedPriors);
    std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
    if( this->m_DebugLevel > 9 )
      {
      this->WriteDebugBlendClippedPriors(0);
      }
    }

  // NOTE:  Labels are only needed if debugging them.
  ComputeLabels<TProbabilityImage>(this->m_WarpedPriors, this->m_PriorIsForegroundPriorVector,
                                   this->m_PriorLabelCodeVector, this->m_NonAirRegion, this->m_DirtyLabels,
                                   this->m_CleanedLabels);
  this->WriteDebugLabels(0);
  this->m_ListOfClassStatistics.resize(0); // Reset this to empty for debugging
                                           // purposes to induce failures when
                                           // being re-used.
  this->m_CorrectedImages =
    CorrectBias(1, 0, SubjectCandidateRegions, this->m_InputImages, this->m_CleanedLabels, this->m_NonAirRegion,
                this->m_WarpedPriors, this->m_PriorUseForBiasVector, this->m_SampleSpacing, this->m_DebugLevel,
                this->m_OutputDebugDir);
  WriteDebugCorrectedImages(this->m_CorrectedImages, 0);
#if 0  // This is probably overkill since the update of NonAirRegion is likely
       // not making much in changes at this level of bias correction! Otsu
       // works well even in hightly biased images
  this->m_NonAirRegion = ComputeTissueRegion<TInputImage, ByteImageType>(this->m_CorrectedImages[0], 0);
  if( this->m_DebugLevel > 5 )
    {
    this->WriteDebugHeadRegion(0);
    }
#endif
  this->m_ListOfClassStatistics =
    ComputeDistributions<TInputImage, TProbabilityImage>(SubjectCandidateRegions, this->m_WarpedPriors,
                                                         this->m_CorrectedImages,
                                                         this->m_DebugLevel);
  this->WritePartitionTable(0);
    {
    // Now check that the intraSubjectOriginalImageList has positive definite
    // covariance matrix.
    // The algorithm is not stable if the covariance matrix is not positive
    // definite, and this
    // occurs when two or more of the images are linearly dependant (i.e. nearly
    // the same image).
    for( unsigned int q = 0; q < this->m_ListOfClassStatistics.size(); q++ )
      {
      try
        {
        ComputeCovarianceDeterminant( this->m_ListOfClassStatistics[q].m_Covariance );
        }
      catch( ... )
        {
        itkExceptionMacro( << "ERROR:\nERROR:\nERROR:\nERROR:"
                           << " Linearly dependant input images detected. "
                           << "Please remove the images in the above table that show very similar values images."
                           << "ERROR:\nERROR:\nERROR:\nERROR:" );
        }
      }
    }

  this->CheckInput();

  // FloatingPrecision logLikelihood = vnl_huge_val(1.0);
  FloatingPrecision logLikelihood = 1.0 / vnl_math::eps;
  FloatingPrecision deltaLogLikelihood = 1.0;

  unsigned int biasdegree = 0;

  // EM loop
  bool   converged = false;
  double priorWeighting = 1.00;       // NOTE:  This turns off blending of
                                      // posteriors and priors when set to 1.0,
                                      // thus short-circuting the system.
  unsigned int CurrentEMIteration = 1;
  while( !converged && ( CurrentEMIteration <= m_MaximumIterations ) )
    {
    // Recompute posteriors, not at full resolution
    this->m_Posteriors =
      ComputePosteriors<TInputImage, TProbabilityImage>(this->m_WarpedPriors, this->m_PriorWeights,
                                                        this->m_CorrectedImages,
                                                        this->m_ListOfClassStatistics);
    std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
    NormalizeProbListInPlace<TProbabilityImage>(this->m_Posteriors);
    std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
    this->WriteDebugPosteriors(CurrentEMIteration);
    ComputeLabels<TProbabilityImage>(this->m_Posteriors, this->m_PriorIsForegroundPriorVector,
                                     this->m_PriorLabelCodeVector, this->m_NonAirRegion, this->m_DirtyLabels,
                                     this->m_CleanedLabels);
    this->WriteDebugLabels(CurrentEMIteration);
    this->m_CorrectedImages =
      CorrectBias(this->m_MaxBiasDegree, CurrentEMIteration, SubjectCandidateRegions, this->m_InputImages,
                  this->m_CleanedLabels, this->m_NonAirRegion, this->m_Posteriors, this->m_PriorUseForBiasVector,
                  this->m_SampleSpacing, this->m_DebugLevel,
                  this->m_OutputDebugDir);
    WriteDebugCorrectedImages(this->m_CorrectedImages, CurrentEMIteration);
    this->m_ListOfClassStatistics.resize(0); // Reset this to empty for
                                             // debugging purposes to induce
                                             // failures when being re-used.
    this->m_ListOfClassStatistics =
      ComputeDistributions<TInputImage, TProbabilityImage>(SubjectCandidateRegions, this->m_Posteriors,
                                                           this->m_CorrectedImages,
                                                           this->m_DebugLevel);
    this->WritePartitionTable(CurrentEMIteration);

    // Now update transformation and estimates of probability regions based on
    // current knowledge.
      {
      this->UpdateTransformation(CurrentEMIteration); // This changes the class
                                                      // value of
                                                      //
                                                      // m_TemplateGenericTransform.
      this->m_WarpedPriors =
        WarpImageList<TProbabilityImage>(this->m_OriginalSpacePriors, this->m_InputImages[0],
                                         this->m_PriorsBackgroundValues,
                                         this->m_TemplateGenericTransform);
      if( this->m_DebugLevel > 9 )
        {
        this->WriteDebugWarpedAtlasPriors(CurrentEMIteration);
        std::vector<typename TInputImage::PixelType> backgroundValues(m_OriginalAtlasImages.size() );
        std::fill(backgroundValues.begin(), backgroundValues.end(), 0);
        this->m_WarpedAtlasImages =
          WarpImageList<TProbabilityImage>(this->m_OriginalAtlasImages, this->m_InputImages[0], backgroundValues,
                                           this->m_TemplateGenericTransform);
        this->WriteDebugWarpedAtlasImages(CurrentEMIteration);
        }
      SubjectCandidateRegions = this->ForceToOne(CurrentEMIteration, this->m_CorrectedImages, this->m_WarpedPriors,
                                                 this->m_CleanedLabels);
        {
        BlendPosteriorsAndPriors<TProbabilityImage>(1.0 - priorWeighting, this->m_Posteriors, this->m_WarpedPriors,
                                                    this->m_WarpedPriors);
        priorWeighting *= priorWeighting;
        std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
        NormalizeProbListInPlace<TProbabilityImage>(this->m_WarpedPriors);
        std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
        this->WriteDebugBlendClippedPriors(CurrentEMIteration);
        }
      }

    FloatingPrecision prevLogLikelihood = ( logLikelihood < vnl_math::eps ) ? vnl_math::eps : logLikelihood;
    // Compute log-likelihood and normalize posteriors
    logLikelihood = this->ComputeLogLikelihood();
    muLogMacro(<< "log(likelihood) = " << logLikelihood <<  std::endl);
    // TODO: move to before prevL update
    deltaLogLikelihood = vcl_fabs( (logLikelihood - prevLogLikelihood) / prevLogLikelihood);
    // (logLikelihood - prevLogLikelihood) / vcl_fabs(prevLogLikelihood);
    CHECK_NAN(deltaLogLikelihood, __FILE__, __LINE__,
              "\n logLikelihood: " << logLikelihood << "\n prevLogLikelihood: " << prevLogLikelihood );
    muLogMacro(
      << "delta vcl_log(likelihood) = " << deltaLogLikelihood << "  Convergence Tolerance: "
      << m_WarpLikelihoodTolerance <<  std::endl);

    // Convergence check
    converged = (CurrentEMIteration >= m_MaximumIterations)
      // Ignore jumps in the vcl_log likelihood
      //    ||
      //    (deltaLogLikelihood < 0)
      ||
      ( (deltaLogLikelihood < m_LikelihoodTolerance)
        &&
        (biasdegree == m_MaxBiasDegree) );

    CurrentEMIteration++;
    const float biasIncrementInterval = (m_MaximumIterations / (m_MaxBiasDegree + 1) );
    CHECK_NAN(biasIncrementInterval, __FILE__, __LINE__,
              "\n m_MaximumIterations: " << m_MaximumIterations << "\n  m_MaxBiasDegree: " << m_MaxBiasDegree );
    // Bias correction
    if( m_MaxBiasDegree > 0 )
      {
      if( (
            (deltaLogLikelihood < m_BiasLikelihoodTolerance)
            || ( CurrentEMIteration > (biasdegree + 1) * biasIncrementInterval) )
          &&
          (biasdegree < m_MaxBiasDegree)
          )
        {
        biasdegree++;
        }
      }
    } // end EM loop

  muLogMacro(<< "Done computing posteriors with " << CurrentEMIteration << " iterations" << std::endl);

  this->m_Posteriors =
    ComputePosteriors<TInputImage, TProbabilityImage>(this->m_WarpedPriors, this->m_PriorWeights,
                                                      this->m_CorrectedImages,
                                                      this->m_ListOfClassStatistics);
  std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
  NormalizeProbListInPlace<TProbabilityImage>(this->m_Posteriors);
  std::cout << "HACK HERE: " << __FILE__ << " " << __LINE__ << std::endl;
  this->WriteDebugPosteriors(CurrentEMIteration + 100);
  ComputeLabels<TProbabilityImage>(this->m_Posteriors, this->m_PriorIsForegroundPriorVector,
                                   this->m_PriorLabelCodeVector, this->m_NonAirRegion, this->m_DirtyLabels,
                                   this->m_CleanedLabels);
  FloatingPrecision inclusionThreshold = 0.75F;
  ComputeLabels<TProbabilityImage>(this->m_Posteriors, this->m_PriorIsForegroundPriorVector,
                                   this->m_PriorLabelCodeVector, this->m_NonAirRegion, this->m_DirtyThresholdedLabels,
                                   this->m_ThresholdedLabels, inclusionThreshold);
  this->WriteDebugLabels(CurrentEMIteration + 100);

  // Bias correction at full resolution, still using downsampled images
  // for computing the bias field coeficients
  if( m_MaxBiasDegree > 0 )
    {
    this->m_CorrectedImages =
      CorrectBias(biasdegree, CurrentEMIteration + 100, SubjectCandidateRegions, this->m_InputImages,
                  this->m_CleanedLabels, this->m_NonAirRegion, this->m_Posteriors, this->m_PriorUseForBiasVector,
                  this->m_SampleSpacing, this->m_DebugLevel,
                  this->m_OutputDebugDir);
    WriteDebugCorrectedImages(this->m_CorrectedImages, CurrentEMIteration + 100);
    this->m_ListOfClassStatistics.resize(0); // Reset this to empty for
                                             // debugging purposes to induce
                                             // failures when being re-used.
    this->m_ListOfClassStatistics =
      ComputeDistributions<TInputImage, TProbabilityImage>(SubjectCandidateRegions, this->m_Posteriors,
                                                           this->m_CorrectedImages,
                                                           this->m_DebugLevel);
#if 0  // This is probably overkill since the update of NonAirRegion is likely
       // not making much in changes at this level of bias correction! Otsu
       // works well even in hightly biased images
    this->m_NonAirRegion = ComputeTissueRegion<TInputImage, ByteImageType>(this->m_CorrectedImages[0], 0);
    if( this->m_DebugLevel > 9 )
      {
      this->WriteDebugHeadRegion(0);
      }
#endif
    this->m_RawCorrectedImages =
      CorrectBias(biasdegree, CurrentEMIteration + 100, SubjectCandidateRegions, this->m_RawInputImages,
                  this->m_CleanedLabels, this->m_NonAirRegion, this->m_Posteriors, this->m_PriorUseForBiasVector,
                  this->m_SampleSpacing, this->m_DebugLevel,
                  this->m_OutputDebugDir);
    }
  this->m_ListOfClassStatistics.resize(0); // Reset this to empty for debugging
                                           // purposes to induce failures when
                                           // being re-used.
  this->m_ListOfClassStatistics =
    ComputeDistributions<TInputImage, TProbabilityImage>(SubjectCandidateRegions, this->m_Posteriors,
                                                         this->m_CorrectedImages,
                                                         this->m_DebugLevel);
  this->WritePartitionTable(0 + 100);
}

#endif
