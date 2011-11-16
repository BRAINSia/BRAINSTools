#include "BRAINSABCUtilities.h"
/*****************************
 * Now call the instantiations
 */
#include "BRAINSABCUtilities.hxx"
#include "LLSBiasCorrector.h"

template std::vector<FloatImageType::Pointer> DuplicateImageList<FloatImageType>(
  const std::vector<FloatImageType::Pointer> & );

template std::vector<ShortImageType::Pointer> DuplicateImageList<ShortImageType>(
  const std::vector<ShortImageType::Pointer> & );

template void ComputeLabels<FloatImageType>( std::vector<FloatImageType::Pointer> &, std::vector<bool> &,
                                             vnl_vector<unsigned int> &, ByteImageType::Pointer &,
                                             ByteImageType::Pointer &, ByteImageType::Pointer & );

template void NormalizeProbListInPlace<FloatImageType>(std::vector<FloatImageType::Pointer> & );

template void ZeroNegativeValuesInPlace<FloatImageType>(  std::vector<FloatImageType::Pointer> & );

std::vector<CorrectIntensityImageType::Pointer> CorrectBias(
  const unsigned int degree,
  const unsigned int CurrentEMIteration,
  const std::vector<ByteImageType::Pointer> &
  CandidateRegions,
  const std::vector<CorrectIntensityImageType::Pointer> &
  inputImages,
  const ByteImageType::Pointer currentBrainMask,
  const ByteImageType::Pointer currentTissueMask,
  const std::vector<FloatImageType::Pointer> & probImages,
  const std::vector<bool> & probUseForBias,
  const FloatingPrecision sampleSpacing,
  const int DebugLevel,
  const std::string& OutputDebugDir
  )
{
  std::vector<CorrectIntensityImageType::Pointer> correctedImages(inputImages.size() );

  if( degree == 0 )
    {
    muLogMacro(<< "Skipping Bias correction, polynomial degree = " << degree <<  std::endl);
    return inputImages;
    }
  muLogMacro(<< "Bias correction, polynomial degree = " << degree <<  std::endl);

  // Perform bias correction
  const unsigned int                   numClasses = probImages.size();
  std::vector<FloatImageType::Pointer> biasPosteriors;
  std::vector<ByteImageType::Pointer>  biasCandidateRegions;
  biasPosteriors.clear();
  biasCandidateRegions.clear();
    {
    for( unsigned int iclass = 0; iclass < numClasses; iclass++ )
      {
      const unsigned iprior = iclass;
      if( probUseForBias[iprior] == 1 )
        {
        // Focus only on FG classes, more accurate if bg classification is bad
        // but sacrifices accuracy in border regions (tend to overcorrect)
        biasPosteriors.push_back(probImages[iclass]);
        biasCandidateRegions.push_back(CandidateRegions[iclass]);
        }
      }
    }

  itk::TimeProbe BiasCorrectorTimer;
  BiasCorrectorTimer.Start();
  typedef LLSBiasCorrector<CorrectIntensityImageType, FloatImageType> BiasCorrectorType;
  typedef BiasCorrectorType::Pointer                                  BiasCorrectorPointer;

  BiasCorrectorPointer biascorr = BiasCorrectorType::New();
  biascorr->SetMaxDegree(degree);
  // biascorr->SetMaximumBiasMagnitude(5.0);
  // biascorr->SetSampleSpacing(2.0*SampleSpacing);
  biascorr->SetSampleSpacing(1);
  biascorr->SetWorkingSpacing(sampleSpacing);
  biascorr->SetForegroundBrainMask(currentBrainMask);
  biascorr->SetAllTissueMask(currentTissueMask);
  biascorr->SetProbabilities(biasPosteriors, biasCandidateRegions);
  biascorr->SetDebugLevel(DebugLevel);
  biascorr->SetOutputDebugDir(OutputDebugDir);

  if( DebugLevel > 0 )
    {
    biascorr->DebugOn();
    }

  biascorr->SetInputImages(inputImages);
  correctedImages = biascorr->CorrectImages(CurrentEMIteration);

  BiasCorrectorTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime = BiasCorrectorTimer.GetTotal();
  muLogMacro(<< "Computing BiasCorrection took " << elapsedTime << " " << BiasCorrectorTimer.GetUnit() << std::endl);

  return correctedImages;
}
