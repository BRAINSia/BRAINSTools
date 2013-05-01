#ifndef __BRAINSABCUtilities__hxx__
#define __BRAINSABCUtilities__hxx__

#include "ExtractSingleLargestRegion.h"

template <class TProbabilityImage>
void ZeroNegativeValuesInPlace(std::vector<typename TProbabilityImage::Pointer> & priors)
{
  const unsigned int numPriors = priors.size();

    {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for
#endif
    // First copy value, and set negative values to zero.
    for( LOOPITERTYPE iprior = 0; iprior < (LOOPITERTYPE)numPriors; iprior++ )
      {
      for( itk::ImageRegionIterator<TProbabilityImage> priorIter(priors[iprior],
                                                                 priors[iprior]->GetLargestPossibleRegion() );
           !priorIter.IsAtEnd(); ++priorIter )
        {
        typename TProbabilityImage::PixelType inputValue(priorIter.Get() );
        if( inputValue < 0.0 )
          {
          priorIter.Set(0.0);
          }
        }
      }
    }
}

template <class TProbabilityImage>
void NormalizeProbListInPlace(std::vector<typename TProbabilityImage::Pointer> & ProbList)
{
  const unsigned int numProbs = ProbList.size();

  const typename TProbabilityImage::SizeType size = ProbList[0]->GetLargestPossibleRegion().GetSize();
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
          FloatingPrecision sumPrior = 0.0;
          for( unsigned int iprior = 0; iprior < numProbs; iprior++ )
            {
            const FloatingPrecision & ProbListValue = ProbList[iprior]->GetPixel(currIndex);
            CHECK_NAN(ProbListValue, __FILE__, __LINE__, "\n  sumPrior: " << sumPrior
                                                                          << "\n  currIndex: " << currIndex << "\n ProbListValue: "
                                                                          << ProbListValue << "\n  iprior: "
                                                                          << iprior );
            sumPrior += ProbListValue;
            }
          if( sumPrior < 1e-20 )
            {
            const FloatingPrecision averageValue = 1.0 / static_cast<FloatingPrecision>(numProbs);
            for( unsigned int iprior = 0; iprior < numProbs; iprior++ )
              {
              ProbList[iprior]->SetPixel(currIndex, averageValue);
              }
            }
          else
            {
            const FloatingPrecision invSumPrior = 1.0 / sumPrior;
            for( unsigned int iprior = 0; iprior < numProbs; iprior++ )
              {
              const FloatingPrecision normValue = ProbList[iprior]->GetPixel(currIndex) * invSumPrior;

              CHECK_NAN(normValue, __FILE__, __LINE__, "\n  sumPrior: "
                        << sumPrior << "\n  iprior: " << iprior << "\n  currIndex: "
                        << currIndex << "\n probList: "
                        << ProbList[iprior]->GetPixel(currIndex)
                        << "\n  invSumPrior: " << invSumPrior );
              ProbList[iprior]->SetPixel(currIndex, normValue);
              }
            }
          }
        }
      }
    }
}

#include "BRAINSComputeLabels.h"

template <class TInputImage>
std::vector<typename TInputImage::Pointer>
DuplicateImageList(const std::vector<typename TInputImage::Pointer> & inputList)
{
  std::vector<typename TInputImage::Pointer> outputList(inputList.size() );
    {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for
#endif
    for( LOOPITERTYPE i = 0; i < (LOOPITERTYPE)inputList.size(); i++ )
      {
      typename itk::ImageDuplicator<TInputImage>::Pointer myDuplicator
        = itk::ImageDuplicator<TInputImage>::New();
      myDuplicator->SetInputImage(inputList[i]);
      myDuplicator->Update();
      outputList[i] = myDuplicator->GetModifiableOutput();
      }
    }

  return outputList;
}


template <class TProbabilityImage>
typename ByteImageType::Pointer ComputeForegroundProbMask(
  const std::vector<typename TProbabilityImage::Pointer> & probList, const std::vector<bool> & IsForegroundPriorVector )
{
  muLogMacro(<< "ComputeForegroundProbMask" << std::endl );
  const unsigned int numPriors = probList.size();
  typename ByteImageType::Pointer currForegroundMask = ByteImageType::New();
  currForegroundMask->CopyInformation(probList[0]);
  currForegroundMask->SetRegions(probList[0]->GetLargestPossibleRegion() );
  currForegroundMask->Allocate();

  const typename TProbabilityImage::SizeType size = probList[0]->GetLargestPossibleRegion().GetSize();
  const typename ByteImageType::PixelType insideMaskValue = 1;
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
          FloatingPrecision tmp = 0.0;
          for( unsigned int iprior = 0; iprior < numPriors; iprior++ )
            {
            const bool fgflag = IsForegroundPriorVector[iprior];
            if( fgflag == true )
              {
              tmp += probList[iprior]->GetPixel(currIndex);
              }
            }
          if( tmp > 0.5 ) // Only include if the sum of the non-background
                          // priors are greater than 50 %
            {
            currForegroundMask->SetPixel(currIndex, static_cast<typename ByteImageType::PixelType>(insideMaskValue) );
            }
          else
            {
            currForegroundMask->SetPixel(currIndex, 0);
            }
          }
        }
      }
    }
  return currForegroundMask;
}

#endif // __BRAINSABCUtilities__hxx__
