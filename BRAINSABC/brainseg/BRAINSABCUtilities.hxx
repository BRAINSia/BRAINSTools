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
            sumPrior += ProbList[iprior]->GetPixel(currIndex);
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
              CHECK_NAN(normValue, __FILE__, __LINE__);
              ProbList[iprior]->SetPixel(currIndex, normValue);
              }
            }
          }
        }
      }
    }
}

// Labeling using maximum a posteriori, also do brain stripping using
// mathematical morphology and connected component
template <class TProbabilityImage>
void ComputeLabels(
  std::vector<typename TProbabilityImage::Pointer> & Posteriors,
  std::vector<bool> & PriorIsForegroundPriorVector,
  vnl_vector<unsigned int> & PriorLabelCodeVector,
  typename ByteImageType::Pointer & NonAirRegion,
  typename ByteImageType::Pointer & DirtyLabels,
  typename ByteImageType::Pointer & CleanedLabels)
{
  muLogMacro(<< "ComputeLabels" << std::endl );
  itk::TimeProbe ComputeLabelsTimer;
  ComputeLabelsTimer.Start();

  const unsigned int numClasses = Posteriors.size();
  const typename TProbabilityImage::RegionType region = Posteriors[0]->GetLargestPossibleRegion();

  DirtyLabels = ByteImageType::New();
  DirtyLabels->CopyInformation(Posteriors[0]);
  DirtyLabels->SetRegions(region);
  DirtyLabels->Allocate();
  DirtyLabels->FillBuffer(0);

  typename ByteImageType::Pointer foregroundMask = ByteImageType::New();
  foregroundMask->CopyInformation(Posteriors[0]);
  foregroundMask->SetRegions(region);
  foregroundMask->Allocate();
  foregroundMask->FillBuffer(0);

  const typename ByteImageType::SizeType size = DirtyLabels->GetLargestPossibleRegion().GetSize();
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
          if( NonAirRegion->GetPixel(currIndex) == 0 ) // If outside the tissue
                                                       // region, then set to
                                                       // zero vIndex!
            {
            // TODO:  May want to specify this explicitly in the XML file for
            // the proper background value
            DirtyLabels->SetPixel(currIndex, 0); // This is implied by the
                                                 // FillBuffer(0) above;
            continue;
            }

          FloatingPrecision maxPosteriorClassValue = Posteriors[0]->GetPixel(currIndex);
          unsigned int      indexMaxPosteriorClassValue = 0;
          for( unsigned int iclass = 1; iclass < numClasses; iclass++ )
            {
            const FloatingPrecision currentPosteriorClassValue = Posteriors[iclass]->GetPixel(currIndex);
            if( currentPosteriorClassValue > maxPosteriorClassValue )
              {
              maxPosteriorClassValue = currentPosteriorClassValue;
              indexMaxPosteriorClassValue = iclass;
              }
            }

            {
            bool         fgflag = PriorIsForegroundPriorVector[indexMaxPosteriorClassValue];
            unsigned int label = PriorLabelCodeVector[indexMaxPosteriorClassValue];
            // Only use non-zero probabilities and foreground classes
            if( !fgflag || ( maxPosteriorClassValue < 0.001 ) )
              {
              fgflag = false; // If priors are zero or negative, then set the
                              // fgflag back to false
              }
            DirtyLabels->SetPixel(currIndex, label);
            foregroundMask->SetPixel(currIndex, fgflag);
            }
          }
        }
      }
    }
  //
  // CleanedLabels=ExtractSingleLargestRegionFromMask(foregroundMask,2,2,1,DirtyLabels);
  CleanedLabels = ExtractSingleLargestRegionFromMask(foregroundMask, 0, 0, 0, DirtyLabels);
  ComputeLabelsTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime =
    ComputeLabelsTimer.GetTotal();
  muLogMacro(<< "Computing Labels took " << elapsedTime << " " << ComputeLabelsTimer.GetUnit() << std::endl);
}

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
      outputList[i] = myDuplicator->GetOutput();
      }
    }

  return outputList;
}

#if 0
template <class ImageType>
void WinnerTakesAll(std::vector<typename ImageType::Pointer> & listOfImages)
{
  itk::ImageRegionIteratorWithIndex<ImageType> it(listOfImages[0], listOfImages[0]->GetLargestPossibleRegion() );

  while( !it.IsAtEnd() )
    {
    const typename ImageType::IndexType currIndex = it.GetIndex();
    unsigned int currMaxIndex = 0;
    typename ImageType::PixelType currMax = listOfImages[0]->GetPixel(currIndex);
    typename ImageType::PixelType sum = currMax;
    listOfImages[0]->SetPixel(currIndex, 0);
    for( unsigned int im = 1; im < listOfImages.size(); im++ )
      {
      const typename ImageType::PixelType currValue = listOfImages[im]->GetPixel(currIndex);
      sum += currValue;
      if( currValue > currMax )
        {
        currMax = currValue;
        currMaxIndex = im;
        }
      listOfImages[im]->SetPixel(currIndex, 0);
      }
    listOfImages[currMaxIndex]->SetPixel(currIndex, sum);
    ++it;
    }

  return;
}

#endif

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
#if 0
    {
    // Pre-Dilate mask
    typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
    typedef
      itk::BinaryDilateImageFilter<ByteImageType, ByteImageType,
                                   StructElementType> DilateType;

    StructElementType structel;
    structel.SetRadius(1);
    structel.CreateStructuringElement();

    typename DilateType::Pointer dil = DilateType::New();
    dil->SetDilateValue(50);
    dil->SetKernel(structel);
    dil->SetInput(currForegroundMask);

    dil->Update();

    ByteImagePointer dilmask = dil->GetOutput();
    // a simple assignment is probably sufficient, test both ways?
    currForegroundMask = CopyImage<ByteImageType>(dilmask);
    }
#endif
  return currForegroundMask;
}

#endif // __BRAINSABCUtilities__hxx__
