#ifndef __ComputeDistributions__h_
#define __ComputeDistributions__h_

#include <vector>

#define EXPP(x) vcl_exp( ( x ) )
#define LOGP(x) vcl_log( ( x ) )

typedef  itk::Image<unsigned char, 3> ByteImageType;

template <class TInputImage, class TProbabilityImage, class MatrixType>
void
CombinedComputeDistributions( const std::vector<typename ByteImageType::Pointer> & SubjectCandidateRegions,
                              const std::vector<typename TInputImage::Pointer> & InputImagesList,
                              const std::vector<typename TProbabilityImage::Pointer> & PosteriorsList,
                              std::vector<RegionStats> & ListOfClassStatistics, //
                                                                                //
                                                                                // This
                                                                                //
                                                                                // is
                                                                                //
                                                                                // an
                                                                                //
                                                                                // output!
                              const unsigned int DebugLevel,
                              const bool logConvertValues
                              )
{
  const LOOPITERTYPE numChannels = InputImagesList.size();
  const LOOPITERTYPE numClasses = PosteriorsList.size();

  ListOfClassStatistics.clear();
  ListOfClassStatistics.resize(numClasses);
  for( LOOPITERTYPE iclass = 0; iclass < numClasses; iclass++ )
    {
    ListOfClassStatistics[iclass].resize(numChannels);
    }

  typename TInputImage::SizeType size
    = PosteriorsList[0]->GetLargestPossibleRegion().GetSize();

  // Compute sum of posteriors for each class
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
  for( LOOPITERTYPE iclass = 0; iclass < numClasses; iclass++ )
    {
    const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
    const typename ByteImageType::ConstPointer     currentCandidateRegion =
      SubjectCandidateRegions[iclass].GetPointer();
    double tmp = 1e-20; // NOTE:  vnl_math:eps is too small vnl_math::eps;
      {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for reduction(+:tmp) default(shared)
#endif
      for( long kk = 0; kk < (long)size[2]; kk++ )
        {
        for( long jj = 0; jj < (long)size[1]; jj++ )
          {
          for( long ii = 0; ii < (long)size[0]; ii++ )
            {
            const typename TProbabilityImage::IndexType currIndex = {{ii, jj, kk}};
            if( currentCandidateRegion->GetPixel(currIndex) )
              {
              const double currentProbValue = currentProbImage->GetPixel(currIndex);
              tmp = tmp + currentProbValue;
              }
            }
          }
        }
      }
    ListOfClassStatistics[iclass].m_Weighting = tmp;
    }
  // Compute the means weighted by the probability of each value.
    {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
    for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
      {
      const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
      const typename ByteImageType::ConstPointer     currentCandidateRegion =
        SubjectCandidateRegions[iclass].GetPointer();
      for( LOOPITERTYPE ichan = 0; ichan < numChannels; ichan++ )
        {
        double muSum = 0.0;
          {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared) reduction(+:muSum)
#endif
          for( long kk = 0; kk < (long)size[2]; kk++ )
            {
            for( long jj = 0; jj < (long)size[1]; jj++ )
              {
              for( long ii = 0; ii < (long)size[0]; ii++ )
                {
                const typename TProbabilityImage::IndexType currIndex = {{ii, jj, kk}};
                if( currentCandidateRegion->GetPixel(currIndex) )
                  {
                  const double currentProbValue = currentProbImage->GetPixel(currIndex);
                  const double currentInputValue = InputImagesList[ichan]->GetPixel(currIndex);
                  if( logConvertValues )
                    {
                    muSum += currentProbValue * LOGP(currentInputValue);
                    }
                  else
                    {
                    muSum += currentProbValue * ( currentInputValue );
                    }
                  }
                }
              }
            }
          }
        ListOfClassStatistics[iclass].m_Means[ichan] =  ( muSum ) / ( ListOfClassStatistics[iclass].m_Weighting );
        }
      }
    }

  // Compute the covariances
  std::vector<MatrixType> oldCovariances(ListOfClassStatistics.size() );
  if( (LOOPITERTYPE)oldCovariances.size() != numClasses )
    {
    oldCovariances.clear();
    oldCovariances.resize(numClasses);
    for( LOOPITERTYPE iclass = 0; iclass < numClasses; iclass++ )
      {
      MatrixType C(numChannels, numChannels);
      C.set_identity();
      C *= 1e-10;
      oldCovariances[iclass] = C;
      }
    }
  else // Copy from previous version.
    {
    for( LOOPITERTYPE iclass = 0; iclass < numClasses; iclass++ )
      {
      oldCovariances[iclass] = ListOfClassStatistics[iclass].m_Covariance;
      }
    }
    {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
    for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
      {
      const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
      const typename ByteImageType::ConstPointer     currentCandidateRegion =
        SubjectCandidateRegions[iclass].GetPointer();
      MatrixType covtmp(numChannels, numChannels, 0.0);
      for( LOOPITERTYPE r = 0; r < numChannels; r++ )
        {
        const double mu1 = ListOfClassStatistics[iclass].m_Means[r];
        typename TInputImage::Pointer img1 = InputImagesList[r];
        for( LOOPITERTYPE c = r; c < numChannels; c++ )
          {
          const double mu2 = ListOfClassStatistics[iclass].m_Means[c];
          typename TInputImage::Pointer img2 = InputImagesList[c];
          double var = 0.0;
            {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared) reduction(+:var)
#endif
            for( long kk = 0; kk < (long)size[2]; kk++ )
              {
              for( long jj = 0; jj < (long)size[1]; jj++ )
                {
                for( long ii = 0; ii < (long)size[0]; ii++ )
                  {
                  const typename TInputImage::IndexType currIndex = {{ii, jj, kk}};
                  if( currentCandidateRegion->GetPixel(currIndex) )
                    {
                    const double currentProbValue = currentProbImage->GetPixel(currIndex);
                    if( logConvertValues )
                      {
                      const double diff1 = LOGP( static_cast<double>( img1->GetPixel(currIndex) ) ) - mu1;
                      const double diff2 = LOGP( static_cast<double>( img2->GetPixel(currIndex) ) ) - mu2;
                      var += currentProbValue * ( diff1 * diff2 );
                      }
                    else
                      {
                      const double diff1 = ( static_cast<double>( img1->GetPixel(currIndex) ) ) - mu1;
                      const double diff2 = ( static_cast<double>( img2->GetPixel(currIndex) ) ) - mu2;
                      var += currentProbValue * ( diff1 * diff2 );
                      }
                    }
                  }
                }
              }
            }
          var /= ListOfClassStatistics[iclass].m_Weighting;

          // Adjust diagonal, to make sure covariance is pos-def
          if( r == c )
            {
            var += 1e-20;
            }

          // Assign value to the covariance matrix (symmetric)
          covtmp(r, c) = var;
          covtmp(c, r) = var;
          }
        }
        {
        ListOfClassStatistics[iclass].m_Covariance = covtmp;
        }
      } // end covariance loop
    }
  if( DebugLevel > 5 )
    {
    for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
      {
      muLogMacro(
        << "DEBUG USING NEW COVARIANCES: " << iclass << "\n" << ListOfClassStatistics[iclass].m_Covariance
        << std::endl );
      }
    }
  if( DebugLevel > 9 )
    {
    std::cout << "=================================================" << std::endl;
    for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
      {
      for( LOOPITERTYPE ichan = 0; ichan < numChannels; ichan++ )
        {
        muLogMacro( << "DEBUG MEAN " << ichan << " : " << iclass << " : \n"
                    << ListOfClassStatistics[iclass].m_Means[ichan] << " \n" << std::endl );
        }
      muLogMacro( << "DEBUG Covariances: " << iclass << "\n" << ListOfClassStatistics[iclass].m_Covariance
                  << std::endl );
      }
    }
}

#endif // __ComputeDistributions__h__
