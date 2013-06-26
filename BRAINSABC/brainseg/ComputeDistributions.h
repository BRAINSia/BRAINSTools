#ifndef __ComputeDistributions__h_
#define __ComputeDistributions__h_
#include "BRAINSABCUtilities.h"
#include <vector>
#include <list>
#include <map>
#define EXPP(x) vcl_exp( ( x ) )
#define LOGP(x) vcl_log( ( x ) )

typedef  itk::Image<unsigned char, 3> ByteImageType;

template <class TInputImage, class TProbabilityImage, class MatrixType>
void
CombinedComputeDistributions( const std::vector<typename ByteImageType::Pointer> & SubjectCandidateRegions,
                              const std::map<std::string,std::vector<typename TInputImage::Pointer> >
                              &InputImageMap,
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
  typedef std::vector<typename TInputImage::Pointer> InputImageVector;
  typedef std::map<std::string,InputImageVector> MapOfInputImageVectors;

  const LOOPITERTYPE numClasses =     PosteriorsList.size();
  const LOOPITERTYPE numModalities = InputImageMap.size();

  ListOfClassStatistics.clear();
  ListOfClassStatistics.resize(numClasses);

  // not sure this is needed -- this sets the size of the
  // covariance matrix, but that is overwritten by assignment below
  // once the covariance has been computed.
  for( LOOPITERTYPE iclass = 0; iclass < numClasses; iclass++ )
    {
    ListOfClassStatistics[iclass].resize(numModalities);
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
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
  for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
    {
    const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
    const typename ByteImageType::ConstPointer     currentCandidateRegion =
      SubjectCandidateRegions[iclass].GetPointer();
    ListOfClassStatistics[iclass].m_Means.clear();

    for(typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
        mapIt != InputImageMap.end(); ++mapIt)
      {
      unsigned meanIndex(0);

      ListOfClassStatistics[iclass].m_Means[mapIt->first] = 0.0;

      for(typename InputImageVector::const_iterator imIt = mapIt->second.begin();
          imIt != mapIt->second.end(); ++imIt, ++meanIndex)
        {
        typename TInputImage::Pointer im1 = *imIt;
        double muSum = 0.0;
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
                const double currentInputValue = im1->GetPixel(currIndex);
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
        double __mean__(( muSum ) / ListOfClassStatistics[iclass].m_Weighting );
        ListOfClassStatistics[iclass].m_Means[mapIt->first] += __mean__;
        }
      // averaging the means of all images of this image modality
      ListOfClassStatistics[iclass].m_Means[mapIt->first] /= mapIt->second.size();
      }
    }

  std::vector<MatrixType> oldCovariances(ListOfClassStatistics.size() );
  if( (LOOPITERTYPE)oldCovariances.size() != numClasses )
    {
    oldCovariances.clear();
    oldCovariances.resize(numClasses);
    for( LOOPITERTYPE iclass = 0; iclass < numClasses; iclass++ )
      {
      MatrixType C(numModalities, numModalities);
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
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
  for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
    {
    const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
    const typename ByteImageType::ConstPointer     currentCandidateRegion =
      SubjectCandidateRegions[iclass].GetPointer();
    //
    // this will end up as a vnl_matrix for assignment to
    // the Class Statistics object after this is computed.
    std::map<std::string, std::map<std::string, double> > TypeCovariance;
    // initialize -- no easy way since it is a map of maps
    for(typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
        mapIt != InputImageMap.end(); ++mapIt)
      {
      for(typename MapOfInputImageVectors::const_iterator mapIt2 = InputImageMap.begin();
          mapIt2 != InputImageMap.end(); ++mapIt2)
        {
        TypeCovariance[mapIt->first][mapIt2->first] = 0.0;
        }
      }
    //
    // compute per-Image Type covariance
    for(typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
        mapIt != InputImageMap.end(); ++mapIt)
      {
      const double mu1 =
        ListOfClassStatistics[iclass].m_Means[mapIt->first];

      for(unsigned i = 0; i < mapIt->second.size(); ++i)
        {
        typename TInputImage::Pointer im1 = mapIt->second[i];

        bool first_through_inner_loop(true);

        for(typename MapOfInputImageVectors::const_iterator mapIt2 = mapIt;
            mapIt2 != InputImageMap.end(); ++mapIt2)
          {
          size_t j = 0;
          if(first_through_inner_loop)
            {
            j = i;
            first_through_inner_loop = false;
            }
          const double mu2 =
            ListOfClassStatistics[iclass].m_Means[mapIt2->first];
          for (; j < mapIt2->second.size(); ++j)
            {
            typename TInputImage::Pointer im2 = mapIt2->second[i];
            double var = 0.0;
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
                      const double diff1 = LOGP( static_cast<double>( im1->GetPixel(currIndex) ) ) - mu1;
                      const double diff2 = LOGP( static_cast<double>( im2->GetPixel(currIndex) ) ) - mu2;
                      var += currentProbValue * ( diff1 * diff2 );
                      }
                    else
                      {
                      const double diff1 = ( static_cast<double>( im1->GetPixel(currIndex) ) ) - mu1;
                      const double diff2 = ( static_cast<double>( im2->GetPixel(currIndex) ) ) - mu2;
                      var += currentProbValue * ( diff1 * diff2 );
                      }
                    }
                  }
                }
              }
            var /= ListOfClassStatistics[iclass].m_Weighting;

            // Adjust diagonal, to make sure covariance is pos-def
            if(mapIt == mapIt2 && i == j )
              {
              var += 1e-20;
              }
            // accumulate
            TypeCovariance[mapIt->first][mapIt2->first] += var;
            TypeCovariance[mapIt2->first][mapIt->first] += var;
            }
          }
        }
      }
    // above loop accumulates covariances per type
    // now divide out # of averaged variances
    // and copy to vnl matrix
    MatrixType covtmp(numModalities,numModalities, 0.0);
    unsigned int i = 0;
    for(typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
        mapIt != InputImageMap.end(); ++mapIt,++i)
      {
      unsigned int j = 0;
      for(typename MapOfInputImageVectors::const_iterator mapIt2 = InputImageMap.begin();
          mapIt2 != InputImageMap.end(); ++mapIt2, ++j)
        {
        covtmp(i,j) = TypeCovariance[mapIt->first][mapIt2->first] /
          static_cast<double>(mapIt->second.size() * mapIt2->second.size());
        }
      }
      ListOfClassStatistics[iclass].m_Covariance = covtmp;
    } // end covariance loop
  if( DebugLevel > 5 )
    {
    for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
      {
      muLogMacro(<< "DEBUG USING NEW COVARIANCES: " << iclass << std::endl
                 << ListOfClassStatistics[iclass].m_Covariance << std::endl );
      }
    }
  if( DebugLevel > 9 )
    {
    std::cout << "=================================================" << std::endl;
    unsigned ichan = 0;
    for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
      {
      for(typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
          mapIt != InputImageMap.end(); ++mapIt)
        {
        muLogMacro( << "DEBUG MEAN " << ichan << " : " << iclass << " : \n"
                    << ListOfClassStatistics[iclass].m_Means[mapIt->first]
                    << " \n" << std::endl );
        }
      muLogMacro( << "DEBUG Covariances: " << iclass << "\n" << ListOfClassStatistics[iclass].m_Covariance
                  << std::endl );
      }
    }
}

#endif // __ComputeDistributions__h__
