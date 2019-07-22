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
#ifndef __BRAINSABCUtilities__hxx__
#define __BRAINSABCUtilities__hxx__

#include "ExtractSingleLargestRegion.h"
#include "tbb/tbb.h"

template < typename TProbabilityImage >
void
ZeroNegativeValuesInPlace( std::vector< typename TProbabilityImage::Pointer > & priors )
{
  const unsigned int numPriors = priors.size();

  {
    tbb::parallel_for( tbb::blocked_range< LOOPITERTYPE >( 0, numPriors, 1 ),
                       [=]( tbb::blocked_range< LOOPITERTYPE > & r ) {
                         // First copy value, and set negative values to zero.
                         for ( LOOPITERTYPE iprior = r.begin(); iprior < r.end(); iprior++ )
                         {
                           for ( itk::ImageRegionIterator< TProbabilityImage > priorIter(
                                   priors[iprior], priors[iprior]->GetLargestPossibleRegion() );
                                 !priorIter.IsAtEnd();
                                 ++priorIter )
                           {
                             typename TProbabilityImage::PixelType inputValue( priorIter.Get() );
                             if ( inputValue < 0.0 )
                             {
                               priorIter.Set( 0.0 );
                             }
                           }
                         }
                       } );
  }
}

template < typename TProbabilityImage >
void
NormalizeProbListInPlace( std::vector< typename TProbabilityImage::Pointer > & ProbList )
{
  const unsigned int numProbs = ProbList.size();

  const typename TProbabilityImage::SizeType size = ProbList[0]->GetLargestPossibleRegion().GetSize();
  {
    tbb::parallel_for(
      tbb::blocked_range3d< LOOPITERTYPE >( 0, size[2], 1, 0, size[1], size[1] / 2, 0, size[0], 512 ),
      [=]( tbb::blocked_range3d< LOOPITERTYPE > & r ) {
        for ( LOOPITERTYPE kk = r.pages().begin(); kk < r.pages().end(); ++kk )
        {
          for ( LOOPITERTYPE jj = r.rows().begin(); jj < r.rows().end(); ++jj )
          {
            for ( LOOPITERTYPE ii = r.cols().begin(); ii < r.cols().end(); ++ii )
            {
              const typename TProbabilityImage::IndexType currIndex = { { ii, jj, kk } };
              FloatingPrecision                           sumPrior = 0.0;
              for ( unsigned int iprior = 0; iprior < numProbs; iprior++ )
              {
                const FloatingPrecision & ProbListValue = ProbList[iprior]->GetPixel( currIndex );
                CHECK_NAN( ProbListValue,
                           __FILE__,
                           __LINE__,
                           "\n  sumPrior: " << sumPrior << "\n  currIndex: " << currIndex
                                            << "\n ProbListValue: " << ProbListValue << "\n  iprior: " << iprior );
                sumPrior += ProbListValue;
              }
              if ( sumPrior < 1e-20 )
              {
                const FloatingPrecision averageValue = 1.0 / static_cast< FloatingPrecision >( numProbs );
                for ( unsigned int iprior = 0; iprior < numProbs; iprior++ )
                {
                  ProbList[iprior]->SetPixel( currIndex, averageValue );
                }
              }
              else
              {
                const FloatingPrecision invSumPrior = 1.0 / sumPrior;
                for ( unsigned int iprior = 0; iprior < numProbs; iprior++ )
                {
                  const FloatingPrecision normValue = ProbList[iprior]->GetPixel( currIndex ) * invSumPrior;

                  CHECK_NAN( normValue,
                             __FILE__,
                             __LINE__,
                             "\n  sumPrior: " << sumPrior << "\n  iprior: " << iprior << "\n  currIndex: " << currIndex
                                              << "\n probList: " << ProbList[iprior]->GetPixel( currIndex )
                                              << "\n  invSumPrior: " << invSumPrior );
                  ProbList[iprior]->SetPixel( currIndex, normValue );
                }
              }
            }
          }
        }
      } );
  }
}

template < typename TInputImage >
std::vector< typename TInputImage::Pointer >
DuplicateImageList( const std::vector< typename TInputImage::Pointer > & inputList )
{
  std::vector< typename TInputImage::Pointer > outputList( inputList.size() );
  {
    tbb::parallel_for( tbb::blocked_range< LOOPITERTYPE >( 0, inputList.size(), 1 ),
                       [=, &outputList]( const tbb::blocked_range< LOOPITERTYPE > & r ) {
                         for ( auto i = r.begin(); i < r.end(); i++ )
                         {
                           typename itk::ImageDuplicator< TInputImage >::Pointer myDuplicator =
                             itk::ImageDuplicator< TInputImage >::New();
                           myDuplicator->SetInputImage( inputList[i] );
                           myDuplicator->Update();
                           outputList[i] = myDuplicator->GetOutput();
                         }
                       } );
  }

  return outputList;
}


template < typename TProbabilityImage >
typename ByteImageType::Pointer
ComputeForegroundProbMask( const std::vector< typename TProbabilityImage::Pointer > & probList,
                           const std::vector< bool > &                                IsForegroundPriorVector )
{
  muLogMacro( << "ComputeForegroundProbMask" << std::endl );
  const unsigned int              numPriors = probList.size();
  typename ByteImageType::Pointer currForegroundMask = ByteImageType::New();
  currForegroundMask->CopyInformation( probList[0] );
  currForegroundMask->SetRegions( probList[0]->GetLargestPossibleRegion() );
  currForegroundMask->Allocate();

  const typename TProbabilityImage::SizeType  size = probList[0]->GetLargestPossibleRegion().GetSize();
  constexpr typename ByteImageType::PixelType insideMaskValue = 1;
  {
    tbb::parallel_for( tbb::blocked_range3d< LOOPITERTYPE >( 0, size[2], 1, 0, size[1], size[1] / 2, 0, size[0], 512 ),
                       [=]( tbb::blocked_range3d< LOOPITERTYPE > & r ) {
                         for ( LOOPITERTYPE kk = r.pages().begin(); kk < r.pages().end(); ++kk )
                         {
                           for ( LOOPITERTYPE jj = r.rows().begin(); jj < r.rows().end(); ++jj )
                           {
                             for ( LOOPITERTYPE ii = r.cols().begin(); ii < r.cols().end(); ++ii )
                             {
                               const typename TProbabilityImage::IndexType currIndex = { { ii, jj, kk } };
                               FloatingPrecision                           tmp = 0.0;
                               for ( unsigned int iprior = 0; iprior < numPriors; iprior++ )
                               {
                                 const bool fgflag = IsForegroundPriorVector[iprior];
                                 if ( fgflag == true )
                                 {
                                   tmp += probList[iprior]->GetPixel( currIndex );
                                 }
                               }
                               if ( tmp > 0.5 ) // Only include if the sum of the non-background
                                                // priors are greater than 50 %
                               {
                                 currForegroundMask->SetPixel(
                                   currIndex, static_cast< typename ByteImageType::PixelType >( insideMaskValue ) );
                               }
                               else
                               {
                                 currForegroundMask->SetPixel( currIndex, 0 );
                               }
                             }
                           }
                         }
                       } );
  }
  return currForegroundMask;
}

#endif // __BRAINSABCUtilities__hxx__
