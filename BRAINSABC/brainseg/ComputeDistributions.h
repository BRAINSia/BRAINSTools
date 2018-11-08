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
#ifndef __ComputeDistributions__h_
#define __ComputeDistributions__h_
#include "BRAINSABCUtilities.h"
#include <vector>
#include <list>
#include <map>
#define EXPP(x) std::exp( ( x ) )
#define LOGP(x) std::log( ( x ) )

typedef  itk::Image<unsigned char, 3> ByteImageType;
typedef itk::CompensatedSummation<double> CompensatedSummationType;

template <typename TInputImage, typename TProbabilityImage, typename MatrixType>
void
CombinedComputeDistributions( const std::vector<typename ByteImageType::Pointer> & SubjectCandidateRegions,
                              const orderedmap<std::string,std::vector<typename TInputImage::Pointer> >
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
  typedef orderedmap<std::string,InputImageVector>   MapOfInputImageVectors;

  typedef itk::NearestNeighborInterpolateImageFunction< TInputImage, double > InputImageNNInterpolationType;

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
  tbb::parallel_for(tbb::blocked_range<LOOPITERTYPE>(0,numClasses,1),
                    [=,&ListOfClassStatistics](const tbb::blocked_range<LOOPITERTYPE> &r) {
                      for (LOOPITERTYPE iclass = r.begin(); iclass < r.end(); ++iclass) {
                        const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
                        const typename ByteImageType::ConstPointer currentCandidateRegion =
                            SubjectCandidateRegions[iclass].GetPointer();

                        // NOTE:  itk::Math:eps is too small itk::Math::eps;
                        CompensatedSummationType tmp_accumC = tbb::parallel_reduce(tbb::blocked_range3d<long>(0,size[2],1,
                                                                                       0,size[1],size[1]/2,
                                                                                       0,size[0],512),
                                               CompensatedSummationType(),
                                               [=](const tbb::blocked_range3d<long> &rng3d, CompensatedSummationType tmp) -> CompensatedSummationType {
                                                 for (long kk = rng3d.pages().begin(); kk < rng3d.pages().end(); ++kk) {
                                                   for (long jj = rng3d.rows().begin(); jj < rng3d.rows().end(); ++jj) {
                                                     for (long ii = rng3d.cols().begin(); ii < rng3d.cols().end(); ++ii) {
                                                       const typename TProbabilityImage::IndexType currIndex = {{ii, jj, kk}};
                                                       // Here pure plugs mask implicitly comes in! as CandidateRegions are multiplied by purePlugsMask!
                                                       if (currentCandidateRegion->GetPixel(currIndex)) {
                                                         const double currentProbValue = currentProbImage->GetPixel( currIndex);
                                                         tmp += currentProbValue;
                                                       }
                                                     }
                                                   }
                                                 }
                                                 return tmp;
                                               },
                                               [] ( CompensatedSummationType a,
                                                    const CompensatedSummationType & b ) ->  CompensatedSummationType {
                                                       a += b.GetSum();
                                                       return a;
                                               }
                        );
                        tmp_accumC += 1e-20;
                        ListOfClassStatistics[iclass].m_Weighting = tmp_accumC.GetSum();
                      }
                    });
  // Compute the means weighted by the probability of each value.
  tbb::parallel_for(tbb::blocked_range<LOOPITERTYPE>(0,numClasses,1),
                    [=,&ListOfClassStatistics](const tbb::blocked_range<LOOPITERTYPE> &r) {
                      for (LOOPITERTYPE iclass = r.begin(); iclass < r.end(); ++iclass) {
                        const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
                        const typename ByteImageType::ConstPointer currentCandidateRegion =
                            SubjectCandidateRegions[iclass].GetPointer();
                        ListOfClassStatistics[iclass].m_Means.clear();

                        for (typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
                             mapIt != InputImageMap.end(); ++mapIt) {
                          unsigned meanIndex(0);

                          ListOfClassStatistics[iclass].m_Means[mapIt->first] = 0.0;

                          for (typename InputImageVector::const_iterator imIt = mapIt->second.begin();
                               imIt != mapIt->second.end(); ++imIt, ++meanIndex) {
                            typename TInputImage::Pointer im1 = *imIt;
                            typename InputImageNNInterpolationType::Pointer im1Interp =
                                InputImageNNInterpolationType::New();
                            im1Interp->SetInputImage(im1);

                            const CompensatedSummationType muSumFinal =
                            tbb::parallel_reduce(tbb::blocked_range3d<long>(0,size[2],1,
                                                                            0,size[1],size[1]/2,
                                                                            0,size[0],512),
                            CompensatedSummationType(),
                            [=](const tbb::blocked_range3d<long> &rng, CompensatedSummationType muSum) -> CompensatedSummationType {
                              typename TProbabilityImage::PointType currPoint;
                              for (long kk = rng.pages().begin(); kk < rng.pages().end(); ++kk) {
                                for (long jj = rng.rows().begin(); jj < rng.rows().end(); ++jj) {
                                  for (long ii = rng.cols().begin(); ii < rng.cols().end(); ++ii) {
                                    const typename TProbabilityImage::IndexType currIndex = {{ii, jj, kk}};
                                    // transform probability image index to physical point
                                    PosteriorsList[0]->TransformIndexToPhysicalPoint(currIndex, currPoint);
                                    // Here pure plugs mask comes in, since CandidateRegions are multiplied by purePlugsMask!
                                    if (currentCandidateRegion->GetPixel(currIndex)) {
                                      const double currentProbValue = currentProbImage->GetPixel(currIndex);
                                      // input volumes may have a different voxel lattice than the probability image
                                      double currentInputValue = 1;
                                      if (im1Interp->IsInsideBuffer(currPoint)) {
                                        currentInputValue = im1Interp->Evaluate(currPoint);
                                      }
                                      if (logConvertValues) {
                                        muSum += currentProbValue * LOGP(currentInputValue);
                                      }
                                      else {
                                        muSum += currentProbValue * (currentInputValue);
                                      }
                                    }
                                  }
                                }
                              }
                              return muSum;
                            }
                              ,
                            [] ( CompensatedSummationType a,
                                 const CompensatedSummationType & b ) ->  CompensatedSummationType {
                                   a += b.GetSum();
                                   return a;
                            }
                              );
                            const double mymean=(muSumFinal.GetSum()) / ListOfClassStatistics[iclass].m_Weighting;
                            ListOfClassStatistics[iclass].m_Means[mapIt->first] += mymean;
                          }
                          // averaging the means of all images of this image modality
                          ListOfClassStatistics[iclass].m_Means[mapIt->first] /= mapIt->second.size();
                        }
                      }
                    });

  // for each prior (posterior) class, different means are computed for each modality channel.

  ////////////////////////////////
  // Now compute covariance matrix( numOfModalities x numOfModalities )
  // e.g. A 2x2 matrix if only T1 and T2 modality channels are involved.
  // Note that we can have several T1s and several T2 images.

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
  tbb::parallel_for(tbb::blocked_range<LOOPITERTYPE>(0,numClasses,1),
                    [=,&ListOfClassStatistics](const tbb::blocked_range<LOOPITERTYPE> &r) {
                      for (LOOPITERTYPE iclass = r.begin(); iclass < r.end(); ++iclass) {
                        const typename TProbabilityImage::ConstPointer currentProbImage = PosteriorsList[iclass].GetPointer();
                        const typename ByteImageType::ConstPointer currentCandidateRegion =
                            SubjectCandidateRegions[iclass].GetPointer();
                        //
                        // this will end up as a vnl_matrix for assignment to
                        // the Class Statistics object after this is computed.
                        orderedmap<std::string, orderedmap<std::string, double> > TypeCovariance;
                        // initialize -- no easy way since it is a map of maps
                        for (typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
                             mapIt != InputImageMap.end(); ++mapIt) {
                          for (typename MapOfInputImageVectors::const_iterator mapIt2 = InputImageMap.begin();
                               mapIt2 != InputImageMap.end(); ++mapIt2) {
                            TypeCovariance[mapIt->first][mapIt2->first] = 0.0;
                          }
                        }
                        //
                        // compute per-Image Type covariance
                        for (typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
                             mapIt != InputImageMap.end(); ++mapIt) {
                          const double mu1 =
                              ListOfClassStatistics[iclass].m_Means[mapIt->first];

                          for (unsigned i = 0; i < mapIt->second.size(); ++i) {
                            typename TInputImage::Pointer im1 = mapIt->second[i];
                            typename InputImageNNInterpolationType::Pointer im1Interp =
                                InputImageNNInterpolationType::New();
                            im1Interp->SetInputImage(im1);

                            bool first_through_inner_loop(true);

                            for (typename MapOfInputImageVectors::const_iterator mapIt2 = mapIt;
                                 mapIt2 != InputImageMap.end(); ++mapIt2) {
                              size_t j = 0;
                              if (first_through_inner_loop) {
                                j = i;
                                first_through_inner_loop = false;
                              }
                              const double mu2 =
                                  ListOfClassStatistics[iclass].m_Means[mapIt2->first];
                              for (; j < mapIt2->second.size(); ++j) {
                                typename TInputImage::Pointer im2 = mapIt2->second[j];
                                typename InputImageNNInterpolationType::Pointer im2Interp =
                                    InputImageNNInterpolationType::New();
                                im2Interp->SetInputImage(im2);

                                CompensatedSummationType reduced_varC = tbb::parallel_reduce
                                    (tbb::blocked_range3d<long>(0,size[2],1,
                                                                0,size[1],size[1]/2,
                                                                0,size[0],512),
                                                     CompensatedSummationType(), /*Initial value of reduction */
                                                     [=](const tbb::blocked_range3d<long> &rng, CompensatedSummationType var )
                                                                        -> CompensatedSummationType {
                                                       for (long kk = rng.pages().begin(); kk < rng.pages().end(); ++kk) {
                                                         for (long jj = rng.rows().begin(); jj < rng.rows().end(); ++jj) {
                                                           for (long ii = rng.cols().begin(); ii < rng.cols().end(); ++ii) {
                                                             const typename TProbabilityImage::IndexType currIndex = {{ii, jj, kk}};
                                                             // transform probability image index to physical point
                                                             typename TProbabilityImage::PointType currPoint;
                                                             PosteriorsList[0]->TransformIndexToPhysicalPoint(currIndex,
                                                                                                              currPoint);
                                                             // Here pure plugs mask comes in, since CandidateRegions are multiplied by purePlugsMask!
                                                             if (currentCandidateRegion->GetPixel(currIndex)) {
                                                               const double currentProbValue = currentProbImage->GetPixel(
                                                                   currIndex);
                                                               // input image values should be evaluated in physical space.
                                                               double inputValue1 = 1;
                                                               double inputValue2 = 1;
                                                               if (im1Interp->IsInsideBuffer(currPoint)) {
                                                                 inputValue1 = im1Interp->Evaluate(currPoint);
                                                               }
                                                               if (im2Interp->IsInsideBuffer(currPoint)) {
                                                                 inputValue2 = im2Interp->Evaluate(currPoint);
                                                               }

                                                               if (logConvertValues) {
                                                                 const double diff1 = LOGP(inputValue1) - mu1;
                                                                 const double diff2 = LOGP(inputValue2) - mu2;
                                                                 var += currentProbValue * (diff1 * diff2);
                                                               }
                                                               else {
                                                                 const double diff1 = inputValue1 - mu1;
                                                                 const double diff2 = inputValue2 - mu2;
                                                                 var += currentProbValue * (diff1 * diff2);
                                                               }
                                                             }
                                                           }
                                                         }
                                                       }
                                                       return var;
                                                     },
                                                     /* Reduction Operator */
                            [] ( CompensatedSummationType a,
                                 const CompensatedSummationType & b ) ->  CompensatedSummationType {
                                   a += b.GetSum();
                                   return a;
                            }
                                    );
                                double reduced_var = reduced_varC.GetSum() / ListOfClassStatistics[iclass].m_Weighting;

                                // Adjust diagonal, to make sure covariance is pos-def
                                if (mapIt == mapIt2 && i == j) {
                                  reduced_var += 1e-20;
                                }
                                // accumulate
                                TypeCovariance[mapIt->first][mapIt2->first] += reduced_var;
                                TypeCovariance[mapIt2->first][mapIt->first] += reduced_var;
                              }
                            }
                          }
                        }
                        // above loop accumulates covariances per type
                        // now divide out # of averaged variances
                        // and copy to vnl matrix
                        MatrixType covtmp(numModalities, numModalities, 0.0);
                        unsigned int i = 0;
                        for (typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
                             mapIt != InputImageMap.end(); ++mapIt, ++i) {
                          unsigned int j = 0;
                          for (typename MapOfInputImageVectors::const_iterator mapIt2 = InputImageMap.begin();
                               mapIt2 != InputImageMap.end(); ++mapIt2, ++j) {
                            covtmp(i, j) = TypeCovariance[mapIt->first][mapIt2->first] /
                                           static_cast<double>(mapIt->second.size() * mapIt2->second.size());
                          }
                        }
                        ListOfClassStatistics[iclass].m_Covariance = covtmp;
                      } // end covariance loop
                    });

  if( DebugLevel > 9 )
    {
    std::cout << "=================================================" << std::endl;
    for( LOOPITERTYPE iclass = 0; iclass < (LOOPITERTYPE)numClasses; iclass++ )
      {
      unsigned ichan = 0;
      for(typename MapOfInputImageVectors::const_iterator mapIt = InputImageMap.begin();
          mapIt != InputImageMap.end(); ++mapIt)
        {
        muLogMacro( << "DEBUG MEAN (channel " << ichan << ", class " << iclass << "): \n"
                    << ListOfClassStatistics[iclass].m_Means[mapIt->first]
                    << " \n" << std::endl );
        ++ichan;
        }
      muLogMacro( << "DEBUG Covariances (class " << iclass << "):\n" << ListOfClassStatistics[iclass].m_Covariance
                  << std::endl );
      }
    }
}

#endif // __ComputeDistributions__h__
