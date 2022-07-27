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
#ifndef __LLSBiasCorrector_hxx
#define __LLSBiasCorrector_hxx

#include <cfloat>
#include <iostream>

#include "Log.h"
#include "StandardizeMaskIntensity.h"
#include "LLSBiasCorrector.h"
#include "itkMath.h"
#include "itkTimeProbe.h"
#include "ComputeDistributions.h"

#include "tbb/blocked_range2d.h"
#include "tbb/parallel_reduce.h"


#define USE_HALF_RESOLUTION 1
#define MIN_SKIP_SIZE 2

// Use the normal equation? Less accurate, but requires much less memory
#define LLSBIAS_USE_NORMAL_EQUATION 1
//
//
// //////////////////////////////////////////////////////////////////////////////
#define mypow(a, b) std::pow((a), static_cast<double>(b))
//
//
// //////////////////////////////////////////////////////////////////////////////

template <typename TInputImage, typename TProbabilityImage>
LLSBiasCorrector<TInputImage, TProbabilityImage>::LLSBiasCorrector()
{
  m_InputImages.clear();

  m_MaxDegree = 4;

  m_DebugLevel = 0;
  m_OutputDebugDir = "";
  m_SampleSpacing = 4.0;

  m_XMu[0] = 0.0;
  m_XMu[1] = 0.0;
  m_XMu[2] = 0.0;

  m_XStd[0] = 1.0;
  m_XStd[1] = 1.0;
  m_XStd[2] = 1.0;
}

template <typename TInputImage, typename TProbabilityImage>
LLSBiasCorrector<TInputImage, TProbabilityImage>::~LLSBiasCorrector()
{
  m_ForegroundBrainMask = nullptr;
}

template <typename TInputImage, typename TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>::CheckInputs()
{
  // if (m_MaxDegree == 0)
  // itkExceptionMacro(<< "Max bias degree is zero" << std::endl );

  if (m_InputImages.empty())
  {
    itkExceptionMacro(<< "No input image specified" << std::endl);
  }

  if (this->GetFirstInputImage()->GetImageDimension() != 3)
  {
    itkExceptionMacro(<< "Input dimension invalid: only supports 3D images" << std::endl);
  }

  if (m_BiasPosteriors.empty())
  {
    itkExceptionMacro(<< "Must have one or more class probabilities" << std::endl);
  }

  const InputImageSizeType size = this->GetFirstInputImage()->GetLargestPossibleRegion().GetSize();
  for (unsigned int i = 1; i < m_InputImages.size(); i++)
  {
    const InputImageSizeType size_i = this->GetFirstInputImage()->GetLargestPossibleRegion().GetSize();
    if (size != size_i)
    {
      itkExceptionMacro(<< "Image sizes do not match" << std::endl);
    }
  }
  for (unsigned int i = 0; i < m_BiasPosteriors.size(); i++)
  {
    if (m_BiasPosteriors[i]->GetImageDimension() != 3)
    {
      itkExceptionMacro(<< "Probability [" << i << "] has invalid dimension: only supports 3D images" << std::endl);
    }
    const ProbabilityImageSizeType psize = m_BiasPosteriors[i]->GetLargestPossibleRegion().GetSize();
    if (size[0] != psize[0] || size[1] != psize[1] || size[2] != psize[2])
    {
      itkExceptionMacro(<< "Image data and probabilities 3D size mismatch" << std::endl);
    }
  }
}

template <typename TInputImage, typename TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>::ComputeDistributions()
{
  muLogMacro(<< "LLSBiasCorrector: Computing means and variances..." << std::endl);
  this->m_ListOfClassStatistics = CombinedComputeDistributions<TInputImage, TProbabilityImage, MatrixType>(
    this->m_CandidateRegions, m_InputImages, m_BiasPosteriors, this->m_DebugLevel, true);
}

template <typename TInputImage, typename TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>::SetMaxDegree(unsigned int n)
{
  muLogMacro(<< "SetMaxDegree" << std::endl);

  m_MaxDegree = n;

  // NOTE:  update basis equations, Need to call the Initialize in order to
  // update
  // the basis equations.
  if (!m_ForegroundBrainMask.IsNull())
  {
    this->Initialize();
  }
#if 1 // HACK
  if (!m_BiasPosteriors.empty())
  {
    this->SetProbabilities(m_BiasPosteriors, this->m_CandidateRegions);
  }
#endif
}

template <typename TInputImage, typename TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>::SetSampleSpacing(double s)
{
  muLogMacro(<< "SetSampleSpacing" << std::endl);

  m_SampleSpacing = s;

#if 1 // HACK
  if (!m_BiasPosteriors.empty())
  {
    this->SetProbabilities(m_BiasPosteriors, this->m_CandidateRegions);
  }
#endif
  // NOTE:  update basis equations, Need to call the Initialize in order to
  // update
  // the basis equations.
  if (m_ForegroundBrainMask.IsNotNull())
  {
    this->Initialize();
  }
}

template <typename TInputImage, typename TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>::SetForegroundBrainMask(ByteImageType * mask)
{
  m_ForegroundBrainMask = mask;
  this->Initialize();
}

template <typename TInputImage, typename TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>::Initialize()
{
  const InputImageSizeType size = m_ForegroundBrainMask->GetLargestPossibleRegion().GetSize();

  // Compute skips along each dimension
#ifdef USE_HALF_RESOLUTION // Need to do full scale to prevent checkerboard
  // images from being created
  const InputImageSpacingType spacing = m_ForegroundBrainMask->GetSpacing();
  unsigned int                skips[3];
  skips[0] = static_cast<int>(m_SampleSpacing / spacing[0]);
  skips[1] = static_cast<int>(m_SampleSpacing / spacing[1]);
  skips[2] = static_cast<int>(m_SampleSpacing / spacing[2]);

  if (skips[0] < MIN_SKIP_SIZE)
  {
    skips[0] = MIN_SKIP_SIZE;
  }
  if (skips[1] < MIN_SKIP_SIZE)
  {
    skips[1] = MIN_SKIP_SIZE;
  }
  if (skips[2] < MIN_SKIP_SIZE)
  {
    skips[2] = MIN_SKIP_SIZE;
  }
  muLogMacro(<< "Sample skips: " << skips[0] << " x " << skips[1] << " x " << skips[2] << std::endl);
#else
  const unsigned int skips[3] = { 1, 1, 1 };
#endif

  const unsigned int numCoefficients = (m_MaxDegree + 1) * (m_MaxDegree + 2) / 2 * (m_MaxDegree + 3) / 3;

  // Number of pixels with non-zero weights, downsampled
  unsigned numEquations = 0;
  // Assume that only 0.125 of the image is part of mask
  m_ValidIndicies.reserve(size[2] * size[1] * size[0] / (8 * skips[0] * skips[1] * skips[2]));

  m_ValidIndicies.resize(0);
  {
    // Not parallizable! ORDER IS IMPORTANT  //#pragma omp parallel for
    // reduction(+:numEquations) default(shared)
    for (long kk = 0; kk < static_cast<long>(size[2]); kk += skips[2])
    {
      for (long jj = 0; jj < static_cast<long>(size[1]); jj += skips[1])
      {
        for (long ii = 0; ii < static_cast<long>(size[0]); ii += skips[0])
        {
          const ProbabilityImageIndexType currProbIndex = { { ii, jj, kk } };
          if (m_ForegroundBrainMask->GetPixel(currProbIndex) != 0)
          {
            m_ValidIndicies.push_back(currProbIndex);
            // numEquations++;
          }
        }
      }
    }
  }
  numEquations = m_ValidIndicies.size();
  muLogMacro(<< "Linear system size = " << numEquations << " x " << numCoefficients << std::endl);

  // Make sure that number of equations >= number of unknowns
  if (numEquations < numCoefficients)
  {
    std::cout << size << std::endl;
    itkExceptionMacro(<< "Number of unknowns exceed number of equations:" << numEquations << " < " << numCoefficients);
  }

  // Create basis matrix

  muLogMacro(<< "Computing polynomial basis functions..." << std::endl);

  m_Basis.set_size(numEquations, numCoefficients);

  using IterType = typename std::vector<ProbabilityImageIndexType>::const_iterator;
  {
    // Coordinate scaling and offset parameters
    vnl_vector_fixed<unsigned long long int, 3> local_XMu = tbb::parallel_reduce(
      tbb::blocked_range<IterType>(m_ValidIndicies.begin(), m_ValidIndicies.end(), 1),

      vnl_vector_fixed<unsigned long long int, 3>(),

      [](const tbb::blocked_range<IterType> &        r,
         vnl_vector_fixed<unsigned long long int, 3> init) -> vnl_vector_fixed<unsigned long long int, 3> {
        for (auto currIndex = r.begin(); currIndex != r.end(); ++currIndex)
        {
          init[0] = init[0] + (*currIndex)[0];
          init[1] = init[1] + (*currIndex)[1];
          init[2] = init[2] + (*currIndex)[2];
        }
        return init;
      },

      [](const vnl_vector_fixed<unsigned long long int, 3> & a, const vnl_vector_fixed<unsigned long long int, 3> & b)
        -> vnl_vector_fixed<unsigned long long int, 3> { return a + b; });

    {
      local_XMu += tbb::parallel_reduce(
        tbb::blocked_range<IterType>(m_ValidIndicies.begin(), m_ValidIndicies.end(), 1),
        vnl_vector_fixed<unsigned long long int, 3>(),
        [=](const tbb::blocked_range<IterType> &        r,
            vnl_vector_fixed<unsigned long long int, 3> newLocal_XMu) -> vnl_vector_fixed<unsigned long long int, 3> {
          for (auto currProbIndex = r.begin(); currProbIndex < r.end(); ++currProbIndex)
          {
            newLocal_XMu[0] += (*currProbIndex)[0];
            newLocal_XMu[1] += (*currProbIndex)[1];
            newLocal_XMu[2] += (*currProbIndex)[2];
          }
          return newLocal_XMu;
        },
        [](const vnl_vector_fixed<unsigned long long int, 3> & a, const vnl_vector_fixed<unsigned long long int, 3> & b)
          -> vnl_vector_fixed<unsigned long long int, 3> { return a + b; });
    }
    const double invNumEquations = 1.0 / static_cast<double>(numEquations);
    m_XMu[0] = static_cast<double>(local_XMu[0]) * invNumEquations;
    m_XMu[1] = static_cast<double>(local_XMu[1]) * invNumEquations;
    m_XMu[2] = static_cast<double>(local_XMu[2]) * invNumEquations;
  }

  {
    const std::vector<CompensatedSummationType> local_XStd_final = tbb::parallel_reduce(
      tbb::blocked_range<IterType>(m_ValidIndicies.begin(), m_ValidIndicies.end(), 1),
      std::vector<CompensatedSummationType>(),
      [=](const tbb::blocked_range<IterType> &  rng,
          std::vector<CompensatedSummationType> local_XStd) -> std::vector<CompensatedSummationType> {
        local_XStd.resize(3);
        for (auto currProbIndex = rng.begin(); currProbIndex < rng.end(); ++currProbIndex)
        {
          const double diff0 = static_cast<double>((*currProbIndex)[0]) - m_XMu[0];
          local_XStd[0] += diff0 * diff0;
          const double diff1 = static_cast<double>((*currProbIndex)[1]) - m_XMu[1];
          local_XStd[1] += diff1 * diff1;
          const double diff2 = static_cast<double>((*currProbIndex)[2]) - m_XMu[2];
          local_XStd[2] += diff2 * diff2;
        }
        return local_XStd;
      },
      [](const std::vector<CompensatedSummationType> & a,
         const std::vector<CompensatedSummationType> & b) -> std::vector<CompensatedSummationType> {
        std::vector<CompensatedSummationType> c(a);
        c[0] += b[0].GetSum();
        c[1] += b[1].GetSum();
        c[2] += b[2].GetSum();
        return c;
      });
    m_XStd[0] = std::sqrt(local_XStd_final[0].GetSum() / numEquations);
    m_XStd[1] = std::sqrt(local_XStd_final[1].GetSum() / numEquations);
    m_XStd[2] = std::sqrt(local_XStd_final[2].GetSum() / numEquations);
  }

  {
    // Row and column indices
    // Fill in polynomial basis values
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numEquations, 1),
                      [=](tbb::blocked_range<unsigned int> & rng) {
                        for (unsigned int r = rng.begin(); r < rng.end(); ++r)
                        {
                          const ProbabilityImageIndexType & currProbIndex = m_ValidIndicies[r];
                          unsigned int                      c = 0;
                          for (unsigned int order = 0; order <= m_MaxDegree; order++)
                          {
                            for (unsigned int xorder = 0; xorder <= order; xorder++)
                            {
                              for (unsigned int yorder = 0; yorder <= (order - xorder); yorder++)
                              {
                                const int zorder = order - xorder - yorder;

                                const double xc = (currProbIndex[0] - m_XMu[0]) / m_XStd[0];
                                const double yc = (currProbIndex[1] - m_XMu[1]) / m_XStd[1];
                                const double zc = (currProbIndex[2] - m_XMu[2]) / m_XStd[2];

                                m_Basis(r, c) = mypow(xc, xorder) * mypow(yc, yorder) * mypow(zc, zorder);
                                ++c;
                              }
                            }
                          }
                        }
                      });
  }
}

template <typename TInputImage, typename TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>::SetProbabilities(
  const std::vector<ProbabilityImagePointer> &         probs,
  const std::vector<typename ByteImageType::Pointer> & candidateRegions)
{
  muLogMacro(<< "SetProbabilities" << std::endl);

  if (probs.empty())
  {
    itkExceptionMacro(<< "Need one or more probabilities" << std::endl);
  }
  if (probs.size() != candidateRegions.size())
  {
    itkExceptionMacro(<< "Probabilities and CandidateRegions size must be the same." << std::endl);
  }
  for (unsigned int i = 0; i < probs.size(); i++)
  {
    if (probs[i].IsNull())
    {
      itkExceptionMacro(<< "One of input probabilities not initialized" << std::endl);
    }
  }
  for (unsigned int i = 0; i < candidateRegions.size(); i++)
  {
    if (candidateRegions[i].IsNull())
    {
      itkExceptionMacro(<< "One of input candidate regions is not initialized" << std::endl);
    }
  }
  m_BiasPosteriors = probs;
  m_CandidateRegions = candidateRegions;
}

template <typename TInputImage, typename TProbabilityImage>
typename LLSBiasCorrector<TInputImage, TProbabilityImage>::MapOfInputImageVectors
LLSBiasCorrector<TInputImage, TProbabilityImage>::CorrectImages(const unsigned int CurrentIterationID)
{
  muLogMacro(<< "\n*** Correct Images ***" << std::endl);
  itk::TimeProbe CorrectImagesTimer;
  CorrectImagesTimer.Start();
  // Verify input
  this->CheckInputs();

  // Compute means and variances
  this->ComputeDistributions();

  unsigned int numModalities = this->m_InputImages.size();

  const unsigned int numClasses = m_BiasPosteriors.size();

  /* if m_MaxDegree = 4/3/2/1, then this is 35/20/10/4 */
  const unsigned int numCoefficients = (m_MaxDegree + 1) * (m_MaxDegree + 2) / 2 * (m_MaxDegree + 3) / 3;

  muLogMacro(<< numClasses << " classes" << std::endl);
  muLogMacro(<< numCoefficients << " coefficients" << std::endl);

  /* compute inverse matrix for each tissue type */
  muLogMacro(<< "Computing inverse covars...\n" << std::endl);
  std::vector<MatrixType> invCovars;
  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  {
    std::cout << this->m_ListOfClassStatistics[iclass].m_Covariance << std::endl;
    MatrixInverseType inverse_temp(this->m_ListOfClassStatistics[iclass].m_Covariance);
    MatrixType        temp = inverse_temp.as_matrix();
    invCovars.push_back(temp);
  }

  // Create matrices and vectors
  // lhs = replicated basis polynomials for each channel, weighted by inv cov
  // rhs = difference image between original and reconstructed mean image

  muLogMacro(<< "Creating matrices for LLS..." << std::endl);

  const unsigned int numEquations = m_Basis.rows(); /*  since m_Basis is [numEqu *
                                                        35 ] matrix */

  muLogMacro(<< numEquations << " equations, " << numCoefficients << " coefficients" << std::endl);

  // Compute  orthogonal transpose component of basis
  muLogMacro(<< "Computing ortho part of basis" << std::endl);
#if 1
  // Note: vnl_qr gives Q mxm and R mxn for A mxn

  MatrixType basisT;
  {
    MatrixQRType qr(m_Basis);

    // Get economy size R (square)
    MatrixType R(numCoefficients, numCoefficients, 0);

    {
      MatrixType Rfull = qr.R(); /* right triangular matrix */
      for (unsigned int r = 0; r < numCoefficients; r++)
      {
        for (unsigned int c = r; c < numCoefficients; c++)
        {
          R(r, c) = Rfull(r, c);
        }
      }
      // Rfull.set_size(1, 1);
    }
    // NOTE: to get mxn Q from vnl_qr, Q'*Q = id nxn
    basisT = m_Basis * MatrixInverseType(R).as_matrix();
  }
  basisT.inplace_transpose(); // basisT = Q'
#else
  // Do this instead for ordinary weighted LSQ
  MatrixType basisT = m_Basis.transpose();
#endif

#if LLSBIAS_USE_NORMAL_EQUATION
  MatrixType lhs(numCoefficients * numModalities, numCoefficients * numModalities);
  MatrixType rhs(numCoefficients * numModalities, 1);
#else
  MatrixType lhs(numEquations * numModalities, numCoefficients * numModalities);
  MatrixType rhs(numEquations * numModalities, 1);
#endif

  muLogMacro(<< "Fill rhs" << std::endl);
  rhs.fill(0.0);

  // Compute ratio between original and flat image, weighted using posterior
  // probability and inverse covariance
  {
    unsigned int modality1 = 0;
    for (typename MapOfInputImageVectors::const_iterator mapIt = this->m_InputImages.begin();
         mapIt != this->m_InputImages.end();
         ++mapIt, ++modality1)
    {
      MatrixType R_i(numEquations, 1, 0.0);

      unsigned int modality2 = 0;
      for (typename MapOfInputImageVectors::const_iterator mapIt2 = this->m_InputImages.begin();
           mapIt2 != this->m_InputImages.end();
           ++mapIt2, ++modality2)
      {
        unsigned int numCurModalityImages = mapIt2->second.size();
        for (unsigned int imIndex = 0; imIndex < numCurModalityImages; ++imIndex)
        {
          typename InputImageNNInterpolationType::Pointer inputImageInterp = InputImageNNInterpolationType::New();
          inputImageInterp->SetInputImage(mapIt2->second[imIndex].GetPointer());
          tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numEquations, 1),
                            [=, &R_i](const tbb::blocked_range<unsigned int> & r) {
                              for (unsigned int eq = r.begin(); eq < r.end(); eq++)
                              {
                                const ProbabilityImageIndexType & currProbIndex = m_ValidIndicies[eq];
                                // Compute reconstructed intensity, weighted by prob * invCov
                                double sumW = DBL_EPSILON;
                                double recon = 0;
                                for (unsigned int iclass = 0; iclass < numClasses; iclass++)
                                {
                                  const MatrixType & invCov = invCovars[iclass];

                                  const double w =
                                    m_BiasPosteriors[iclass]->GetPixel(currProbIndex) * invCov(modality1, modality2);
                                  sumW += w;
                                  recon += w * this->m_ListOfClassStatistics[iclass].m_Means[mapIt2->first];
                                }
                                recon /= sumW;

                                // transform probability image index to physical point
                                typename ProbabilityImageType::PointType currProbPoint;
                                m_BiasPosteriors[0]->TransformIndexToPhysicalPoint(currProbIndex, currProbPoint);

                                typename InputImageNNInterpolationType::OutputType inputImageValue =
                                  1; // default value must be 1
                                if (inputImageInterp->IsInsideBuffer(currProbPoint))
                                {
                                  inputImageValue = inputImageInterp->Evaluate(currProbPoint);
                                }

                                const double bias = LOGP(inputImageValue) - recon;
                                //,&R_,&R_,&R_iii divide by # of images of current modality -- in essence
                                // you're averaging them.
                                R_i(eq, 0) += (sumW * bias) / numCurModalityImages;
                              }
                            });
        }
      } // for jchan

#if LLSBIAS_USE_NORMAL_EQUATION
      R_i = basisT * R_i;
      for (unsigned int row = 0; row < numCoefficients; row++)
      {
        rhs(modality1 * numCoefficients + row, 0) = R_i(row, 0);
      }
#else
      for (unsigned int row = 0; row < numEquations; row++)
      {
        rhs(modality1 * numEquations + row, 0) = R_i(row, 0);
      }
#endif
    }
  }

  muLogMacro(<< "Fill lhs" << std::endl);

  // Compute LHS using replicated basis entries, weighted using posterior
  // probability and inverse covariance
  tbb::parallel_for(tbb::blocked_range2d<LOOPITERTYPE>(0, numModalities, 0, numModalities),
                    [=, &lhs](const tbb::blocked_range2d<LOOPITERTYPE> & r) {
                      for (LOOPITERTYPE ichan = r.rows().begin(); ichan < r.rows().end(); ichan++)
                      {
                        for (LOOPITERTYPE jchan = r.cols().begin(); jchan < r.cols().end(); jchan++)
                        {
                          MatrixType Wij_A(numEquations, numCoefficients, 0.0);
                          {
                            for (LOOPITERTYPE eq = 0; eq < numEquations; eq++)
                            {
                              const ProbabilityImageIndexType & currProbIndex = m_ValidIndicies[eq];
                              double                            sumW = DBL_EPSILON;
                              for (unsigned int iclass = 0; iclass < numClasses; iclass++)
                              {
                                const MatrixType & invCov = invCovars[iclass];
                                double w = m_BiasPosteriors[iclass]->GetPixel(currProbIndex) * invCov(ichan, jchan);
                                sumW += w;
                              }
                              for (unsigned int col = 0; col < numCoefficients; col++)
                              {
                                Wij_A(eq, col) = sumW * m_Basis(eq, col);
                              }
                            }
                          }
#if LLSBIAS_USE_NORMAL_EQUATION
                          MatrixType lhs_ij = basisT * Wij_A;
                          for (unsigned int row = 0; row < numCoefficients; row++)
                          {
                            for (unsigned int col = 0; col < numCoefficients; col++)
                            {
                              lhs(row + ichan * numCoefficients, col + jchan * numCoefficients) = lhs_ij(row, col);
                            }
                          }
#else
                            for( unsigned int row = 0; row < numEquations; row++ )
                              {
                              for( unsigned int col = 0; col < numCoefficients; col++ )
                                {
                                lhs(row + ichan * numEquations, col + jchan * numCoefficients)
                                  = Wij_A(row, col);
                                }
                              }
#endif
                        } // for jchan
                      }   // for ichan
                    });

  muLogMacro(<< "Solve " << lhs.rows() << " x " << lhs.columns() << std::endl);

  // Use VNL to solve linear system
  MatrixType coeffs;
  {
    // MatrixQRType qr(lhs);
    // coeffs = qr.solve(rhs);
    // SVD more expensive, should be more accurate
    MatrixSVDType svd(lhs);

    svd.zero_out_absolute(1e-8);
    coeffs = svd.solve(rhs);
  }
#ifndef WIN32
  if (!std::isfinite(static_cast<double>(coeffs[0][0])))
#else
  if (coeffs[0][0] != std::numeric_limits::infinity())
#endif
  {
    itkExceptionMacro(<< "\ncoeffs: \n"
                      << coeffs
                      // << "\nlhs_ij: \n" << lhs_ij
                      << "\nbasisT: \n"
                      << basisT
                      // << "\nWij_A: \n" << Wij_A
                      << "\nlhs: \n"
                      << lhs << "\nrhs: \n"
                      << rhs);
  }
  // Clear memory for the basis transpose
  basisT.set_size(0, 0);

  if (this->m_DebugLevel > 9)
  {
    muLogMacro(<< "Bias field coeffs after LLS:" << std::endl << coeffs);
  }

  // Remove bias
  muLogMacro(<< "Correcting input images..." << std::endl);

  MapOfInputImageVectors outputs;
  /*
   * Each bias corrected output image has the same physical/voxel space as its corresponding input image.
   * Note that input images may have different voxel lattices although they are all aligned in physical space.
   * Also, notice that brain mask and posteriors are at the same voxel lattice of the very first image of
     input images map (the first image of the first modality channel), so output image "i" may have a different
     voxel space than the brain mask and posteriors.
   */
  {
    typename MaskNNInterpolationType::Pointer foregroundBrainMaskInterp = MaskNNInterpolationType::New();
    foregroundBrainMaskInterp->SetInputImage(this->m_ForegroundBrainMask.GetPointer());

    typename MaskNNInterpolationType::Pointer allTissueMaskInterp = MaskNNInterpolationType::New();
    allTissueMaskInterp->SetInputImage(this->m_AllTissueMask.GetPointer());

    unsigned int ichan = 0;
    for (typename MapOfInputImageVectors::const_iterator mapIt = this->m_InputImages.begin();
         mapIt != this->m_InputImages.end();
         ++mapIt, ++ichan)
    {
      for (auto imIt = mapIt->second.begin(); imIt != mapIt->second.end(); ++imIt)
      {
        InternalImagePointer curOutput = InternalImageType::New();
        curOutput->CopyInformation((*imIt));
        curOutput->SetRegions((*imIt)->GetLargestPossibleRegion());
        curOutput->Allocate();
        curOutput->FillBuffer(0);

        // Compute the LOGP transformed bias field
        InternalImagePointer biasIntensityScaleFactor = InternalImageType::New();
        biasIntensityScaleFactor->CopyInformation(curOutput);
        biasIntensityScaleFactor->SetRegions(curOutput->GetLargestPossibleRegion());
        biasIntensityScaleFactor->Allocate();
        biasIntensityScaleFactor->FillBuffer(0);

        double maxBiasInForegroundMask = std::numeric_limits<double>::min();
        double minBiasInForegroundMask = std::numeric_limits<double>::max();

        const InputImageSizeType outsize = curOutput->GetLargestPossibleRegion().GetSize();
        tbb::parallel_for(tbb::blocked_range3d<long>(0, outsize[2], 0, outsize[1], 0, outsize[0]),
                          [=, &maxBiasInForegroundMask, &minBiasInForegroundMask](tbb::blocked_range3d<long> & r) {
                            for (long kk = r.pages().begin(); kk < r.pages().end(); ++kk)
                            {
                              for (long jj = r.rows().begin(); jj < r.rows().end(); ++jj)
                              {
                                for (long ii = r.cols().begin(); ii < r.cols().end(); ++ii)
                                {
                                  const InternalImageIndexType currOutIndex = { { ii, jj, kk } }; // index of currOutput
                                  // Masks and probability image should be evaluated in physical space, since
                                  // they may not have the same voxel lattice as the current output image.
                                  typename InternalImageType::PointType currOutPoint;
                                  curOutput->TransformIndexToPhysicalPoint(currOutIndex, currOutPoint);

                                  double       logFitValue = 0.0;
                                  unsigned int c = ichan * numCoefficients;
                                  for (unsigned int order = 0; order <= m_MaxDegree; order++)
                                  {
                                    for (unsigned int xorder = 0; xorder <= order; xorder++)
                                    {
                                      for (unsigned int yorder = 0; yorder <= (order - xorder); yorder++)
                                      {
                                        const int zorder = order - xorder - yorder;

                                        const double xc = (currOutIndex[0] - m_XMu[0]) / m_XStd[0];
                                        const double yc = (currOutIndex[1] - m_XMu[1]) / m_XStd[1];
                                        const double zc = (currOutIndex[2] - m_XMu[2]) / m_XStd[2];

                                        const double poly = mypow(xc, xorder) * mypow(yc, yorder) * mypow(zc, zorder);

                                        // logFitValue += coeffs[c] * poly;
                                        logFitValue += coeffs(c, 0) * poly;
                                        ++c;
                                      }
                                    }
                                  }

                                  ByteImagePixelType maskValue = 0;
                                  if (foregroundBrainMaskInterp->IsInsideBuffer(currOutPoint))
                                  {
                                    maskValue = foregroundBrainMaskInterp->Evaluate(currOutPoint);
                                  }
                                  /* NOTE:  For regions listed as background, clamp the outputs[ichan
                                   */
                                  if (std::isnan(logFitValue) || std::isinf(logFitValue))
                                  {
                                    std::cout << "WARNING:  Bad Scale Value Computed!" << std::endl;
                                    logFitValue = 0.0;
                                  }
                                  const double denom = EXPP(logFitValue);
                                  const double multiplicitiveBiasCorrectionFactor = 1.0 / denom;
                                  if (maskValue != 0)
                                  {
                                    if (multiplicitiveBiasCorrectionFactor > maxBiasInForegroundMask)
                                    {
                                      maxBiasInForegroundMask = multiplicitiveBiasCorrectionFactor;
                                    }
                                    if (multiplicitiveBiasCorrectionFactor < minBiasInForegroundMask)
                                    {
                                      minBiasInForegroundMask = multiplicitiveBiasCorrectionFactor;
                                    }
                                  }
                                  biasIntensityScaleFactor->SetPixel(
                                    currOutIndex, static_cast<InternalImagePixelType>(multiplicitiveBiasCorrectionFactor));
                                } // for currOutIndex[0]
                              }
                            }
                          });
        std::cout << "Foreground Mask Bias Correction MIN: " << minBiasInForegroundMask
                  << " MAX: " << maxBiasInForegroundMask << std::endl;
        // Correct image using (clamped) bias field)
        tbb::parallel_for(
          tbb::blocked_range3d<long>(0, outsize[2], 0, outsize[1], 0, outsize[0]), [=](tbb::blocked_range3d<long> & r) {
            for (auto kk = r.pages().begin(); kk < r.pages().end(); ++kk)
            {
              for (auto jj = r.rows().begin(); jj < r.rows().end(); ++jj)
              {
                for (auto ii = r.cols().begin(); ii < r.cols().end(); ++ii)
                {
                  const InternalImageIndexType          currOutIndex = { { ii, jj, kk } }; // index of currOutput
                  typename InternalImageType::PointType currOutPoint;
                  curOutput->TransformIndexToPhysicalPoint(currOutIndex, currOutPoint);

                  double multiplicitiveBiasCorrectionFactor = biasIntensityScaleFactor->GetPixel(currOutIndex);
                  if (multiplicitiveBiasCorrectionFactor > maxBiasInForegroundMask) //
                                                                                    //
                                                                                    // CLAMP
                  {
                    multiplicitiveBiasCorrectionFactor = maxBiasInForegroundMask;
                    biasIntensityScaleFactor->SetPixel(currOutIndex,
                                                       static_cast<InternalImagePixelType>(multiplicitiveBiasCorrectionFactor));
                  }
                  else if (multiplicitiveBiasCorrectionFactor < minBiasInForegroundMask) //
                                                                                         //
                                                                                         // CLAMP
                  {
                    multiplicitiveBiasCorrectionFactor = minBiasInForegroundMask;
                    biasIntensityScaleFactor->SetPixel(currOutIndex,
                                                       static_cast<InternalImagePixelType>(multiplicitiveBiasCorrectionFactor));
                  }

                  typename MaskNNInterpolationType::OutputType allTissueMaskValue = 0;
                  if (allTissueMaskInterp->IsInsideBuffer(currOutPoint))
                  {
                    allTissueMaskValue = allTissueMaskInterp->Evaluate(currOutPoint);
                  }

                  if (allTissueMaskValue == 0)
                  {
                    // Now clamp intensities outside the probability mask region to
                    // the min and
                    // max of inside the mask region.
                    curOutput->SetPixel(currOutIndex, 0);
                  }
                  else
                  {
                    const double originalPixelValue = (*imIt)->GetPixel(currOutIndex);
                    const double correctedPixelValue = originalPixelValue * multiplicitiveBiasCorrectionFactor;
                    curOutput->SetPixel(currOutIndex, (InputImagePixelType)correctedPixelValue);
                  }
                } // for currOutIndex[0]
              }
            }
          });

        std::cout << "Standardizing Bias Corrected Intensities: ..." << std::endl;
        curOutput = StandardizeMaskIntensity<InputImageType, ByteImageType>(curOutput,
                                                                            this->m_ForegroundBrainMask,
                                                                            0.0005,
                                                                            1.0 - 0.0005,
                                                                            1,
                                                                            0.95 * MAX_IMAGE_OUTPUT_VALUE,
                                                                            0,
                                                                            MAX_IMAGE_OUTPUT_VALUE);
        std::cout << "done." << std::endl;

        outputs[mapIt->first].push_back(curOutput);

        if (this->m_DebugLevel > 7)
        { // DEBUG:  This code is for debugging purposes only;
          using WriterType = itk::ImageFileWriter<InputImageType>;
          typename WriterType::Pointer writer = WriterType::New();
          writer->UseCompressionOn();

          std::stringstream CurrentIterationID_stream("");
          CurrentIterationID_stream << CurrentIterationID;
          std::stringstream template_index_stream("");
          template_index_stream << mapIt->first << std::distance(mapIt->second.begin(), imIt);
          const std::string fn = this->m_OutputDebugDir + "/BIAS_INDEX_" + template_index_stream.str() + "_LEVEL_" +
                                 CurrentIterationID_stream.str() + ".nii.gz";
          writer->SetInput(biasIntensityScaleFactor);
          writer->SetFileName(fn.c_str());
          writer->Update();
          muLogMacro(<< "DEBUG:  Wrote image " << fn << std::endl);
        }
        //
      } // for ichan
    }
  }
  CorrectImagesTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime = CorrectImagesTimer.GetTotal();
  muLogMacro(<< "Correcting Images took " << elapsedTime << " " << CorrectImagesTimer.GetUnit() << std::endl);
  return outputs;
}

#endif
