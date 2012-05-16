#ifndef __LLSBiasCorrector_hxx
#define __LLSBiasCorrector_hxx

#include <float.h>
#include <math.h>
#include <cmath>
#include <iostream>

#include "Log.h"
#include "StandardizeMaskIntensity.h"
#include "LLSBiasCorrector.h"
#include "vnl/vnl_math.h"

#include "ComputeDistributions.h"

#define USE_HALF_RESOLUTION 1
#define MIN_SKIP_SIZE 2

// Use the normal equation? Less accurate, but requires much less memory
#define LLSBIAS_USE_NORMAL_EQUATION 1
//
//
// //////////////////////////////////////////////////////////////////////////////
#define mypow(a, b) vcl_pow( ( a ), static_cast<double>(b) )
//
//
// //////////////////////////////////////////////////////////////////////////////

template <class TInputImage, class TProbabilityImage>
LLSBiasCorrector<TInputImage, TProbabilityImage>
::LLSBiasCorrector()
{
  m_InputImages.clear();

  m_MaxDegree = 4;

  m_DebugLevel = 0;
  m_OutputDebugDir = "";
  m_SampleSpacing = 4.0;
  m_WorkingSpacing = 1.0;

  // m_MaximumBiasMagnitude = .1;

  m_XMu[0] = 0.0;
  m_XMu[1] = 0.0;
  m_XMu[2] = 0.0;

  m_XStd[0] = 1.0;
  m_XStd[1] = 1.0;
  m_XStd[2] = 1.0;
}

template <class TInputImage, class TProbabilityImage>
LLSBiasCorrector<TInputImage, TProbabilityImage>
::~LLSBiasCorrector()
{
  m_ForegroundBrainMask = 0;
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::CheckInputs()
{
  // if (m_MaxDegree == 0)
  // itkExceptionMacro(<< "Max bias degree is zero" << std::endl );

  if( m_InputImages.size() == 0 )
    {
    itkExceptionMacro(<< "No input image specified" << std::endl );
    }

  if( m_InputImages[0]->GetImageDimension() != 3 )
    {
    itkExceptionMacro(<< "Input dimension invalid: only supports 3D images" << std::endl );
    }

  if( m_BiasPosteriors.size() < 1 )
    {
    itkExceptionMacro(<< "Must have one or more class probabilities" << std::endl );
    }

  const InputImageSizeType size = m_InputImages[0]->GetLargestPossibleRegion().GetSize();
  for( unsigned int i = 1; i < m_InputImages.size(); i++ )
    {
    const InputImageSizeType size_i = m_InputImages[i]->GetLargestPossibleRegion().GetSize();
    if( size != size_i )
      {
      itkExceptionMacro(<< "Image sizes do not match" << std::endl );
      }
    }
  for( unsigned int i = 0; i < m_BiasPosteriors.size(); i++ )
    {
    if( m_BiasPosteriors[i]->GetImageDimension() != 3 )
      {
      itkExceptionMacro(<< "Probability [" << i << "] has invalid dimension: only supports 3D images" << std::endl );
      }
    const ProbabilityImageSizeType psize = m_BiasPosteriors[i]->GetLargestPossibleRegion().GetSize();
    if( size[0] != psize[0] || size[1] != psize[1] || size[2] != psize[2] )
      {
      itkExceptionMacro(<< "Image data and probabilities 3D size mismatch" << std::endl );
      }
    }
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::ComputeDistributions()
{
  muLogMacro(<< "LLSBiasCorrector: Computing means and variances..." << std::endl );

  std::vector<double> ClassProbabilityWeightings;
  CombinedComputeDistributions<TInputImage, TProbabilityImage, MatrixType>(this->m_CandidateRegions, m_InputImages,
                                                                           m_BiasPosteriors,
                                                                           this->m_ListOfClassStatistics, //
                                                                                                          //
                                                                                                          // This
                                                                                                          //
                                                                                                          // is
                                                                                                          //
                                                                                                          // an
                                                                                                          //
                                                                                                          // output!
                                                                           this->m_DebugLevel,
                                                                           true);
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::SetListOfClassStatistics(const std::vector<RegionStats> & regionStats)
{
  this->m_ListOfClassStatistics = regionStats;
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::SetMaxDegree(unsigned int n)
{
  muLogMacro(<< "SetMaxDegree" << std::endl );

  m_MaxDegree = n;

  // NOTE:  update basis equations, Need to call the Initialize in order to
  // update
  // the basis equations.
  if( !m_ForegroundBrainMask.IsNull() )
    {
    this->Initialize();
    }
#if 1 // HACK
  if( m_BiasPosteriors.size() > 0 )
    {
    this->SetProbabilities(m_BiasPosteriors, this->m_CandidateRegions);
    }
#endif
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::SetSampleSpacing(double s)
{
  muLogMacro(<< "SetSampleSpacing" << std::endl );

  m_SampleSpacing = s;

#if 1 // HACK
  if( m_BiasPosteriors.size() > 0 )
    {
    this->SetProbabilities(m_BiasPosteriors, this->m_CandidateRegions);
    }
#endif
  // NOTE:  update basis equations, Need to call the Initialize in order to
  // update
  // the basis equations.
  if( m_ForegroundBrainMask.IsNotNull() )
    {
    this->Initialize();
    }
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::SetForegroundBrainMask(ByteImageType *mask)
{
  m_ForegroundBrainMask = mask;
  this->Initialize();
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::Initialize()
{
  const InputImageSizeType size
    = m_ForegroundBrainMask->GetLargestPossibleRegion().GetSize();

  // Compute skips along each dimension
#ifdef USE_HALF_RESOLUTION // Need to do full scale to prevent checkerboard
  // images from being created
  const InputImageSpacingType spacing = m_ForegroundBrainMask->GetSpacing();
  unsigned int                skips[3];
  skips[0] = (unsigned int)( m_SampleSpacing / spacing[0] );
  skips[1] = (unsigned int)( m_SampleSpacing / spacing[1] );
  skips[2] = (unsigned int)( m_SampleSpacing / spacing[2] );

  if( skips[0] < MIN_SKIP_SIZE )
    {
    skips[0] = MIN_SKIP_SIZE;
    }
  if( skips[1] < MIN_SKIP_SIZE )
    {
    skips[1] = MIN_SKIP_SIZE;
    }
  if( skips[2] < MIN_SKIP_SIZE )
    {
    skips[2] = MIN_SKIP_SIZE;
    }
  muLogMacro(<< "Sample skips: " << skips[0] << " x " << skips[1] << " x " << skips[2] << std::endl);
#else
  const unsigned int skips[3] = {1, 1, 1};
#endif

  const unsigned int numCoefficients
    = ( m_MaxDegree + 1 ) * ( m_MaxDegree + 2 ) / 2 * ( m_MaxDegree + 3 ) / 3;

  // Number of pixels with non-zero weights, downsampled
  unsigned numEquations = 0;
  // Assume that only 0.125 of the image is part of mask
  m_ValidIndicies.reserve(size[2] * size[1] * size[0] / (8 * skips[0] * skips[1] * skips[2]) );

  m_ValidIndicies.resize(0);
    {
    // Not parallizable! ORDER IS IMPORTANT  //#pragma omp parallel for
    // reduction(+:numEquations) default(shared)
    for( long kk = 0; kk < (long)size[2]; kk += skips[2] )
      {
      for( long jj = 0; jj < (long)size[1]; jj += skips[1] )
        {
        for( long ii = 0; ii < (long)size[0]; ii += skips[0] )
          {
          const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
          if( m_ForegroundBrainMask->GetPixel(currIndex) != 0 )
            {
            m_ValidIndicies.push_back(currIndex);
            // numEquations++;
            }
          }
        }
      }
    }
  numEquations = m_ValidIndicies.size();
  muLogMacro(<< "Linear system size = " << numEquations << " x " << numCoefficients << std::endl);

  // Make sure that number of equations >= number of unknowns
  if( numEquations < numCoefficients )
    {
    std::cout << size << std::endl;
    itkExceptionMacro(<< "Number of unknowns exceed number of equations:" << numEquations << " < " << numCoefficients);
    }

  // Create basis matrix

  muLogMacro(<< "Computing polynomial basis functions..." << std::endl );

  m_Basis.set_size(numEquations, numCoefficients);

    {
    // Coordinate scaling and offset parameters
    unsigned long long int local_XMu_x = 0;
    unsigned long long int local_XMu_y = 0;
    unsigned long long int local_XMu_z = 0;
      {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared) reduction(+:local_XMu_x,local_XMu_y,local_XMu_z)
#endif
      for( unsigned int kk = 0; kk < numEquations; kk++ )
        {
        const ProbabilityImageIndexType & currIndex = m_ValidIndicies[kk];
        local_XMu_x += currIndex[0];
        local_XMu_y += currIndex[1];
        local_XMu_z += currIndex[2];
        }
      }
    const double invNumEquations = 1.0 / static_cast<double>(numEquations);
    m_XMu[0] = static_cast<double>(local_XMu_x) * invNumEquations;
    m_XMu[1] = static_cast<double>(local_XMu_y) * invNumEquations;
    m_XMu[2] = static_cast<double>(local_XMu_z) * invNumEquations;
    }

    {
    double local_XStd_x = 0.0;
    double local_XStd_y = 0.0;
    double local_XStd_z = 0.0;
      {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared) reduction(+:local_XStd_x,local_XStd_y,local_XStd_z)
#endif
      for( unsigned int kk = 0; kk < numEquations; kk++ )
        {
        const ProbabilityImageIndexType & currIndex = m_ValidIndicies[kk];
        const double                      diff0 = static_cast<double>(currIndex[0]) - m_XMu[0];
        local_XStd_x += diff0 * diff0;
        const double diff1 = static_cast<double>(currIndex[1]) - m_XMu[1];
        local_XStd_y += diff1 * diff1;
        const double diff2 = static_cast<double>(currIndex[2]) - m_XMu[2];
        local_XStd_z += diff2 * diff2;
        }
      }
    m_XStd[0] = vcl_sqrt(local_XStd_x / numEquations);
    m_XStd[1] = vcl_sqrt(local_XStd_y / numEquations);
    m_XStd[2] = vcl_sqrt(local_XStd_z / numEquations);
    }

    {
    // Row and column indices
    // Fill in polynomial basis values
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
    for( unsigned int r = 0; r < numEquations; r++ )
      {
      const ProbabilityImageIndexType & currIndex = m_ValidIndicies[r];
      unsigned int                      c = 0;
      for( unsigned int order = 0; order <= m_MaxDegree; order++ )
        {
        for( unsigned int xorder = 0; xorder <= order; xorder++ )
          {
          for( unsigned int yorder = 0; yorder <= ( order - xorder ); yorder++ )
            {
            const int zorder = order - xorder - yorder;

            const double xc = ( currIndex[0] - m_XMu[0] ) / m_XStd[0];
            const double yc = ( currIndex[1] - m_XMu[1] ) / m_XStd[1];
            const double zc = ( currIndex[2] - m_XMu[2] ) / m_XStd[2];

            m_Basis(r, c)
              = mypow(xc, xorder) * mypow(yc, yorder) * mypow(zc, zorder);
            c++;
            }
          }
        }
      }
    }
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector<TInputImage, TProbabilityImage>
::SetProbabilities(const std::vector<ProbabilityImagePointer> & probs,
                   const std::vector<typename ByteImageType::Pointer> & candidateRegions )
{
  muLogMacro(<< "SetProbabilities" << std::endl );

  if( probs.size() < 1 )
    {
    itkExceptionMacro(<< "Need one or more probabilities" << std::endl );
    }
  if( probs.size() != candidateRegions.size() )
    {
    itkExceptionMacro(<< "Probabilities and CandidateRegions size must be the same." << std::endl );
    }
  for( unsigned int i = 0; i < probs.size(); i++ )
    {
    if( probs[i].IsNull() )
      {
      itkExceptionMacro(<< "One of input probabilities not initialized" << std::endl );
      }
    }
  for( unsigned int i = 0; i < candidateRegions.size(); i++ )
    {
    if( candidateRegions[i].IsNull() )
      {
      itkExceptionMacro(<< "One of input candidate regions is not initialized" << std::endl );
      }
    }
  m_BiasPosteriors = probs;
  m_CandidateRegions = candidateRegions;
}

template <class TInputImage, class TProbabilityImage>
std::vector<typename TInputImage::Pointer>
LLSBiasCorrector<TInputImage, TProbabilityImage>
::CorrectImages(const unsigned int CurrentIterationID)
{
  muLogMacro(<< "Correct Images" << std::endl );
  itk::TimeProbe CorrectImagesTimer;
  CorrectImagesTimer.Start();
  // Verify input
  this->CheckInputs();

  const InputImageSizeType size = m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  // Compute means and variances
  this->ComputeDistributions();

#ifdef USE_HALF_RESOLUTION
  // Compute skips along each dimension
  const InputImageSpacingType spacing = m_InputImages[0]->GetSpacing();

  unsigned int sampleofft[3];
  sampleofft[0] = (unsigned int)vcl_floor(m_SampleSpacing / spacing[0]);
  sampleofft[1] = (unsigned int)vcl_floor(m_SampleSpacing / spacing[1]);
  sampleofft[2] = (unsigned int)vcl_floor(m_SampleSpacing / spacing[2]);

  if( sampleofft[0] < MIN_SKIP_SIZE )
    {
    sampleofft[0] = MIN_SKIP_SIZE;
    }
  if( sampleofft[1] < MIN_SKIP_SIZE )
    {
    sampleofft[1] = MIN_SKIP_SIZE;
    }
  if( sampleofft[2] < MIN_SKIP_SIZE )
    {
    sampleofft[2] = MIN_SKIP_SIZE;
    }

  muLogMacro(
    << "Sample offsets: " << sampleofft[0] << " x " << sampleofft[1] << " x " << sampleofft[2] << std::endl);
#else
  const unsigned int sampleofft[3] = {1, 1, 1};
#endif

#ifdef USE_HALF_RESOLUTION  // Need more accurate value the downsampling was
  // causing images with zeros to be produced, and the
  // bspline registrations were not doing very much
  // because of this.
  unsigned int workingofft[3];
  workingofft[0] = (unsigned int)vcl_floor(m_WorkingSpacing / spacing[0]);
  workingofft[1] = (unsigned int)vcl_floor(m_WorkingSpacing / spacing[1]);
  workingofft[2] = (unsigned int)vcl_floor(m_WorkingSpacing / spacing[2]);

  if( workingofft[0] < MIN_SKIP_SIZE )
    {
    workingofft[0] = MIN_SKIP_SIZE;
    }
  if( workingofft[1] < MIN_SKIP_SIZE )
    {
    workingofft[1] = MIN_SKIP_SIZE;
    }
  if( workingofft[2] < MIN_SKIP_SIZE )
    {
    workingofft[2] = MIN_SKIP_SIZE;
    }
  muLogMacro(
    << "Working offsets: " << workingofft[0] << " x " << workingofft[1] << " x " << workingofft[2] << std::endl);
#else
  //  const unsigned int workingofft[3] ={ {1,1,1} };
#endif

  const unsigned int numChannels = m_InputImages.size();
  const unsigned int numClasses = m_BiasPosteriors.size();

  /* if m_MaxDegree = 4/3/2, then this is 35/20/10 */
  const unsigned int numCoefficients
    = ( m_MaxDegree + 1 ) * ( m_MaxDegree + 2 ) / 2 * ( m_MaxDegree + 3 ) / 3;

  muLogMacro(<< numClasses << " classes\n" << std::endl );
  muLogMacro(<< numCoefficients << " coefficients\n" << std::endl );

  /* compute inverse matrix for each tissue type */
  muLogMacro(<< "Computing inverse covars...\n" << std::endl );
  std::vector<MatrixType> invCovars;
  for( unsigned int iclass = 0; iclass < numClasses; iclass++ )
    {
    invCovars.push_back( MatrixInverseType(this->m_ListOfClassStatistics[iclass].m_Covariance) );
    }

  // Create matrices and vectors
  // lhs = replicated basis polynomials for each channel, weighted by inv cov
  // rhs = difference image between original and reconstructed mean image

  muLogMacro(<< "Creating matrices for LLS..." << std::endl );

  const unsigned int numEquations = m_Basis.rows();   /*  since m_Basis is [numEqu *
                                                          35 ] matrix */

  muLogMacro(
    << numEquations << " equations, " << numCoefficients << " coefficients" << std::endl );

  // Compute  orthogonal transpose component of basis
  muLogMacro(<< "Computing ortho part of basis" << std::endl );
#if 1
  // Note: vnl_qr gives Q mxm and R mxn for A mxn

  MatrixType basisT;
    {
    MatrixQRType qr(m_Basis);

    // Get economy size R (square)
    MatrixType R(numCoefficients, numCoefficients, 0);

      {
      MatrixType Rfull = qr.R();  /* right triangular matrix */
      for( unsigned int r = 0; r < numCoefficients; r++ )
        {
        for( unsigned int c = r; c < numCoefficients; c++ )
          {
          R(r, c) = Rfull(r, c);
          }
        }
      // Rfull.set_size(1, 1);
      }
    // NOTE: to get mxn Q from vnl_qr, Q'*Q = id nxn
    basisT = m_Basis * MatrixInverseType(R);
    }
  basisT.inplace_transpose(); // basisT = Q'
#else
  // Do this instead for ordinary weighted LSQ
  MatrixType basisT = m_Basis.transpose();
#endif

#if LLSBIAS_USE_NORMAL_EQUATION
  MatrixType lhs(numCoefficients * numChannels, numCoefficients * numChannels);

  MatrixType rhs(numCoefficients * numChannels, 1);

#else
  MatrixType lhs(numEquations * numChannels, numCoefficients * numChannels);

  MatrixType rhs(numEquations * numChannels, 1);

#endif

  muLogMacro(<< "Fill rhs" << std::endl );
  rhs.fill(0.0);

    {
    // Compute ratio between original and flat image, weighted using posterior
    // probability and inverse covariance
    for( LOOPITERTYPE ichan = 0; ichan < (LOOPITERTYPE) numChannels; ichan++ )
      {
      MatrixType R_i(numEquations, 1, 0.0);
      for( unsigned int jchan = 0; jchan < numChannels; jchan++ )
        {
          {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
          for( unsigned int eq = 0; eq < numEquations; eq++ )
            {
            const ProbabilityImageIndexType & currIndex = m_ValidIndicies[eq];
            // Compute reconstructed intensity, weighted by prob * invCov
            double sumW = DBL_EPSILON;
            double recon = 0;
            for( unsigned int iclass = 0; iclass < numClasses; iclass++ )
              {
              const MatrixType & invCov = invCovars[iclass];
              const double       w = m_BiasPosteriors[iclass]->GetPixel(currIndex)
                * invCov(ichan, jchan);
              sumW += w;
              recon += w * this->m_ListOfClassStatistics[iclass].m_Means[jchan];
              }
            recon /= sumW;

            const double bias = LOGP( m_InputImages[jchan]->GetPixel(currIndex) ) - recon;
            R_i(eq, 0) += sumW * bias;
            }
          }
        } // for jchan

#if LLSBIAS_USE_NORMAL_EQUATION
      R_i = basisT * R_i;
      for( unsigned int row = 0; row < numCoefficients; row++ )
        {
        rhs(ichan * numCoefficients + row, 0) = R_i(row, 0);
        }
#else
      for( unsigned int row = 0; row < numEquations; row++ )
        {
        rhs(ichan * numEquations + row, 0) = R_i(row, 0);
        }
#endif
      } // for ichan
    }

  muLogMacro(<< "Fill lhs" << std::endl );

  // Compute LHS using replicated basis entries, weighted using posterior
  // probability and inverse covariance
    {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
    for( LOOPITERTYPE ichan = 0; ichan < (LOOPITERTYPE)numChannels; ichan++ )
      {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
      for( LOOPITERTYPE jchan = 0; jchan < (LOOPITERTYPE)numChannels; jchan++ )
        {
        MatrixType Wij_A(numEquations, numCoefficients, 0.0);
          {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for shared(Wij_A) // default(shared)
#endif
          for( unsigned int eq = 0; eq < numEquations; eq++ )
            {
            const ProbabilityImageIndexType & currIndex = m_ValidIndicies[eq];
            double                            sumW = DBL_EPSILON;
            for( unsigned int iclass = 0; iclass < numClasses; iclass++ )
              {
              const MatrixType & invCov = invCovars[iclass];
              double             w = m_BiasPosteriors[iclass]->GetPixel(currIndex)
                * invCov(ichan, jchan);
              sumW += w;
              }
            for( unsigned int col = 0; col < numCoefficients; col++ )
              {
              Wij_A(eq, col) = sumW * m_Basis(eq, col);
              }
            }
          }

#if LLSBIAS_USE_NORMAL_EQUATION
        MatrixType lhs_ij = basisT * Wij_A;
        for( unsigned int row = 0; row < numCoefficients; row++ )
          {
          for( unsigned int col = 0; col < numCoefficients; col++ )
            {
            lhs(row + ichan * numCoefficients, col + jchan * numCoefficients)
              = lhs_ij(row, col);
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
    }

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
  if( !std::isfinite( (double)coeffs[0][0]) )
#else
  if( coeffs[0][0] != std::numeric_limits::infinity() )
#endif
    {
    itkExceptionMacro(
      << "\ncoeffs: \n" << coeffs
      // << "\nlhs_ij: \n" << lhs_ij
      << "\nbasisT: \n" << basisT
      // << "\nWij_A: \n" << Wij_A
      << "\nlhs: \n" << lhs
      << "\nrhs: \n" << rhs
      );
    }
  // Clear memory for the basis transpose
  basisT.set_size(0, 0);

  if( this->m_DebugLevel > 9 )
    {
    muLogMacro(<< "Bias field coeffs after LLS:" << std::endl  << coeffs);
    }

  // Remove bias
  muLogMacro(<< "Correcting input images..." << std::endl );

  std::vector<InputImagePointer> outputs(numChannels);
  for( unsigned int ichan = 0; ichan < numChannels; ichan++ )
    {
    outputs[ichan] = InternalImageType::New();
    outputs[ichan]->CopyInformation(this->m_InputImages[ichan]);
    outputs[ichan]->SetRegions( this->m_InputImages[ichan]->GetLargestPossibleRegion() );
    outputs[ichan]->Allocate();
    outputs[ichan]->FillBuffer(0);

    // Compute the vcl_log transformed bias field
    InternalImagePointer biasIntensityScaleFactor = InternalImageType::New();
    biasIntensityScaleFactor->CopyInformation(this->m_InputImages[ichan]);
    biasIntensityScaleFactor->SetRegions( this->m_InputImages[ichan]->GetLargestPossibleRegion() );
    biasIntensityScaleFactor->Allocate();
    biasIntensityScaleFactor->FillBuffer(0);

    double maxBiasInForegroundMask = vcl_numeric_limits<double>::min();
    double minBiasInForegroundMask = vcl_numeric_limits<double>::max();

      {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for shared(maxBiasInForegroundMask,minBiasInForegroundMask) default(shared)
#endif
      for( long kk = 0; kk < (long)size[2]; kk++ )
        {
        for( long jj = 0; jj < (long)size[1]; jj++ )
          {
          for( long ii = 0; ii < (long)size[0]; ii++ )
            {
            const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
            double                          logFitValue = 0.0;
            unsigned int                    c = ichan * numCoefficients;
            for( unsigned int order = 0; order <= m_MaxDegree; order++ )
              {
              for( unsigned int xorder = 0; xorder <= order; xorder++ )
                {
                for( unsigned int yorder = 0; yorder <= ( order - xorder ); yorder++ )
                  {
                  const int zorder = order - xorder - yorder;

                  const double xc = ( currIndex[0] - m_XMu[0] ) / m_XStd[0];
                  const double yc = ( currIndex[1] - m_XMu[1] ) / m_XStd[1];
                  const double zc = ( currIndex[2] - m_XMu[2] ) / m_XStd[2];

                  const double poly
                    = mypow(xc, xorder) * mypow(yc, yorder) * mypow(zc, zorder);

                  // logFitValue += coeffs[c] * poly;
                  logFitValue += coeffs(c, 0) * poly;
                  c++;
                  }
                }
              }

            const ByteImagePixelType maskValue = m_ForegroundBrainMask->GetPixel(currIndex);
            /* NOTE:  For regions listed as background, clamp the outputs[ichan]
              */
            if( vnl_math_isnan(logFitValue) || vnl_math_isinf(logFitValue) )
              {
              std::cout << "WARNING:  Bad Scale Value Computed!" << std::endl;
              logFitValue = 0.0;
              }
            const double multiplicitiveBiasCorrectionFactor = 1.0 / EXPP(logFitValue);
            if( maskValue != 0 )
              {
              if( multiplicitiveBiasCorrectionFactor > maxBiasInForegroundMask )
                {
                maxBiasInForegroundMask = multiplicitiveBiasCorrectionFactor;
                }
              if( multiplicitiveBiasCorrectionFactor < minBiasInForegroundMask )
                {
                minBiasInForegroundMask = multiplicitiveBiasCorrectionFactor;
                }
              }
            biasIntensityScaleFactor->SetPixel(currIndex, (InternalImagePixelType)multiplicitiveBiasCorrectionFactor);
            } // for currIndex[0]
          }
        }
      }
    std::cout << "Foreground Mask Bias Correction MIN: " << minBiasInForegroundMask << " MAX: "
              << maxBiasInForegroundMask << std::endl;
    // Correct image using (clamped) bias field)
      {
#if defined(LOCAL_USE_OPEN_MP)
#pragma omp parallel for default(shared)
#endif
      for( long kk = 0; kk < (long)size[2]; kk++ )
        {
        for( long jj = 0; jj < (long)size[1]; jj++ )
          {
          for( long ii = 0; ii < (long)size[0]; ii++ )
            {
            const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
            double                          multiplicitiveBiasCorrectionFactor = biasIntensityScaleFactor->GetPixel(
                currIndex);
            if( multiplicitiveBiasCorrectionFactor > maxBiasInForegroundMask )  //
                                                                                //
                                                                                // CLAMP
              {
              multiplicitiveBiasCorrectionFactor = maxBiasInForegroundMask;
              biasIntensityScaleFactor->SetPixel(currIndex,
                                                 (InternalImagePixelType) multiplicitiveBiasCorrectionFactor );
              }
            else if( multiplicitiveBiasCorrectionFactor < minBiasInForegroundMask )  //
                                                                                     //
                                                                                     // CLAMP
              {
              multiplicitiveBiasCorrectionFactor = minBiasInForegroundMask;
              biasIntensityScaleFactor->SetPixel(currIndex,
                                                 (InternalImagePixelType) multiplicitiveBiasCorrectionFactor );
              }
            if( this->m_AllTissueMask->GetPixel(currIndex) == 0 )
              {
              // Now clamp intensities outside the probability mask region to
              // the min and
              // max of inside the mask region.
              outputs[ichan]->SetPixel(currIndex, 0 );
              }
            else
              {
              const double originalPixelValue = m_InputImages[ichan]->GetPixel(currIndex);
              const double correctedPixelValue = originalPixelValue * multiplicitiveBiasCorrectionFactor;
              outputs[ichan]->SetPixel(currIndex, (InputImagePixelType)correctedPixelValue);
              }
            } // for currIndex[0]
          }
        }
      }

      {
      std::cout << "Standardizing Bias Corrected Intensities: ...";
      outputs[ichan]
        = StandardizeMaskIntensity<InputImageType, ByteImageType>(
            outputs[ichan],
            this->m_ForegroundBrainMask,
            0.0005, 1.0 - 0.0005,
            1, 0.95 * MAX_IMAGE_OUTPUT_VALUE,
            0, MAX_IMAGE_OUTPUT_VALUE);
      std::cout << "done." << std::endl;
      }
    if( this->m_DebugLevel > 7 )
      { // DEBUG:  This code is for debugging purposes only;
      typedef itk::ImageFileWriter<InputImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();

      std::stringstream CurrentIterationID_stream("");
      CurrentIterationID_stream << CurrentIterationID;
      std::stringstream template_index_stream("");
      template_index_stream << ichan;
      const std::string fn = this->m_OutputDebugDir + "/BIAS_INDEX_" + template_index_stream.str() + "_LEVEL_"
        + CurrentIterationID_stream.str() + ".nii.gz";
      writer->SetInput(biasIntensityScaleFactor);
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
      }
    //
    } // for ichan
  // Remove internal references to input images when done
  m_InputImages.clear();
  CorrectImagesTimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime = CorrectImagesTimer.GetTotal();
  muLogMacro(<< "Correcting Images took " << elapsedTime << " " << CorrectImagesTimer.GetUnit() << std::endl);
  return outputs;
}

#endif
