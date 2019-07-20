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
//
//
// //////////////////////////////////////////////////////////////////////////////
//
// Bias field estimation using linear least squares polynomial fitting.
// Requires input image intensities and probability estimates for each voxel.
//
// Corrections are done either on subsampled grid or full resolution.
// This is designed to be run iteratively. Corrections should not be
// accumulated. Always use original image as the input, otherwise may get
// strange results.
//
// Van Leemput K, Maes F, Vandermeulen D, Suetens P. Automated model based
// bias field correction of MR images of the brain. IEEE TMI 1999; 18:885-896.
//
//
//
// //////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2004

#ifndef __LLSBiasCorrector_h
#define __LLSBiasCorrector_h

#include "itkImage.h"
#include "itkObject.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_qr.h"
#include "vnl/algo/vnl_svd.h"

#include <vector>
#include <list>
#include <map>
/** \class LLSBiasCorrector
 */
template < typename TInputImage, typename TProbabilityImage >
class LLSBiasCorrector : public itk::Object
{
public:
  /** Standard class type alias. */
  using Self = LLSBiasCorrector;
  using Pointer = itk::SmartPointer< Self >;
  using ConstPointer = itk::SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** The dimension of the image. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  // Image types
  using InputImageType = TInputImage;
  using InputImagePointer = typename TInputImage::Pointer;
  using InputImageIndexType = typename TInputImage::IndexType;
  using InputImagePixelType = typename TInputImage::PixelType;
  using InputImageRegionType = typename TInputImage::RegionType;
  using InputImageSizeType = typename TInputImage::SizeType;
  using InputImageSpacingType = typename TInputImage::SpacingType;

  using InputImageVector = std::vector< InputImagePointer >;
  using MapOfInputImageVectors = orderedmap< std::string, InputImageVector >;

  using ByteImageType = itk::Image< unsigned char, Self::ImageDimension >;
  using ByteImagePointer = typename ByteImageType::Pointer;
  using ByteImageIndexType = typename ByteImageType::IndexType;
  using ByteImageOffsetType = typename ByteImageType::OffsetType;
  using ByteImagePixelType = typename ByteImageType::PixelType;
  using ByteImageRegionType = typename ByteImageType::RegionType;
  using ByteImageSizeType = typename ByteImageType::SizeType;

  using ProbabilityImageType = TProbabilityImage;
  using ProbabilityImagePointer = typename ProbabilityImageType::Pointer;
  using ProbabilityImageIndexType = typename ProbabilityImageType::IndexType;
  using ProbabilityImagePixelType = typename ProbabilityImageType::PixelType;
  using ProbabilityImageRegionType = typename ProbabilityImageType::RegionType;
  using ProbabilityImageSizeType = typename ProbabilityImageType::SizeType;

  using InternalImageType = itk::Image< float, 3 >;
  using InternalImagePointer = InternalImageType::Pointer;
  using InternalImageIndexType = InternalImageType::IndexType;
  using InternalImagePixelType = InternalImageType::PixelType;
  using InternalImageRegionType = InternalImageType::RegionType;
  using InternalImageSizeType = InternalImageType::SizeType;

  using InputImageNNInterpolationType = itk::NearestNeighborInterpolateImageFunction< InputImageType, double >;
  using MaskNNInterpolationType = itk::NearestNeighborInterpolateImageFunction< ByteImageType, double >;

  using ScalarType = double;

  using MatrixType = vnl_matrix< ScalarType >;
  using VectorType = vnl_vector< ScalarType >;

  using MatrixInverseType = vnl_matrix_inverse< ScalarType >;
  using MatrixQRType = vnl_qr< ScalarType >;
  using MatrixSVDType = vnl_svd< ScalarType >;

  // The maximum polynomial degree of the bias field estimate
  void
  SetMaxDegree( unsigned int );

  itkGetMacro( MaxDegree, unsigned int );

  // Spacing for determining voxels in LLS
  void
  SetSampleSpacing( double s );

  itkGetMacro( SampleSpacing, double );

  // Spacing for determining which voxels need to be updated
  // if correction is not done at full resolution
  itkSetMacro( WorkingSpacing, double );
  itkGetMacro( WorkingSpacing, double );

  // Bias field max magnitude
  // itkSetMacro(MaximumBiasMagnitude, double);
  // itkGetMacro(MaximumBiasMagnitude, double);

  void
  Initialize();

  itkSetObjectMacro( AllTissueMask, ByteImageType );
  itkGetConstObjectMacro( AllTissueMask, ByteImageType );

  void
  SetForegroundBrainMask( ByteImageType * mask );

  void
  SetInputImages( MapOfInputImageVectors inputs )
  {
    this->m_InputImages = inputs;
    this->Modified();
  }

  void
  SetProbabilities( const std::vector< ProbabilityImagePointer > &         probs,
                    const std::vector< typename ByteImageType::Pointer > & candidateRegions );

  void
  SetListOfClassStatistics( const std::vector< RegionStats > & regionStats );

  // Set/Get the Debugging level for filter verboseness
  itkSetMacro( DebugLevel, unsigned int );
  itkGetMacro( DebugLevel, unsigned int );

  // Set/Get the Debugging level for filter verboseness
  itkSetMacro( OutputDebugDir, std::string );
  itkGetMacro( OutputDebugDir, std::string );
  // Correct input images and write it to the designated output
  // fullRes flag selects whether to correct whole image or just grid points
  // defined by WorkingSpacing
  MapOfInputImageVectors
  CorrectImages( const unsigned int CurrentIterationID );

protected:
  LLSBiasCorrector();
  ~LLSBiasCorrector();

  void
  CheckInputs();

  void
  ComputeDistributions();

private:
  InputImagePointer
  GetFirstInputImage()
  {
    return GetMapVectorFirstElement( this->m_InputImages );
  }
  MapOfInputImageVectors                   m_InputImages;
  std::vector< ProbabilityImageIndexType > m_ValidIndicies;
  ByteImagePointer                         m_ForegroundBrainMask;
  ByteImagePointer                         m_AllTissueMask;

  std::vector< ProbabilityImagePointer >         m_BiasPosteriors;
  std::vector< typename ByteImageType::Pointer > m_CandidateRegions;

  unsigned int m_DebugLevel;
  std::string  m_OutputDebugDir;

  unsigned int m_MaxDegree;

  double m_SampleSpacing;
  double m_WorkingSpacing;

  // double m_MaximumBiasMagnitude;

  std::vector< RegionStats > m_ListOfClassStatistics;
  MatrixType                 m_Basis;

  // Coordinate scaling and offset, computed from input probabilities
  // for preconditioning the polynomial basis equations
  double m_XMu[3];
  double m_XStd[3];
};

#ifndef MU_MANUAL_INSTANTIATION
#  include "LLSBiasCorrector.hxx"
#endif

#endif
