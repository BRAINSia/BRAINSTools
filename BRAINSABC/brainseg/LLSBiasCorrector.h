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
/** \class LLSBiasCorrector
 */
template <class TInputImage, class TProbabilityImage>
class LLSBiasCorrector : public itk::Object
{
public:

  /** Standard class typedefs. */
  typedef LLSBiasCorrector              Self;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** The dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  // Image types
  typedef TInputImage                       InputImageType;
  typedef typename TInputImage::Pointer     InputImagePointer;
  typedef typename TInputImage::IndexType   InputImageIndexType;
  typedef typename TInputImage::PixelType   InputImagePixelType;
  typedef typename TInputImage::RegionType  InputImageRegionType;
  typedef typename TInputImage::SizeType    InputImageSizeType;
  typedef typename TInputImage::SpacingType InputImageSpacingType;

  typedef itk::Image<unsigned char, itkGetStaticConstMacro(ImageDimension)> ByteImageType;
  typedef typename ByteImageType::Pointer                                   ByteImagePointer;
  typedef typename ByteImageType::IndexType                                 ByteImageIndexType;
  typedef typename ByteImageType::OffsetType                                ByteImageOffsetType;
  typedef typename ByteImageType::PixelType                                 ByteImagePixelType;
  typedef typename ByteImageType::RegionType                                ByteImageRegionType;
  typedef typename ByteImageType::SizeType                                  ByteImageSizeType;

  typedef TProbabilityImage                         ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer    ProbabilityImagePointer;
  typedef typename ProbabilityImageType::IndexType  ProbabilityImageIndexType;
  typedef typename ProbabilityImageType::PixelType  ProbabilityImagePixelType;
  typedef typename ProbabilityImageType::RegionType ProbabilityImageRegionType;
  typedef typename ProbabilityImageType::SizeType   ProbabilityImageSizeType;

  typedef itk::Image<float, 3>          InternalImageType;
  typedef InternalImageType::Pointer    InternalImagePointer;
  typedef InternalImageType::IndexType  InternalImageIndexType;
  typedef InternalImageType::PixelType  InternalImagePixelType;
  typedef InternalImageType::RegionType InternalImageRegionType;
  typedef InternalImageType::SizeType   InternalImageSizeType;

  typedef double ScalarType;

  typedef vnl_matrix<ScalarType> MatrixType;
  typedef vnl_vector<ScalarType> VectorType;

  typedef vnl_matrix_inverse<ScalarType> MatrixInverseType;
  typedef vnl_qr<ScalarType>             MatrixQRType;
  typedef vnl_svd<ScalarType>            MatrixSVDType;

  // The maximum polynomial degree of the bias field estimate
  void SetMaxDegree(unsigned int);

  itkGetMacro(MaxDegree, unsigned int);

  // Spacing for determining voxels in LLS
  void SetSampleSpacing(double s);

  itkGetMacro(SampleSpacing, double);

  // Spacing for determining which voxels need to be updated
  // if correction is not done at full resolution
  itkSetMacro(WorkingSpacing, double);
  itkGetMacro(WorkingSpacing, double);

  // Bias field max magnitude
  // itkSetMacro(MaximumBiasMagnitude, double);
  // itkGetMacro(MaximumBiasMagnitude, double);

  void Initialize();

  itkSetObjectMacro(AllTissueMask, ByteImageType);
  itkGetConstObjectMacro(AllTissueMask, ByteImageType);

  void SetForegroundBrainMask(ByteImageType *mask);

  void SetInputImages(const std::vector<InputImagePointer> & inputs)
  {
    this->m_InputImages = inputs;
    this->Modified();
  }

  void SetProbabilities(const std::vector<ProbabilityImagePointer> & probs,
                        const std::vector<typename ByteImageType::Pointer> & candidateRegions);

  void SetListOfClassStatistics(const std::vector<RegionStats> & regionStats);

  // Set/Get the Debugging level for filter verboseness
  itkSetMacro(DebugLevel, unsigned int);
  itkGetMacro(DebugLevel, unsigned int);

  // Set/Get the Debugging level for filter verboseness
  itkSetMacro(OutputDebugDir, std::string);
  itkGetMacro(OutputDebugDir, std::string);
  // Correct input images and write it to the designated output
  // fullRes flag selects whether to correct whole image or just grid points
  // defined by WorkingSpacing
  std::vector<InputImagePointer> CorrectImages(const unsigned int CurrentIterationID);

protected:

  LLSBiasCorrector();
  ~LLSBiasCorrector();

  void CheckInputs();

  void ComputeDistributions();

private:

  std::vector<InputImagePointer>         m_InputImages;
  std::vector<ProbabilityImageIndexType> m_ValidIndicies;
  ByteImagePointer                       m_ForegroundBrainMask;
  ByteImagePointer                       m_AllTissueMask;

  std::vector<ProbabilityImagePointer>         m_BiasPosteriors;
  std::vector<typename ByteImageType::Pointer> m_CandidateRegions;

  unsigned int m_DebugLevel;
  std::string  m_OutputDebugDir;

  unsigned int m_MaxDegree;

  double m_SampleSpacing;
  double m_WorkingSpacing;

  // double m_MaximumBiasMagnitude;

  std::vector<RegionStats> m_ListOfClassStatistics;
  MatrixType               m_Basis;

  // Coordinate scaling and offset, computed from input probabilities
  // for preconditioning the polynomial basis equations
  double m_XMu[3];
  double m_XStd[3];
};

#ifndef MU_MANUAL_INSTANTIATION
#include "LLSBiasCorrector.hxx"
#endif

#endif
