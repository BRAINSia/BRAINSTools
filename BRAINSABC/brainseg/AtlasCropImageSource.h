//
//
// //////////////////////////////////////////////////////////////////////////////
//
// Generate cropped images using ROI obtained from probabilities
//
//
//
// //////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 11/2003

#ifndef __AtlasCropImageSource_h
#define __AtlasCropImageSource_h

#include "itkObject.h"

#include <vector>

/** \class AtlasCropImageSource
 */
template <class TInputImage, class TProbabilityImage>
class AtlasCropImageSource : public itk::Object
{
public:

  /** Standard class typedefs. */
  typedef AtlasCropImageSource          Self;
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
  typedef typename TInputImage::OffsetType  InputImageOffsetType;
  typedef typename TInputImage::PixelType   InputImagePixelType;
  typedef typename TInputImage::PointType   InputImagePointType;
  typedef typename TInputImage::RegionType  InputImageRegionType;
  typedef typename TInputImage::SizeType    InputImageSizeType;
  typedef typename TInputImage::SpacingType InputImageSpacingType;

  typedef TProbabilityImage                          ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer     ProbabilityImagePointer;
  typedef typename ProbabilityImageType::IndexType   ProbabilityImageIndexType;
  typedef typename ProbabilityImageType::OffsetType  ProbabilityImageOffsetType;
  typedef typename ProbabilityImageType::PixelType   ProbabilityImagePixelType;
  typedef typename ProbabilityImageType::RegionType  ProbabilityImageRegionType;
  typedef typename ProbabilityImageType::SizeType    ProbabilityImageSizeType;
  typedef typename ProbabilityImageType::SpacingType ProbabilityImageSpacingType;

  typedef std::vector<ProbabilityImagePointer> ProbabilityImageList;

  typedef struct
    {
    InputImageIndexType offset;
    InputImageSizeType cropped_size;
    InputImageSizeType original_size;
    } CropInfoType;

  // Set/get output image padding, in mm
  itkGetMacro(Padding, double);
  itkSetMacro(Padding, double);

  bool CheckBounds();

  void UseProbabilities(ProbabilityImageList probs);

  // Create new images (either cropped or padded)
  InputImagePointer Restore(InputImagePointer img);

  InputImagePointer Crop(InputImagePointer img);

  // Crop region information
  void SetCropInfo(const CropInfoType & info)
  {
    m_CropInfo = info;
  }

  const CropInfoType & GetCropInfo()
  {
    return m_CropInfo;
  }

  // For debugging, generate slabs in last dim with top and bottom parts removed
  itkSetMacro(SlabMode, bool);
  itkGetConstMacro(SlabMode, bool);
  itkBooleanMacro(SlabMode);
protected:

  AtlasCropImageSource();
  ~AtlasCropImageSource()
  {
  }

  double m_Padding;

  InputImageIndexType m_LowerBound;
  InputImageIndexType m_UpperBound;

  InputImagePointType m_InputOrigin;
  InputImagePointType m_CropOrigin;

  InputImageSizeType m_OriginalSize;

  CropInfoType m_CropInfo;

  bool m_SlabMode;
};

#ifndef MU_MANUAL_INSTANTIATION
#include "AtlasCropImageSource.hxx"
#endif

#endif
