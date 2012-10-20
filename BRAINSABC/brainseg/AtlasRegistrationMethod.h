//
//
// //////////////////////////////////////////////////////////////////////////////
//
//  Registration of a dataset to an atlas using affine transformation and
//  MI image match metric
//
//  Only for 3D!
//
//  Given a list of filenames for atlas template and probabilities along with
//  the dataset, this class generate images that are in the space of the first
//  image (all data and probability images).
//
//
//
// //////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 10/2003

#ifndef __AtlasRegistrationMethod_h
#define __AtlasRegistrationMethod_h

#include "itkAffineTransform.h"
#include "itkArray.h"
#include "itkImage.h"
#include "itkObject.h"

#include <vector>

#include "itkAffineTransform.h"
#include "BRAINSFitHelper.h"

#include <string>

/** \class AtlasRegistrationMethod
 */
template <class TOutputPixel, class TProbabilityPixel>
class AtlasRegistrationMethod : public itk::Object
{
public:

  /** Standard class typedefs. */
  typedef AtlasRegistrationMethod       Self;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  // Image types
  typedef itk::Image<TOutputPixel, 3>          OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::IndexType  OutputImageIndexType;
  typedef typename OutputImageType::OffsetType OutputImageOffsetType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;
  typedef typename OutputImageType::SizeType   OutputImageSizeType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef itk::Image<TProbabilityPixel, 3>          ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer    ProbabilityImagePointer;
  typedef typename ProbabilityImageType::IndexType  ProbabilityImageIndexType;
  typedef typename ProbabilityImageType::OffsetType ProbabilityImageOffsetType;
  typedef typename ProbabilityImageType::PixelType  ProbabilityImagePixelType;
  typedef typename ProbabilityImageType::SizeType   ProbabilityImageSizeType;
  typedef typename ProbabilityImageType::RegionType ProbabilityImageRegionType;

  typedef itk::Image<float, 3>                   InternalImageType;
  typedef typename InternalImageType::Pointer    InternalImagePointer;
  typedef typename InternalImageType::IndexType  InternalImageIndexType;
  typedef typename InternalImageType::OffsetType InternalImageOffsetType;
  typedef typename InternalImageType::PixelType  InternalImagePixelType;
  typedef typename InternalImageType::RegionType InternalImageRegionType;
  typedef typename InternalImageType::SizeType   InternalImageSizeType;

  typedef itk::Image<unsigned char, 3>       ByteImageType;
  typedef typename ByteImageType::Pointer    ByteImagePointer;
  typedef typename ByteImageType::IndexType  ByteImageIndexType;
  typedef typename ByteImageType::OffsetType ByteImageOffsetType;
  typedef typename ByteImageType::PixelType  ByteImagePixelType;
  typedef typename ByteImageType::RegionType ByteImageRegionType;
  typedef typename ByteImageType::SizeType   ByteImageSizeType;

  typedef std::vector<ProbabilityImagePointer> ProbabilityImageList;
  typedef std::vector<OutputImagePointer>      OutputImageList;

  typedef itk::Array<unsigned char> FlagArrayType;

  void SetSuffix(std::string suffix);

  itkGetConstMacro(OutputDebugDir, std::string);
  itkSetMacro(OutputDebugDir, std::string);

  void SetAtlasOriginalImageList(std::vector<InternalImagePointer> & NewAtlasList);

  void SetIntraSubjectOriginalImageList(std::vector<InternalImagePointer> & NewImageList);

  // itkSetMacro( IntraSubjectTransformFileNames, std::vector<std::string> );
  itkSetMacro( AtlasToSubjectTransformFileName, std::string );

  // TODO: KENT:  Move all code from class definition to the .hxx file outside the class definition
  void SetIntraSubjectTransformFileNames(std::vector<std::string> userlist)
  {
    m_IntraSubjectTransformFileNames = userlist;
    m_RegistrationUpdateNeeded = true;
  }

  void RegisterImages();

  GenericTransformType::Pointer GetAtlasToSubjectTransform()
  {
    return m_AtlasToSubjectTransform;
  }

  std::vector<GenericTransformType::Pointer> GetIntraSubjectTransforms()
  {
    return m_IntraSubjectTransforms;
  }

  // Set/Get the Debugging level for filter verboseness
  itkSetMacro(DebugLevel, unsigned int);
  itkGetMacro(DebugLevel, unsigned int);

  itkGetMacro(UseNonLinearInterpolation, bool);
  itkSetMacro(UseNonLinearInterpolation, bool);

  void SetAtlasLinearTransformChoice(const std::string & c)
  {
    m_AtlasLinearTransformChoice = c;
    m_RegistrationUpdateNeeded = true;
  }

  void SetImageLinearTransformChoice(const std::string & c)
  {
    m_ImageLinearTransformChoice = c;
    m_RegistrationUpdateNeeded = true;
  }

  void SetWarpGrid(const unsigned int gx, const unsigned int gy, const unsigned int gz)
  {
    m_WarpGrid.resize(3); m_WarpGrid[0] = gx; m_WarpGrid[1] = gy; m_WarpGrid[2] = gz;
    m_RegistrationUpdateNeeded = true;
  }

  // Get and set the input volume image types.  i.e. T1 or T2 or PD
  void SetInputVolumeTypes(const std::vector<std::string> & newInputVolumeTypes)
  {
    this->m_InputVolumeTypes = newInputVolumeTypes;
    m_RegistrationUpdateNeeded = true;
  }

  void SetAtlasToSubjectInitialTransform( const GenericTransformType::Pointer atlasToSubjectInitialTransform)
  {
    if( this->m_AtlasToSubjectInitialTransform != atlasToSubjectInitialTransform )
      {
      this->m_AtlasToSubjectInitialTransform = atlasToSubjectInitialTransform;
      m_RegistrationUpdateNeeded = true;
      }
  }

  std::vector<std::string> GetInputVolumeTypes(void) const
  {
    return this->m_InputVolumeTypes;
  }

  void Update();

protected:
  void RegisterIntraSubjectImages(void);

  void RegisterAtlasToSubjectImages(void);

  AtlasRegistrationMethod();
  ~AtlasRegistrationMethod();

  OutputImagePointer CopyOutputImage(InternalImagePointer img);

  ProbabilityImagePointer CopyProbabilityImage(InternalImagePointer img);

private:

  std::string m_Suffix;
  std::string m_OutputDebugDir;

  //  ByteImagePointer                  m_AtlasOriginalMask;
  std::vector<InternalImagePointer> m_AtlasOriginalImageList;
  std::vector<InternalImagePointer> m_IntraSubjectOriginalImageList;
  ByteImagePointer                  m_InputImageTissueRegion;
  ImageMaskPointer                  m_InputSpatialObjectTissueRegion;

  std::vector<unsigned int> m_WarpGrid;
  std::vector<std::string>  m_IntraSubjectTransformFileNames;
  std::string               m_AtlasToSubjectTransformFileName;

  std::vector<std::string> m_InputVolumeTypes;

  GenericTransformType::Pointer              m_AtlasToSubjectTransform;
  GenericTransformType::Pointer              m_AtlasToSubjectInitialTransform;
  std::vector<GenericTransformType::Pointer> m_IntraSubjectTransforms;

  bool m_UseNonLinearInterpolation;
  bool m_DoneRegistration;
  bool m_RegistrationUpdateNeeded; // TODO: KENT: The m_RegistrationUpdateNeeded is a hack to replicate the behavior
                                   // that should come from using the modified times of the itk::Object class
                                   //            All the Get/Set functions should use the itkSetMacro so that the
                                   // itk::Object->Modified times are updated correctly, then we can just use that
  //            modify status to determine when re-running is necessary.

  std::string m_AtlasLinearTransformChoice;
  std::string m_ImageLinearTransformChoice;

  unsigned int m_DebugLevel;
};

#ifndef MU_MANUAL_INSTANTIATION
#include "AtlasRegistrationMethod.hxx"
#endif

#endif
