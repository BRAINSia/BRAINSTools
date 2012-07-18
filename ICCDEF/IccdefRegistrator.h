#ifndef _iccdefRegistrator_h
#define _iccdefRegistrator_h

#include "itkObject.h"
#include "itkVector.h"
#include "itkImage.h"
#include "itkMultiResolutionICCDeformableRegistration.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkICCDeformableRegistrationFilter.h"

#include "itkArray.h"

namespace itk
{
/** \class IccdefRegistrator
 *
 * This component computes the transform to register a
 * moving image onto a fixed image.
 *
 * In particular, it uses the deformable demons registration
 * algorithm.
 *
 * The registration is done using a multiresolution strategy.
 * At each resolution level, the downsampled images are obtained
 * using a RecursiveMultiResolutionPyramidImageFilter.
 *
 * \warning This class requires both images to be 3D.
 * It can write out the deformation field and the checker board image
 * of the fixed and output image.
 *
 * The registration process is activated by method Execute().
 *
 * Inputs:
 *   - pointer to fixed image
 *   - pointer to moving image
 *   - number of resolution levels
 *   - number of optimization iterations at each level
 *   - the initial rigid (quaternion) transform parameters
 *   - the coarest level shrink factors for the fixed image
 *   - the coarest level shrink factors for the moving image
 *
 * Outputs:
 *   - output deformation field
 *   - output image
 *   - Checkerboard image
 *   - x,y,z components of displacement fields.
 */
template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue = ITK_TYPENAME TRealImage::PixelType
  >
class IccdefRegistrator : public Object
{
public:

  /** Standard class typedefs. */
  typedef IccdefRegistrator        Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(IccdefRegistrator, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Fixed Image Type. */
  typedef TRealImage RealImageType;

  /** Moving Image Type. */
  typedef TOutputImage                        OutputImageType;
  typedef typename OutputImageType::PixelType PixelType;

  /** Image dimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TRealImage::ImageDimension);

  /** Deformation field value type. */
  typedef TFieldValue FieldValueType;

  /** Deformation field pixel type. */
  typedef Vector<FieldValueType,
                 itkGetStaticConstMacro(ImageDimension)> FieldPixelType;

  /** Deformation field type. */
  typedef Image<FieldPixelType,
                itkGetStaticConstMacro(ImageDimension)> TDeformationField;

  /** Fixed Image Pyramid Type. */
  typedef RecursiveMultiResolutionPyramidImageFilter<
      RealImageType,
      RealImageType>    FixedImagePyramidType;

  /** Moving Image Pyramid Type. */
  typedef RecursiveMultiResolutionPyramidImageFilter<
      RealImageType,
      RealImageType>   MovingImagePyramidType;

  /** Registration Method. */
  typedef MultiResolutionICCDeformableRegistration<
      RealImageType,
      RealImageType,
      TDeformationField>    RegistrationType;

  /** UnsignedIntArray type. */
  typedef Array<unsigned int> UnsignedIntArray;

  /** ShrinkFactorsArray type. */
  typedef FixedArray<unsigned int,
                     itkGetStaticConstMacro(ImageDimension)> ShrinkFactorsArray;

  /** Set the intial deformation field **/
  itkSetObjectMacro(InitialDeformationField, TDeformationField);

  /** Set the fixed image. */
  itkSetObjectMacro( FixedImage, RealImageType );

  /** Set the moving image. */
  itkSetObjectMacro( MovingImage, RealImageType );

  /** Set the Unnormalized moving image. */
  itkSetObjectMacro( UnNormalizedMovingImage, RealImageType );

  /** Set the Unnormalized moving image. */
  itkSetObjectMacro( UnNormalizedFixedImage, RealImageType );

  /** Set the number of resolution levels. */
  itkSetClampMacro( NumberOfLevels, unsigned short, 1,
                    NumericTraits<unsigned short>::max() );

  /** Set the number of iterations per level. */
  itkSetMacro( NumberOfIterations, UnsignedIntArray );

  itkSetStringMacro(InitialFixedDeformationFieldFilename);
  itkGetStringMacro(InitialFixedDeformationFieldFilename);

  itkSetStringMacro(InitialMovingDeformationFieldFilename);
  itkGetStringMacro(InitialMovingDeformationFieldFilename);

  itkSetMacro(OutputJacobianImage, bool);
  itkGetConstMacro(OutputJacobianImage, bool);

  itkSetMacro(OutputDeformationField, bool);
  itkGetConstMacro(OutputDeformationField, bool);

  itkSetMacro(OutputDisplacement, bool);
  itkGetConstMacro(OutputDisplacement, bool);

  /** Set WarpedImageName */
  itkSetStringMacro( OutputPrefix );
  itkGetStringMacro( OutputPrefix );

  /** Set Deformation field output file Name */
  itkSetStringMacro(ForwardDeformationFieldOutputName);
  itkGetStringMacro(ForwardDeformationFieldOutputName);

  itkSetStringMacro(BackwardDeformationFieldOutputName);
  itkGetStringMacro(BackwardDeformationFieldOutputName);

  /**Set histogram matching*/
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);

  /** Method to execute the registration. */
  virtual void Execute();

  /** Get the deformation field. */
  itkGetObjectMacro( DeformationField, TDeformationField );

  /** Initialize registration at the start of new level. */
  void StartNewLevel();

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool );
  itkGetConstMacro( OutDebug,  bool );

  itkSetMacro(DefaultPixelValue, typename RealImageType::PixelType);
  itkGetMacro(DefaultPixelValue, typename RealImageType::PixelType);

  typedef ICCDeformableRegistrationFilter<RealImageType, RealImageType,
                                          TDeformationField> BaseRegistrationFilterType;
  void SetRegistrationFilter(
    BaseRegistrationFilterType * filter)
  {
    this->m_Registration->SetRegistrationFilter(filter);
  }

  RegistrationType * GetRegistrationType(void)
  {
    return m_Registration;
  }

protected:
  IccdefRegistrator();
  ~IccdefRegistrator();
private:
  IccdefRegistrator( const Self & );  // purposely not implemented
  void operator=( const Self & );     // purposely not implemented

  void WriteDisplacementComponents(TDeformationField *, std::string);

  typename TDeformationField::Pointer m_InitialDeformationField;
  typename RealImageType::Pointer m_FixedImage;
  typename RealImageType::Pointer m_MovingImage;
  typename RealImageType::Pointer m_UnNormalizedMovingImage;
  typename RealImageType::Pointer m_UnNormalizedFixedImage;
  typename FixedImagePyramidType::Pointer m_FixedImagePyramid;
  typename MovingImagePyramidType::Pointer m_MovingImagePyramid;
  typename RegistrationType::Pointer m_Registration;
  typename RealImageType::PixelType m_DefaultPixelValue;

  unsigned short   m_NumberOfLevels;
  UnsignedIntArray m_NumberOfIterations;

  typename TDeformationField::Pointer m_DeformationField;
  typename TDeformationField::Pointer m_BackwardDeformationField;
  unsigned long m_Tag;
  std::string   m_DisplacementBaseName;
  std::string   m_OutputPrefix;
  std::string   m_ForwardDeformationFieldOutputName;
  std::string   m_BackwardDeformationFieldOutputName;
  bool          m_OutDebug;
  bool          m_UseHistogramMatching;
  std::string   m_InitialMovingDeformationFieldFilename;
  std::string   m_InitialFixedDeformationFieldFilename;
  bool          m_OutputJacobianImage;
  bool          m_OutputDisplacement;
  bool          m_OutputDeformationField;
  std::string   m_ForwardDir;
  std::string   m_BackwardDir;
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "IccdefRegistrator.txx"
#endif

#endif
