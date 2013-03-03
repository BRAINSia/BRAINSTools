#ifndef _IccdefPreprocessor_h
#define _IccdefPreprocessor_h

#include "itkObject.h"
#include "itkAffineTransform.h"
#include "itkScaleVersor3DTransform.h"

#include "itkVersorRigid3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkBSplineDeformableTransform.h"

namespace itk
{
/** \class IccdefPreprocessor
 *
 * This component pre-processes the moving and fixed image before
 * registration.
 * If the fixed image dimensions are different from the moving image it will     * resample the moving image to match the fixed image dimensions.
 * Histogram matching is done to solve the intensity mismatch problem.
 *
 * The preprocessor also called the skull-stripping filter itkBOBF
 * if an atlas and subject whole brain masks are specified.
 *
 * The preprocessing is activatived by method Execute().
 *
 * Inputs:
 *    - pointer to original fixed image
 *    - pointer original moving image
 *    - number of histogram levels
 *    - number of match points
 *
 * Outputs:
 *    - pointer to processed fixed image
 *    - pointer to processed moving image
 *    - the minimum value of original fixed image
 *    - the minimum value of original moving image
 *
 */
template <typename TInputImage, typename TOutputImage>
class ITK_EXPORT IccdefPreprocessor : public Object
{
public:

  /** Standard class typedefs. */
  typedef IccdefPreprocessor       Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(IccdefPreprocessor, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Input Image Type. */
  typedef TInputImage InputImageType;
  /** Output Image Type. */
  typedef TOutputImage OutputImageType;

  /** Input image pixel type. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType PixelType;
  typedef typename OutputImageType::IndexType IndexType;
  typedef typename OutputImageType::SizeType  SizeType;

  /** Image dimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Set the atlas patient ID. */
  itkSetStringMacro(TheMovingImageFilename);
  itkGetStringMacro(TheMovingImageFilename);

  /** Set the subject patient ID. */
  itkSetStringMacro(TheFixedImageFilename);
  itkGetStringMacro(TheFixedImageFilename);

  /** Deformation field value type. */
  typedef float FieldValueType;

  typedef itk::Image<unsigned char, itkGetStaticConstMacro(ImageDimension)> TImageType;

  /** Deformation field pixel type. */
  typedef Vector<FieldValueType,
                 itkGetStaticConstMacro(ImageDimension)> FieldPixelType;

  /** Deformation field type. */
  typedef Image<FieldPixelType,
                itkGetStaticConstMacro(ImageDimension)> TDisplacementField;

  /** Set the initial Deformation Field. */
  itkSetObjectMacro( InitialDisplacementField, TDisplacementField );
  itkGetModifiableObjectMacro( InitialDisplacementField, TDisplacementField );

  /** Set the number of histogram levels to use. */
  itkSetMacro( NumberOfHistogramLevels, unsigned long );

  /** Set the number of match points to use. */
  itkSetMacro( NumberOfMatchPoints, unsigned long );

  /** Method to execute the preprocessing. */
  virtual void Execute();

  /** Get the output fixed image. */
  itkGetModifiableObjectMacro( OutputFixedImage, OutputImageType );
  itkGetModifiableObjectMacro( OutputMovingImage, OutputImageType );

  /** Get the output moving image. */
  itkGetModifiableObjectMacro( UnNormalizedMovingImage, OutputImageType );
  itkGetModifiableObjectMacro( UnNormalizedFixedImage, OutputImageType );

  /** Get minimum value of original fixed image. */
  itkGetMacro( FixedImageMinimum, InputPixelType );

  /** Get minimum value of original moving image. */
  itkGetMacro( MovingImageMinimum, InputPixelType );

  /** Transform Types. */
  typedef VersorRigid3DTransform<double>                                  VersorRigid3DTransformType;
  typedef ScaleSkewVersor3DTransform<double>                              ScaleSkewVersor3DTransformType;
  typedef AffineTransform<double, itkGetStaticConstMacro(ImageDimension)> AffineTransformType;

  /** Image dimension enumeration. */
  itkStaticConstMacro(SplineOrder, unsigned int, 3);

  typedef double CoordinateRepType;
  typedef typename itk::BSplineDeformableTransform<
      CoordinateRepType,
      itkGetStaticConstMacro(ImageDimension),
      itkGetStaticConstMacro(SplineOrder)> BSplineTransformType;

  /*BOBF macros*/
  /**Set Target Mask filename*/
  itkSetStringMacro(FixedBinaryVolume );
  itkGetStringMacro( FixedBinaryVolume );

  /**Set Template Mask filename*/
  itkSetStringMacro(MovingBinaryVolume );
  itkGetStringMacro( MovingBinaryVolume );

#if 0
  /** Set/Get the lower threshold. The default is 0. */
  itkSetMacro(Lower, PixelType);
  itkGetMacro(Lower, PixelType);

  /** Set/Get the upper threshold. The default is 70 */
  itkSetMacro(Upper, PixelType);
  itkGetMacro(Upper, PixelType);

  /** Set the radius of the neighborhood used for a mask. */
  itkSetMacro(Radius, SizeType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Radius, SizeType);

  /** Set the Seed of the neighborhood used for a mask. */
  itkSetMacro(Seed, IndexType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Seed, IndexType);
#endif

  itkSetMacro(DefaultPixelValue,  PixelType);
  itkGetMacro(DefaultPixelValue,  PixelType);

  itkSetMacro(MedianFilterSize,  SizeType);
  itkGetMacro(MedianFilterSize,  SizeType);

  itkSetStringMacro( InitialTransformFilename );

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool);
  itkGetConstMacro(OutDebug, bool);

  /**Set histogram matching*/
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);
protected:
  IccdefPreprocessor();
  ~IccdefPreprocessor()
  {
  };
private:
  IccdefPreprocessor( const Self & );  // purposely not implemented
  void operator=( const Self & );      // purposely not implemented

  typename InputImageType::Pointer m_InputFixedImage;
  typename InputImageType::Pointer m_InputMovingImage;
  typename OutputImageType::Pointer m_OutputFixedImage;
  typename OutputImageType::Pointer m_OutputMovingImage;
  typename OutputImageType::Pointer m_UnNormalizedMovingImage;
  typename OutputImageType::Pointer m_UnNormalizedFixedImage;
  typename TDisplacementField::Pointer m_InitialDisplacementField;

  unsigned long m_NumberOfHistogramLevels;
  unsigned long m_NumberOfMatchPoints;

  InputPixelType m_FixedImageMinimum;
  InputPixelType m_MovingImageMinimum;

  std::string m_FixedBinaryVolume;
  std::string m_MovingBinaryVolume;
//  IndexType   m_Seed;
//  PixelType   m_Lower;
//  PixelType   m_Upper;
//  SizeType    m_Radius;
  PixelType m_DefaultPixelValue;
  bool      m_OutDebug;
  SizeType  m_MedianFilterSize;
  typedef typename OutputImageType::Pointer OutputImagePointer;
  typedef typename InputImageType::Pointer  InputImagePointer;

  bool        m_UseHistogramMatching;
  std::string m_TheMovingImageFilename;
  std::string m_TheFixedImageFilename;
  std::string m_InitialTransformFilename;

  /*MakeBOBF function takes in a brain image and a whole brain mask and strips
    the skull of the image.*/
  InputImagePointer MakeBOBFImage( InputImagePointer input, std::string MaskName );
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "IccdefPreprocessor.hxx"
#endif

#endif // _IccdefPreprocessor_h
