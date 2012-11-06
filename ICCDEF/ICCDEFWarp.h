#ifndef _ICCDEFWarp_h
#define _ICCDEFWarp_h

#include <string>
#include "IccdefPreprocessor.h"
#include "IccdefRegistrator.h"
#include "ICCApplicationBase.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkICCDeformableRegistrationFilter.h"

namespace itk
{
/*This file defines ICCDEF registration class which initializes the input
  preprocessor and the registrator. */

template <typename TImage,
          typename TRealImage,
          typename TOutputImage
          >
class ICCDEFWarp : public ICCApplicationBase<
    IccdefPreprocessor<TImage, TRealImage>,
    IccdefRegistrator<TRealImage, TOutputImage, typename TRealImage::PixelType>
    >
{
public:

  /** Standard class typedefs. */
  typedef ICCDEFWarp Self;
  typedef ICCApplicationBase<
      IccdefPreprocessor<TImage, TRealImage>,
      IccdefRegistrator<TRealImage, TRealImage,
                        typename TRealImage::PixelType>
      > Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Deformation field pixel type. */
  typedef float                     FieldValueType;
  typedef Vector<FieldValueType, 3> FieldPixelType;
  typedef Image<FieldPixelType, 3>  TDisplacementField;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ICCDEFWarp, ICCApplicationBase);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image types. */
  typedef TImage     ImageType;
  typedef TRealImage RealImageType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImage::ImageDimension);

  /** Type to hold the number of checker boxes per dimension */
  typedef FixedArray<unsigned int, TImage::ImageDimension> PatternArrayType;

  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::SizeType  SizeType;

  /** ShrinkFactors type. */
  typedef FixedArray<unsigned int, TImage::ImageDimension> ShrinkFactorsType;

  /** IterationArray type. */
  typedef Array<unsigned int> IterationsArrayType;

  /** Set the atlas patient ID. */
  itkSetStringMacro(TheMovingImageFilename);
  itkGetStringMacro(TheMovingImageFilename);

  /** Set the subject patient ID. */
  itkSetStringMacro(TheFixedImageFilename);
  itkGetStringMacro(TheFixedImageFilename);

  /** Set the initial Displacement Field one of 3 ways. */
//  itkSetStringMacro (InitialDisplacementFieldFilename);
//  itkGetStringMacro (InitialDisplacementFieldFilename);

  itkSetStringMacro(InitialFixedDisplacementFieldFilename);
  itkGetStringMacro(InitialFixedDisplacementFieldFilename);

  itkSetStringMacro(InitialMovingDisplacementFieldFilename);
  itkGetStringMacro(InitialMovingDisplacementFieldFilename);

  itkSetStringMacro(InitialCoefficientFilename);
  itkGetStringMacro(InitialCoefficientFilename);

  itkSetStringMacro(InitialTransformFilename);
  itkGetStringMacro(InitialTransformFilename);

  /** Set Displacementname */
//  itkSetStringMacro (DisplacementBaseName);
//  itkGetStringMacro (DisplacementBaseName);
/** Set WarpedImageName */
  itkSetStringMacro(OutputPrefix);
  itkGetStringMacro(OutputPrefix);

  /** Set input parameter file */
  // itkSetStringMacro (ParameterFilename);

  /** Set output transformation filename. */
  // itkSetStringMacro (OutputFilename);

  /**Set Deformation field output filename*/
  itkSetStringMacro(ForwardDisplacementFieldOutputName);
  itkGetStringMacro(ForwardDisplacementFieldOutputName);

  itkSetStringMacro(BackwardDisplacementFieldOutputName);
  itkGetStringMacro(BackwardDisplacementFieldOutputName);
  /**Set Jacobian Image prefix name*/
  itkSetMacro(OutputJacobianImage, bool);
  itkGetConstMacro(OutputJacobianImage, bool);

  itkSetMacro(OutputDisplacementField, bool);
  itkGetConstMacro(OutputDisplacementField, bool);

  itkSetMacro(OutputDisplacement, bool);
  itkGetConstMacro(OutputDisplacement, bool);

  /** Set append output file boolean. */
  itkSetMacro(AppendOutputFile, bool);
  itkGetMacro(AppendOutputFile, bool);
  itkBooleanMacro(AppendOutputFile);

  /**Force Centered Image.*/
  itkSetMacro(ForceCoronalZeroOrigin, bool);
  itkGetConstMacro(ForceCoronalZeroOrigin, bool);

  /*BOBF macros */
  /**Set Target Mask filename*/
  itkSetStringMacro(FixedBinaryVolume);
  itkGetStringMacro(FixedBinaryVolume);

  /**Set Template Mask filename*/
  itkSetStringMacro(MovingBinaryVolume);
  itkGetStringMacro(MovingBinaryVolume);

#if 0
  /** Set/Get the lower threshold. The default is 0. */
  itkSetMacro(Lower, PixelType);
  itkGetMacro(Lower, PixelType);

  /** Set/Get the upper threshold. The default is 70 */
  itkSetMacro(Upper, PixelType);
  itkGetMacro(Upper, PixelType);

  /** Set the Seed of the neighborhood used for a mask. */
  itkSetMacro(Seed, IndexType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Seed, IndexType);

  /** Set the radius of the neighborhood used for a mask. */
  itkSetMacro(Radius, SizeType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Radius, SizeType);
#endif

  /** Set/Get value to replace thresholded pixels. Pixels that lie *
   *  within Lower and Upper (inclusive) will be replaced with this
   *  value. The default is 1. */
  itkSetMacro(DefaultPixelValue, PixelType);
  itkGetMacro(DefaultPixelValue, PixelType);

  itkSetMacro(MedianFilterSize,  SizeType);
  itkGetMacro(MedianFilterSize,  SizeType);

  /** Set the initial deformation field to prime registration */
  //    itkSetObjectMacro(InitialDisplacementField,TDisplacementField);
  /** Set the Input Landmark Filename*/

  /**Set histogram matching*/
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);

  /** Get the number of histogram bins. */
  itkGetConstMacro( NumberOfHistogramLevels, unsigned long );
  itkSetMacro( NumberOfHistogramLevels, unsigned long );

  /** Get the number of match points. */
  itkGetConstMacro( NumberOfMatchPoints, unsigned long );
  itkSetMacro( NumberOfMatchPoints, unsigned long );

  /** Get the number of levels. */
  itkGetMacro( NumberOfLevels, unsigned short );
  itkSetMacro( NumberOfLevels, unsigned short );

  /** Get the number of iterations at each level. */
  itkGetConstReferenceMacro( NumberOfIterations, IterationsArrayType );
  void SetNumberOfIterations(const IterationsArrayType & iterations)
  {
    m_NumberOfIterations = iterations;
  }

  typedef ICCDeformableRegistrationFilter<RealImageType, RealImageType,
                                          TDisplacementField>
    BaseRegistrationFilterType;
  void SetRegistrationFilter(
    BaseRegistrationFilterType * filter)
  {
    this->m_Registrator->SetRegistrationFilter(filter);
  }

protected:

  ICCDEFWarp();
  virtual ~ICCDEFWarp()
  {
  }

  /*** Initialize the preprocessor */
  virtual void InitializePreprocessor();

  /*** Initialize the registrator  */
  virtual void InitializeRegistrator();

private:

  std::string m_TheMovingImageFilename;
  std::string m_TheFixedImageFilename;
  std::string m_InitialDisplacementFieldFilename;
  std::string m_InitialMovingDisplacementFieldFilename;
  std::string m_InitialFixedDisplacementFieldFilename;
  std::string m_InitialCoefficientFilename;
  std::string m_InitialTransformFilename;

//  std::string m_WarpedImageName;
  std::string m_OutputPrefix;

  std::string m_FixedBinaryVolume;
  std::string m_MovingBinaryVolume;

  // std::string m_ParameterFilename;
  bool m_ForceCoronalZeroOrigin;
  bool m_UseHistogramMatching;
  bool m_OutputDisplacement;
  bool m_OutputDisplacementField;
  bool m_OutputJacobianImage;

  std::string m_OutNormalized;
  //    std::string m_OutDebug;
  std::string m_OutputFilename;
  std::string m_ForwardDisplacementFieldOutputName;
  std::string m_BackwardDisplacementFieldOutputName;
  bool        m_AppendOutputFile;

  PixelType m_DefaultPixelValue;
  SizeType  m_MedianFilterSize;

  // typename TDisplacementField::Pointer m_InitialDisplacementField;
  unsigned long       m_NumberOfHistogramLevels;
  unsigned long       m_NumberOfMatchPoints;
  unsigned short      m_NumberOfLevels;
  IterationsArrayType m_NumberOfIterations;
};
}          // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ICCDEFWarp.hxx"
#endif
#endif
