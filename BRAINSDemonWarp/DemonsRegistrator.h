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
#ifndef __DemonsRegistrator_h
#define __DemonsRegistrator_h

#include "itkObject.h"
#include "itkVector.h"
#include "itkImage.h"
#include "itkMultiResolutionPDEDeformableRegistration.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"

#include "itkArray.h"

namespace itk
{
/** \class DemonsRegistrator
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
  class TFieldValue = typename TRealImage::PixelType
  >
class DemonsRegistrator : public Object
{
public:

  /** Standard class typedefs. */
  typedef DemonsRegistrator        Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(DemonsRegistrator, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Fixed Image Type. */
  typedef TRealImage RealImageType;

  /** Moving Image Type. */
  typedef typename TOutputImage::PixelType PixelType;

  /** Image dimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TRealImage::ImageDimension);

  /** Type to hold the number of checker boxes per dimension */
  typedef FixedArray<unsigned int, TRealImage::ImageDimension> PatternArrayType;

  /** Set Checker pattern */
  itkSetMacro(CheckerBoardPattern, PatternArrayType);
  itkGetConstReferenceMacro(CheckerBoardPattern, PatternArrayType);

  /** Displacement field value type. */
  typedef TFieldValue FieldValueType;

  /** Displacement field pixel type. */
  typedef Vector<FieldValueType,
                 itkGetStaticConstMacro(ImageDimension)> FieldPixelType;

  /** Displacement field type. */
  typedef Image<FieldPixelType,
                itkGetStaticConstMacro(ImageDimension)> TDisplacementField;

  /** Fixed Image Pyramid Type. */
  typedef RecursiveMultiResolutionPyramidImageFilter<
      RealImageType,
      RealImageType>    FixedImagePyramidType;

  /** Moving Image Pyramid Type. */
  typedef RecursiveMultiResolutionPyramidImageFilter<
      RealImageType,
      RealImageType>   MovingImagePyramidType;

  /** Registration Method. */
  typedef MultiResolutionPDEDeformableRegistration<
      RealImageType,
      RealImageType,
      TDisplacementField, float>    RegistrationType;

  /** UnsignedIntArray type. */
  typedef Array<unsigned int> UnsignedIntArray;

  /** ShrinkFactorsArray type. */
  typedef FixedArray<unsigned int,
                     itkGetStaticConstMacro(ImageDimension)> ShrinkFactorsArray;

  /** Set the intial deformation field */
  itkSetObjectMacro(InitialDisplacementField, TDisplacementField);

  /** Set the fixed image. */
  itkSetObjectMacro(FixedImage, RealImageType);

  /** Set the moving image. */
  itkSetObjectMacro(MovingImage, RealImageType);

  /** Set the Unnormalized moving image. */
  itkSetObjectMacro(UnNormalizedMovingImage, RealImageType);

  /** Set the Unnormalized moving image. */
  itkSetObjectMacro(UnNormalizedFixedImage, RealImageType);

  /** Set the number of resolution levels. */
  itkSetClampMacro( NumberOfLevels, unsigned short, 1,
                    NumericTraits<unsigned short>::max() );

  /** Set the number of iterations per level. */
  itkSetMacro(NumberOfIterations, UnsignedIntArray);

  /** Set the fixed and moving image shrink factors. */
  itkSetMacro(FixedImageShrinkFactors, ShrinkFactorsArray);
  itkSetMacro(MovingImageShrinkFactors, ShrinkFactorsArray);

  /** Set Displacementname */
  itkSetStringMacro(DisplacementBaseName);
  itkGetStringMacro(DisplacementBaseName);

  /** Set WarpedImageName */
  itkSetStringMacro(WarpedImageName);
  itkGetStringMacro(WarpedImageName);

  /** Set CheckerBoard ImageName */
  itkSetStringMacro(CheckerBoardFilename);
  itkGetStringMacro(CheckerBoardFilename);

  /** Set Displacement field output file Name */
  itkSetStringMacro(DisplacementFieldOutputName);
  itkGetStringMacro(DisplacementFieldOutputName);

  /**Set histogram matching */
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);

  /** Method to execute the registration. */
  virtual void Execute();

  /** Get the deformation field. */
  itkGetConstObjectMacro(DisplacementField, TDisplacementField);

  /** Initialize registration at the start of new level. */
  void StartNewLevel();

  /** Output Normalized Image.*/
  itkSetStringMacro(OutNormalized);
  itkGetStringMacro(OutNormalized);

  /** Set Debug mode */
  itkSetMacro(OutDebug, bool);
  itkGetConstMacro(OutDebug,  bool);

  itkSetStringMacro(FixedLandmarkFilename);
  itkGetStringMacro(FixedLandmarkFilename);
  itkSetStringMacro(MovingLandmarkFilename);
  itkGetStringMacro(MovingLandmarkFilename);

  itkSetMacro(DefaultPixelValue, typename RealImageType::PixelType);
  itkGetMacro(DefaultPixelValue, typename RealImageType::PixelType);

  /** Get the interpolation Mode. */
  itkGetMacro(InterpolationMode, std::string);
  itkSetMacro(InterpolationMode, std::string);

  typedef itk::PDEDeformableRegistrationFilter<RealImageType, RealImageType,
                                               TDisplacementField> BaseRegistrationFilterType;
  void SetRegistrationFilter(BaseRegistrationFilterType *filter)
  {
    this->m_Registration->SetRegistrationFilter(filter);
  }

  RegistrationType * GetRegistrationType(void)
  {
    return m_Registration;
  }

protected:
  DemonsRegistrator();
  ~DemonsRegistrator();
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DemonsRegistrator);

  void WriteDisplacementComponents();

  typename TDisplacementField::Pointer m_InitialDisplacementField;
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

  ShrinkFactorsArray m_MovingImageShrinkFactors;
  ShrinkFactorsArray m_FixedImageShrinkFactors;

  typename TDisplacementField::Pointer m_DisplacementField;
  std::string      m_FixedLandmarkFilename;
  std::string      m_MovingLandmarkFilename;
  unsigned long    m_Tag;
  std::string      m_DisplacementBaseName;
  std::string      m_WarpedImageName;
  std::string      m_CheckerBoardFilename;
  std::string      m_DisplacementFieldOutputName;
  PatternArrayType m_CheckerBoardPattern;
  std::string      m_OutNormalized;
  bool             m_OutDebug;
  bool             m_UseHistogramMatching;
  std::string      m_InterpolationMode;
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "DemonsRegistrator.hxx"
#endif

#endif
