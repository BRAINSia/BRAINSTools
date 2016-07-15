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
  class TFieldValue = typename TRealImage::PixelType
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
  typedef MultiResolutionICCDeformableRegistration<
      RealImageType,
      RealImageType,
      TDisplacementField>    RegistrationType;

  /** UnsignedIntArray type. */
  typedef Array<unsigned int> UnsignedIntArray;

  /** ShrinkFactorsArray type. */
  typedef FixedArray<unsigned int,
                     itkGetStaticConstMacro(ImageDimension)> ShrinkFactorsArray;

  /** Set the intial deformation field **/
  itkSetObjectMacro(InitialDisplacementField, TDisplacementField);

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

  itkSetStringMacro(InitialFixedDisplacementFieldFilename);
  itkGetStringMacro(InitialFixedDisplacementFieldFilename);

  itkSetStringMacro(InitialMovingDisplacementFieldFilename);
  itkGetStringMacro(InitialMovingDisplacementFieldFilename);

  itkSetMacro(OutputJacobianImage, bool);
  itkGetConstMacro(OutputJacobianImage, bool);

  itkSetMacro(OutputDisplacementField, bool);
  itkGetConstMacro(OutputDisplacementField, bool);

  itkSetMacro(OutputDisplacement, bool);
  itkGetConstMacro(OutputDisplacement, bool);

  /** Set WarpedImageName */
  itkSetStringMacro( OutputPrefix );
  itkGetStringMacro( OutputPrefix );

  /** Set Deformation field output file Name */
  itkSetStringMacro(ForwardDisplacementFieldOutputName);
  itkGetStringMacro(ForwardDisplacementFieldOutputName);

  itkSetStringMacro(BackwardDisplacementFieldOutputName);
  itkGetStringMacro(BackwardDisplacementFieldOutputName);

  /**Set histogram matching*/
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);

  /** Method to execute the registration. */
  virtual void Execute();

  /** Get the deformation field. */
  itkGetConstObjectMacro( DisplacementField, TDisplacementField );

  /** Initialize registration at the start of new level. */
  void StartNewLevel();

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool );
  itkGetConstMacro( OutDebug,  bool );

  itkSetMacro(DefaultPixelValue, typename RealImageType::PixelType);
  itkGetMacro(DefaultPixelValue, typename RealImageType::PixelType);

  typedef ICCDeformableRegistrationFilter<RealImageType, RealImageType,
                                          TDisplacementField> BaseRegistrationFilterType;
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
  ITK_DISALLOW_COPY_AND_ASSIGN(IccdefRegistrator);

  void WriteDisplacementComponents(TDisplacementField *, std::string);

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

  typename TDisplacementField::Pointer m_DisplacementField;
  typename TDisplacementField::Pointer m_BackwardDisplacementField;
  unsigned long m_Tag;
  std::string   m_DisplacementBaseName;
  std::string   m_OutputPrefix;
  std::string   m_ForwardDisplacementFieldOutputName;
  std::string   m_BackwardDisplacementFieldOutputName;
  bool          m_OutDebug;
  bool          m_UseHistogramMatching;
  std::string   m_InitialMovingDisplacementFieldFilename;
  std::string   m_InitialFixedDisplacementFieldFilename;
  bool          m_OutputJacobianImage;
  bool          m_OutputDisplacement;
  bool          m_OutputDisplacementField;
  std::string   m_ForwardDir;
  std::string   m_BackwardDir;
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "IccdefRegistrator.hxx"
#endif

#endif
