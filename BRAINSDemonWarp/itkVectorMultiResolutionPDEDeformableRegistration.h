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
#ifndef __itkVectorMultiResolutionPDEDeformableRegistration_h
#define __itkVectorMultiResolutionPDEDeformableRegistration_h

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkPDEDeformableRegistrationFilter.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkVectorDiffeomorphicDemonsRegistrationFilter.h"
#include "itkVectorImageToImageAdaptor.h"

#include <vector>

namespace itk
{
/**
 * \class VectorMultiResolutionPDEDeformableRegistration
 * \brief Framework for performing multi-resolution PDE
 * deformable registration ( with linear interpolation and NN extrapolation )
 *
 * VectorMultiResolutionPDEDeformableRegistration3 provides a generic framework
 * to peform multi-resolution deformable registration.
 *
 * At each resolution level a PDEDeformableRegistrationFilter is used
 * to register two images by computing the deformation field which will
 * map a moving image onto a fixed image.
 *
 * A deformation field is represented as an image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the fixed image. The vector type must support element access via operator
 * []. It is assumed that the vector elements behave like floating point
 * scalars.
 *
 * The internal PDEDeformationRegistrationFilter can be set using
 * SetRegistrationFilter. By default a DemonsRegistrationFilter is used.
 *
 * The input fixed and moving images are set via methods SetFixedImage
 * and SetMovingImage respectively. An initial deformation field maybe set via
 * SetInitialDisplacementField or SetInput. If no initial field is set
 * a zero field is used as the initial condition.
 *
 * RecursiveMultiResolutionPyramidImageFilters are used to downsample the fixed
 * and moving images. A VectorExpandImageFilter is used to upsample
 * the deformation as we move from a coarse to fine solution.
 *
 * This class is templated over the fixed image type, the moving image type,
 * and the Displacement Field type.
 *
 *
 * \sa PDEDeformableRegistrationFilter
 * \sa DemonsRegistrationFilter
 * \sa RecursiveMultiResolutionPyramidImageFilter
 * \sa VectorExpandImageFilter
 *
 * The current implementation of this class does not support streaming.
 *
 * \ingroup DeformableImageRegistration
 */
template < typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType = float >
class VectorMultiResolutionPDEDeformableRegistration
  : public ImageToImageFilter< TDisplacementField, TDisplacementField >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( VectorMultiResolutionPDEDeformableRegistration );

  /** Standard class type alias */
  using Self = VectorMultiResolutionPDEDeformableRegistration;
  using Superclass = ImageToImageFilter< TDisplacementField, TDisplacementField >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( VectorMultiResolutionPDEDeformableRegistration, ImageToImageFilter );

  /** Fixed image type. */
  using FixedImageType = TFixedImage;
  using FixedImagePointer = typename FixedImageType::Pointer;
  using FixedImageConstPointer = typename FixedImageType::ConstPointer;
  //  using PixelType = FixedImageType::PixelType;
  //  using VectorImageType = typename itk::VectorImage<TRealType ,3>;

  /** Moving image type. */
  using MovingImageType = TMovingImage;
  using MovingImagePointer = typename MovingImageType::Pointer;
  using MovingImageConstPointer = typename MovingImageType::ConstPointer;

  /** Displacement field image type. */
  using DisplacementFieldType = TDisplacementField;
  using DisplacementFieldPointer = typename DisplacementFieldType::Pointer;

  /** ImageDimension. */
  static constexpr unsigned int ImageDimension = FixedImageType::ImageDimension;

  /** Internal float image type. */
  using FloatImageType = Image< TRealType, Self::ImageDimension >;
  using AdaptorType = itk::VectorImageToImageAdaptor< TRealType, Self::ImageDimension >;

  /** The internal registration type. */
  //  typedef DiffeomorphicDemonsRegistrationFilter2<
  //    FixedImageType, MovingImageType, DisplacementFieldType >
  // RegistrationType;
  //  using RegistrationPointer = typename RegistrationType::Pointer;

  using RegistrationType = PDEDeformableRegistrationFilter< FixedImageType, MovingImageType, DisplacementFieldType >;
  using RegistrationPointer = typename RegistrationType::Pointer;

  using DefaultRegistrationType =
    VectorDiffeomorphicDemonsRegistrationFilter< FixedImageType, MovingImageType, DisplacementFieldType >;
  //  using RegistrationPointer = typename RegistrationType::Pointer;

  /** The default registration type. */

  /** The fixed multi-resolution image pyramid type. */
  using FixedImagePyramidType = RecursiveMultiResolutionPyramidImageFilter< FloatImageType, FloatImageType >;
  using FixedImagePyramidPointer = typename FixedImagePyramidType::Pointer;

  /** The moving multi-resolution image pyramid type. */
  using MovingImagePyramidType = RecursiveMultiResolutionPyramidImageFilter< FloatImageType, FloatImageType >;
  using MovingImagePyramidPointer = typename MovingImagePyramidType::Pointer;

  /** The deformation field expander type. */
  using FieldExpanderType = VectorResampleImageFilter< DisplacementFieldType, DisplacementFieldType >;
  using FieldExpanderPointer = typename FieldExpanderType::Pointer;

  /** Set the fixed image. */
  virtual void
  SetFixedImage( const FixedImageType * ptr );

  //  virtual void SetVectorFixedImage(const VectorImageType * ptr);

  /** Get the fixed image. */
  const FixedImageType *
  GetFixedImage( void ) const;

  //  const VectorImageType * GetVectorFixedImage(void) const;

  /** Set the moving image. */
  virtual void
  SetMovingImage( const MovingImageType * ptr );

  //  virtual void SetVectorMovingImage(const VectorImageType * ptr);

  /** Get the moving image. */
  const MovingImageType *
  GetMovingImage( void ) const;

  //  const VectorImageType * GetVectorMovingImage(void) const;

  /** Set initial deformation field. */
  virtual void
  SetInitialDisplacementField( DisplacementFieldType * ptr )
  {
    // itkExceptionMacro( << "This feature not implemented yet"  );
    // this->SetInput( ptr );
    this->m_InitialDisplacementField = ptr;
  }

  /** Set initial deformation field. No assumption is made on the
   *  input. It will therefore be smoothed and resampled to match the
   *  images characteristics at the coarsest level of the pyramid. */
  virtual void
  SetArbitraryInitialDisplacementField( DisplacementFieldType * ptr )
  {
    this->SetNthInput( 2, ptr );
  }

  /** Get output deformation field. */
  const DisplacementFieldType *
  GetDisplacementField( void )
  {
    return this->GetOutput();
  }

  /** Get the number of valid inputs.  For
   * MultiResolutionPDEDeformableRegistration2, this checks whether the
   * fixed and moving images have been set. While
   * MultiResolutionPDEDeformableRegistration2 can take a third input
   * as an initial deformation field, this input is not a required input.
   */
  std::vector< SmartPointer< DataObject > >::size_type
  GetNumberOfValidRequiredInputs() const override;

  /** Set the internal registrator. */
  itkSetObjectMacro( RegistrationFilter, RegistrationType );

  /** Get the internal registrator. */
  itkGetConstObjectMacro( RegistrationFilter, RegistrationType );

  /** Set the fixed image pyramid. */
  itkSetObjectMacro( FixedImagePyramid, FixedImagePyramidType );

  /** Get the fixed image pyramid. */
  itkGetConstObjectMacro( FixedImagePyramid, FixedImagePyramidType );

  /** Set the moving image pyramid. */
  itkSetObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Get the moving image pyramid. */
  itkGetConstObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Set number of multi-resolution levels. */
  virtual void
  SetNumberOfLevels( unsigned int num );

  /** Get number of multi-resolution levels. */
  itkGetConstReferenceMacro( NumberOfLevels, unsigned int );

  /** Get the current resolution level being processed. */
  itkGetConstReferenceMacro( CurrentLevel, unsigned int );

  /** Set number of iterations per multi-resolution levels. */
  itkSetVectorMacro( NumberOfIterations, unsigned int, m_NumberOfLevels );

  /** Set the moving image pyramid. */
  itkSetObjectMacro( FieldExpander, FieldExpanderType );
  itkGetModifiableObjectMacro( FieldExpander, FieldExpanderType );

  /** Get number of iterations per multi-resolution levels. */
  virtual const unsigned int *
  GetNumberOfIterations() const
  {
    return &( m_NumberOfIterations[0] );
  }

  /** Stop the registration after the current iteration. */
  virtual void
  StopRegistration();

  void
  VerifyInputInformation() const override;

protected:
  VectorMultiResolutionPDEDeformableRegistration();
  ~VectorMultiResolutionPDEDeformableRegistration() override {}

  void
  PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Generate output data by performing the registration
   * at each resolution level. */
  void
  GenerateData() override;

  /** The current implementation of this class does not support
   * streaming. As such it requires the largest possible region
   * for the moving, fixed and input deformation field. */
  void
  GenerateInputRequestedRegion() override;

  /** By default, the output deformation field has the same
   * spacing, origin and LargestPossibleRegion as the input/initial
   * deformation field.
   *
   * If the initial deformation field is not set, the output
   * information is copied from the fixed image. */
  void
  GenerateOutputInformation() override;

  /** The current implementation of this class does not supprot
   * streaming. As such it produces the output for the largest
   * possible region. */
  void
  EnlargeOutputRequestedRegion( DataObject * ptr ) override;

  /** This method returns true to indicate that the registration should
   * terminate at the current resolution level. */
  virtual bool
  Halt();

private:
  RegistrationPointer                      m_RegistrationFilter;
  FixedImagePyramidPointer                 m_FixedImagePyramid;
  MovingImagePyramidPointer                m_MovingImagePyramid;
  std::vector< FixedImagePyramidPointer >  m_FixedVectorImagePyramid;
  std::vector< MovingImagePyramidPointer > m_MovingVectorImagePyramid;
  FieldExpanderPointer                     m_FieldExpander;

  unsigned int                m_NumberOfLevels;
  unsigned int                m_CurrentLevel;
  std::vector< unsigned int > m_NumberOfIterations;

  /** Flag to indicate user stop registration request. */
  bool                     m_StopRegistrationFlag;
  DisplacementFieldPointer m_InitialDisplacementField;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkVectorMultiResolutionPDEDeformableRegistration.hxx"
#endif

#endif
