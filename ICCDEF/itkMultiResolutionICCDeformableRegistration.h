/*
 *  itkMultiResolutionICCDeformableRegistration.h
 *  iccdefRegistrationNew
 *
 *  Created by Yongqiang Zhao on 5/8/09.
 *  Copyright 2009 UI. All rights reserved.
 *
 */

#ifndef __itkMultiResolutionICCDeformableRegistration_h
#define __itkMultiResolutionICCDeformableRegistration_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkPDEDeformableRegistrationFilter.h"
#include "itkICCDeformableRegistrationFilter.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkVectorResampleImageFilter.h"

#include <vector>

namespace itk
{
/**
 * \class MultiResolutionICCDeformableRegistration
 * \brief Framework for perfoming multi-resolution PDE deformable registration.
 *
 * MultiResolutionICCDeformableRegistration provides a generic framework
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
 * MultiResolutionPyramidImageFilters are used to downsample the fixed
 * and moving images. A VectorExpandImageFilter is used to upsample
 * the deformation as we move from a coarse to fine solution.
 *
 * This class is templated over the fixed image type, the moving image type,
 * and the Deformation Field type.
 *
 * \warning This class assumes that the fixed, moving and deformation
 * field image types all have the same number of dimensions.
 *
 * \sa PDEDeformableRegistrationFilter
 * \sa DemonsRegistrationFilter
 * \sa MultiResolutionPyramidImageFilter
 * \sa VectorExpandImageFilter
 *
 * The current implementation of this class does not support streaming.
 *
 * \ingroup DeformableImageRegistration
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
class ITK_EXPORT MultiResolutionICCDeformableRegistration :
  public         ImageToImageFilter<TDisplacementField, TDisplacementField>
{
public:
  /** Standard class typedefs */
  typedef MultiResolutionICCDeformableRegistration Self;
  typedef ImageToImageFilter<TDisplacementField, TDisplacementField>
    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiResolutionICCDeformableRegistration,
                ImageToImageFilter );

  /** Fixed image type. */
  typedef TFixedImage                           FixedImageType;
  typedef typename FixedImageType::Pointer      FixedImagePointer;
  typedef typename FixedImageType::ConstPointer FixedImageConstPointer;

  /** Moving image type. */
  typedef TMovingImage                           MovingImageType;
  typedef typename MovingImageType::Pointer      MovingImagePointer;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointer;

  /** Deformation field image type. */
  typedef TDisplacementField                      DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer DisplacementFieldPointer;

  typedef typename DataObject::Pointer DataObjectPointer;

  /** ImageDimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      FixedImageType::ImageDimension);

  /** Internal float image type. */
  typedef Image<float, itkGetStaticConstMacro(ImageDimension)> FloatImageType;

  /** The internal registration type. */
  typedef ICCDeformableRegistrationFilter<
      FloatImageType, FloatImageType, DisplacementFieldType> RegistrationType;
  typedef typename RegistrationType::Pointer RegistrationPointer;

  /** The default registration type. */
  typedef DemonsRegistrationFilter<
      FloatImageType, FloatImageType, DisplacementFieldType> DefaultRegistrationType;

  /** The fixed multi-resolution image pyramid type. */
  typedef MultiResolutionPyramidImageFilter<
      FixedImageType, FloatImageType> FixedImagePyramidType;
  typedef typename FixedImagePyramidType::Pointer FixedImagePyramidPointer;

  /** The moving multi-resolution image pyramid type. */
  typedef MultiResolutionPyramidImageFilter<
      MovingImageType, FloatImageType> MovingImagePyramidType;
  typedef typename MovingImagePyramidType::Pointer MovingImagePyramidPointer;

  /** The deformation field expander type. */
  typedef VectorResampleImageFilter<
      DisplacementFieldType, DisplacementFieldType> FieldExpanderType;
  typedef typename FieldExpanderType::Pointer FieldExpanderPointer;

  /** Set the fixed image. */
  virtual void SetFixedImage( const FixedImageType * ptr );

  /** Get the fixed image. */
  const FixedImageType * GetFixedImage(void) const;

  /** Set the moving image. */
  virtual void SetMovingImage( const MovingImageType * ptr );

  /** Get the moving image. */
  const MovingImageType * GetMovingImage(void) const;

  /** Set initial deformation field. */
  virtual void SetInitialDisplacementField( DisplacementFieldType * ptr )
  {
    this->m_InitialDisplacementField = ptr;
    // itkExceptionMacro( << "This feature not implemented yet"  );
  }

  virtual void SetInitialFixedDisplacementField( DisplacementFieldType * ptr )
  {
    this->m_InitialFixedDisplacementField = ptr;
    // itkExceptionMacro( << "This feature not implemented yet"  );
  }

  virtual void SetInitialMovingDisplacementField( DisplacementFieldType * ptr )
  {
    this->m_InitialMovingDisplacementField = ptr;
    // itkExceptionMacro( << "This feature not implemented yet"  );
  }

  /** Get output deformation field. */
  const DisplacementFieldType * GetDisplacementField(void)
  {
    return this->GetOutput();
  }

  DisplacementFieldPointer GetBackwardDisplacementField(void)
  {
    return m_BackwardDisplacementField;
  }

  /** Get the number of valid inputs.  For
   * MultiResolutionICCDeformableRegistration, this checks whether the
   * fixed and moving images have been set. While
   * MultiResolutionICCDeformableRegistration can take a third input
   * as an initial deformation field, this input is not a required input.
   */
  virtual std::vector<SmartPointer<DataObject> >::size_type GetNumberOfValidRequiredInputs() const;

  /** Set the internal registrator. */
  itkSetObjectMacro( RegistrationFilter, RegistrationType );

  /** Get the internal registrator. */
  itkGetObjectMacro( RegistrationFilter, RegistrationType );

  /** Set the fixed image pyramid. */
  itkSetObjectMacro( FixedImagePyramid, FixedImagePyramidType );

  /** Get the fixed image pyramid. */
  itkGetObjectMacro( FixedImagePyramid, FixedImagePyramidType );

  /** Set the moving image pyramid. */
  itkSetObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Get the moving image pyramid. */
  itkGetObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Set number of multi-resolution levels. */
  virtual void SetNumberOfLevels( unsigned int num );

  /** Get number of multi-resolution levels. */
  itkGetConstReferenceMacro( NumberOfLevels, unsigned int );

  /** Get the current resolution level being processed. */
  itkGetConstReferenceMacro( CurrentLevel, unsigned int );

  /** Set number of iterations per multi-resolution levels. */
  itkSetVectorMacro( NumberOfIterations, unsigned int, m_NumberOfLevels );

  /** Get number of iterations per multi-resolution levels. */
  virtual const unsigned int * GetNumberOfIterations() const
  {
    return &(m_NumberOfIterations[0]);
  }

  /** Stop the registration after the current iteration. */
  virtual void StopRegistration();

  /** Set the moving image pyramid. */
  itkSetObjectMacro( FieldExpander12, FieldExpanderType );

  /** Get the moving image pyramid. */
  itkGetObjectMacro( FieldExpander12, FieldExpanderType );

  /** Set the moving image pyramid. */
  itkSetObjectMacro( FieldExpander21, FieldExpanderType );

  /** Get the moving image pyramid. */
  itkGetObjectMacro( FieldExpander21, FieldExpanderType );

  using Superclass::MakeOutput;
  virtual ProcessObject::DataObjectPointer MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx);

  itkSetStringMacro(DisplacementFieldOutputNamePrefix);
  itkGetStringMacro(DisplacementFieldOutputNamePrefix);
protected:
  MultiResolutionICCDeformableRegistration();
  ~MultiResolutionICCDeformableRegistration()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Generate output data by performing the registration
   * at each resolution level. */
  virtual void GenerateData();

  /** The current implementation of this class does not support
   * streaming. As such it requires the largest possible region
   * for the moving, fixed and input deformation field. */
  virtual void GenerateInputRequestedRegion();

  /** By default, the output deformation field has the same
   * spacing, origin and LargestPossibleRegion as the input/initial
   * deformation field.
   *
   * If the initial deformation field is not set, the output
   * information is copied from the fixed image. */
  virtual void GenerateOutputInformation();

  /** The current implementation of this class does not supprot
   * streaming. As such it produces the output for the largest
   * possible region. */
  virtual void EnlargeOutputRequestedRegion( DataObject *ptr );

  /** This method returns true to indicate that the registration should
   * terminate at the current resolution level. */
  virtual bool Halt();

private:
  MultiResolutionICCDeformableRegistration(const Self &); // purposely not implemented
  void operator=(const Self &);                           // purposely not implemented

  RegistrationPointer       m_RegistrationFilter;
  FixedImagePyramidPointer  m_FixedImagePyramid;
  MovingImagePyramidPointer m_MovingImagePyramid;
  FieldExpanderPointer      m_FieldExpander12, m_FieldExpander21;
  DisplacementFieldPointer  m_InitialDisplacementField;
  DisplacementFieldPointer  m_InitialFixedDisplacementField;
  DisplacementFieldPointer  m_InitialMovingDisplacementField;
  DisplacementFieldPointer  m_BackwardDisplacementField;

  unsigned int              m_NumberOfLevels;
  unsigned int              m_CurrentLevel;
  std::vector<unsigned int> m_NumberOfIterations;

  /** Flag to indicate user stop registration request. */
  bool        m_StopRegistrationFlag;
  std::string m_DisplacementFieldOutputNamePrefix;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiResolutionICCDeformableRegistration.hxx"
#endif

#endif
