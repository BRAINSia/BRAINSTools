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
#ifndef __CrossOverAffineSystem_h
#define __CrossOverAffineSystem_h

#include <itkLightProcessObject.h>
#include <itkMatrix.h>
#include <itkVector.h>
#include <itkAffineTransform.h>
#include <itkVersorTransform.h>
#include <itkVersorRigid3DTransform.h>
#include "itkScaleVersor3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkMacro.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_svd.h>
#include <vcl_compiler.h>
#include <iostream>
#include <cmath>
#include <cstdio>

/*
 *  The "Coordinate System" problem is closely related to generation
 *  of reliable choices for an affine transform.  AffineTransforms need
 *  to be generated in terms of a graphical vector space to map from,
 *  and a possibly different graphical vector space to map to.
 *
 *
 */

template < typename TCoordinateType, unsigned int NDimensions = 3 >
class CrossOverAffineSystem : public itk::LightProcessObject
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( CrossOverAffineSystem );

  /** Standard class type alias. */
  using Self = CrossOverAffineSystem;
  using Superclass = itk::LightProcessObject;
  using Pointer = itk::SmartPointer< Self >;
  using ConstPointer = itk::SmartPointer< const Self >;

  /** New method for creating an object using a factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CrossOverAffineSystem, itk::LightProcessObject );

  /** Dimension of the domain space. */
  static constexpr unsigned int SpaceDimension = NDimensions;
  static constexpr unsigned int AffineDimension = NDimensions + 1;

  /** Type of the scalar representing coordinate and vector elements. */
  using ScalarType = TCoordinateType;

  typedef vnl_matrix_fixed< TCoordinateType, NDimensions + 1, NDimensions + 1 > VnlTransformMatrixType44;
  // typedef vnl_matrix_fixed<TCoordinateType, NDimensions+1, NDimensions+1>
  //  VnlTransformMatrixType33;

  /** Affine conversion type for this class */
  using AffineTransformType = itk::AffineTransform< TCoordinateType, NDimensions >;
  using AffineTransformPointer = typename AffineTransformType::Pointer;
  using MatrixType = typename AffineTransformType::MatrixType;
  using PointType = typename AffineTransformType::InputPointType;
  using VectorType = typename AffineTransformType::OutputVectorType;
  using ValueType = typename VectorType::ValueType;

  /** Quaternion conversion types for this class */
  using VersorTransformType = itk::VersorTransform< TCoordinateType >;
  using VersorTransformPointer = typename VersorTransformType::Pointer;
  using VersorParametersType = typename VersorTransformType::ParametersType;

  using VersorRigid3DTransformType = itk::VersorRigid3DTransform< TCoordinateType >;
  typedef typename VersorRigid3DTransformType::Pointer        VersorRigid3DTransformPointer;
  typedef typename VersorRigid3DTransformType::ParametersType VersorRigid3DParametersType;

  using ScaleVersor3DTransformType = itk::ScaleVersor3DTransform< TCoordinateType >;
  typedef typename ScaleVersor3DTransformType::Pointer        ScaleVersor3DTransformPointer;
  typedef typename ScaleVersor3DTransformType::ParametersType ScaleVersor3DParametersType;

  using ScaleSkewVersor3DTransformType = itk::ScaleSkewVersor3DTransform< TCoordinateType >;
  typedef typename ScaleSkewVersor3DTransformType::Pointer        ScaleSkewVersor3DTransformPointer;
  typedef typename ScaleSkewVersor3DTransformType::ParametersType ScaleSkewVersor3DParametersType;

  /** Get the four coordinated AffineTransform conversions. */
  itkGetMacro( InhaleEncodeConversion, AffineTransformPointer );
  itkGetMacro( InhaleDecodeConversion, AffineTransformPointer );
  itkGetMacro( ExhaleEncodeConversion, AffineTransformPointer );
  itkGetMacro( ExhaleDecodeConversion, AffineTransformPointer );

  /** Generate the four coordinated AffineTransform conversions. */
  void
  EncloseInScaling( const VectorType & EncodeScale, const VectorType & DecodeScale );

  void
  EncloseInTranslation( const VectorType & EncodeShift, const VectorType & DecodeShift );

  void
  EncloseInCentering( const PointType & EncodeCenter, const PointType & DecodeCenter );

  void
  EncloseInAffineTransforms( AffineTransformPointer EncodeAffineTransform,
                             AffineTransformPointer DecodeAffineTransform );

protected:
  /** Set the four coordinated AffineTransform conversions. */
  itkSetMacro( InhaleEncodeConversion, AffineTransformPointer );
  itkSetMacro( InhaleDecodeConversion, AffineTransformPointer );
  itkSetMacro( ExhaleEncodeConversion, AffineTransformPointer );
  itkSetMacro( ExhaleDecodeConversion, AffineTransformPointer );

  CrossOverAffineSystem();
  ~CrossOverAffineSystem() override;

  mutable AffineTransformPointer m_InhaleEncodeConversion;
  mutable AffineTransformPointer m_InhaleDecodeConversion;
  mutable AffineTransformPointer m_ExhaleEncodeConversion;
  mutable AffineTransformPointer m_ExhaleDecodeConversion;
};

#ifndef ITK_MANUAL_INSTANTIATION
#  include "CrossOverAffineSystem.hxx"
#endif

#endif
