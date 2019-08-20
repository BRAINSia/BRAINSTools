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
#ifndef __TransformAdaptor_h
#define __TransformAdaptor_h

#include "CrossOverAffineSystem.h"

namespace itk
{
/** \class TransformAdaptor
 *
 * This component converts transform formats required for
 * input and output for use in registration.
 *
 * The preprocessing is activatived by method ExecuteInput().
 *
 * The postprocessing is activatived by method ExecuteOutput().
 *
 * Inputs:
 *    - pointer to original fixed image
 *    - pointer original moving image
 *
 * Outputs:
 *    - pointer to transform representing the pre-transform
 *    - pointer to transform representing the post-transform
 *
 * After registration, the overall transform is obtained by
 * composing pre-transform, the registration transform and
 * the post-transform.
 *
 */
template < typename TCoordinateType, unsigned int NDimensions, typename TInputImage >
class TransformAdaptor : public LightProcessObject
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( TransformAdaptor );

  /** Standard class type alias. */
  using Self = TransformAdaptor;
  using Superclass = Object;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods). */
  itkTypeMacro( TransformAdaptor, Object );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Input Image Type. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageSizeType = typename InputImageType::SizeType;

  /** Image dimension enumeration. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Type of the scalar representing coordinate and vector elements. */
  typedef TCoordinateType ScalarType;

  typedef vnl_matrix_fixed< TCoordinateType, NDimensions + 1, NDimensions + 1 > VnlTransformMatrixType44;

  /** Affine transform type. */
  using AffineTransformType = AffineTransform< TCoordinateType, Self::ImageDimension >;
  typedef typename AffineTransformType::Pointer          AffineTransformPointer;
  typedef typename AffineTransformType::MatrixType       MatrixType;
  typedef typename AffineTransformType::InputPointType   PointType;
  typedef typename AffineTransformType::OutputVectorType VectorType;
  typedef typename VectorType::ValueType                 ValueType;

  /** CrossOverAffineSystem type. */
  using CrossOverAffineSystemType = CrossOverAffineSystem< TCoordinateType, Self::ImageDimension >;
  typedef typename CrossOverAffineSystemType::Pointer CrossOverAffineSystemPointer;

  /** Set the input fixed image. */
  iplSetMacro( FixedImage, InputImagePointer );
  iplGetMacro( FixedImage, InputImagePointer );

  /** Set the input moving image. */
  iplSetMacro( MovingImage, InputImagePointer );
  iplGetMacro( MovingImage, InputImagePointer );

  /** Methods to execute the transform processing. */
  void
  ExecuteInput();

  void
  ExecuteOutput();

  /** Methods to define and perform Air16 AffineTransform conversion. */
  void
  EstablishCrossOverSystemForAir16( void );

  void
  EstablishCrossOverSystemForB2xfrm( void );

  void
  ConvertInputAffineToITKAffine( void );

  void
  ConvertITKAffineToOutputAffine( void );

  iplSetMacro( CenterMovingAffineTransform, AffineTransformPointer );
  iplGetMacro( CenterMovingAffineTransform, AffineTransformPointer );
  iplSetMacro( DeCenterMovingAffineTransform, AffineTransformPointer );
  iplGetMacro( DeCenterMovingAffineTransform, AffineTransformPointer );
  iplSetMacro( CenterFixedAffineTransform, AffineTransformPointer );
  iplGetMacro( CenterFixedAffineTransform, AffineTransformPointer );
  iplSetMacro( DeCenterFixedAffineTransform, AffineTransformPointer );
  iplGetMacro( DeCenterFixedAffineTransform, AffineTransformPointer );

  iplGetMacro( InputAffineTransformFilename, std::string );
  iplSetMacro( InputAffineTransformFilename, std::string );
  iplGetMacro( OutputAffineTransformFilename, std::string );
  iplSetMacro( OutputAffineTransformFilename, std::string );

  iplSetMacro( InputAffineTransform, AffineTransformPointer );
  iplGetMacro( InputAffineTransform, AffineTransformPointer );
  iplSetMacro( ITKAffineTransform, AffineTransformPointer );
  iplGetMacro( ITKAffineTransform, AffineTransformPointer );
  iplSetMacro( OutputAffineTransform, AffineTransformPointer );
  iplGetMacro( OutputAffineTransform, AffineTransformPointer );

  iplSetMacro( CrossOverAffineSystem, CrossOverAffineSystemPointer );
  iplGetMacro( CrossOverAffineSystem, CrossOverAffineSystemPointer );

protected:
  TransformAdaptor();
  ~TransformAdaptor() {}

private:
  InputImagePointer m_FixedImage;
  InputImagePointer m_MovingImage;

  AffineTransformPointer m_CenterFixedAffineTransform;
  AffineTransformPointer m_CenterMovingAffineTransform;

  AffineTransformPointer m_DeCenterFixedAffineTransform;
  AffineTransformPointer m_DeCenterMovingAffineTransform;

  std::string m_InputAffineTransformFilename;
  std::string m_OutputAffineTransformFilename;

  AffineTransformPointer m_ITKAffineTransform;
  AffineTransformPointer m_InputAffineTransform;
  AffineTransformPointer m_OutputAffineTransform;

  CrossOverAffineSystemPointer m_CrossOverAffineSystem;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "TransformAdaptor.hxx"
#endif

#endif
