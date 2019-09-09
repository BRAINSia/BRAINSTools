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
#ifndef _VValidationInputParser_hxx
#define _VValidationInputParser_hxx

#include "VValidationInputParser.h"
#include "itkMetaDataObject.h"
#include "itkNumericTraits.h"
#include "itkMath.h"
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIO.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkSpatialOrientation.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include <itkIO.h>
#include <metaCommand.h>
#include "itkImageRegionIterator.h"
#ifdef __USE_BRAINS2_INTEGRATION
#  include "TransformToDisplacementField.h"
#endif
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <fstream>

#ifdef __USE_BRAINS2_INTEGRATION
#  include "iccdefdeformation/Deformation.h"
#  include "iccdefdeformation/Utils.h"
#  include "b2Affine_rw.h"
#  include "iccdefdeformation/HarmonicArrayIO.h"
#endif

namespace itk
{
template < typename TImage >
VValidationInputParser< TImage >::VValidationInputParser()
{
  //      m_TheMovingImageFilename = "";
  //      m_TheFixedImageFilename = "";

  m_ParameterFilename = "";

  m_TheMovingImages.reserve( 10 );
  m_TheFixedImages.reserve( 10 );

  m_NumberOfHistogramLevels = 1024;
  m_NumberOfMatchPoints = 7;

  m_NumberOfLevels = 1;
  m_TheMovingImageShrinkFactors.Fill( 1 );
  m_TheFixedImageShrinkFactors.Fill( 1 );

  m_NumberOfIterations = IterationsArrayType( 1 );
  m_NumberOfIterations.Fill( 10 );

  m_OutDebug = false;
  m_ForceCoronalZeroOrigin = false;
}

template < typename TImage >
void
VValidationInputParser< TImage >::Execute()
{
  /*************************
   * Read in the images
   */
  if ( this->m_ForceCoronalZeroOrigin == true )
  {
    itkGenericExceptionMacro( << "---Forcing Brains2 Orientation not yet implemented" );
  }
  else
  {
    for ( unsigned int i = 0; i < m_TheFixedImageFilename.size(); ++i )
    {
      m_TheFixedImages.push_back( itkUtil::ReadImage< TImage >( m_TheFixedImageFilename[i] ) );
      m_TheMovingImages.push_back( itkUtil::ReadImage< TImage >( m_TheMovingImageFilename[i] ) );
    }
  }
  // HACK:  INFO:  Need to ensure that the fixed and moving images have the same
  // orientations.

  // INFO:  Need to figure out how to read in the initial deformation field.
  // std::cerr << "About to check for deformation field file " <<
  // m_InitialDisplacementFieldFilename << std::endl;
  // std::cerr << "About to check for transform file " <<
  // m_InitialTransformFilename << std::endl;
  // std::cerr << "About to check for Coefficient file" <<
  // m_InitialCoefficientFilename << std::endl;
  if ( m_InitialDisplacementFieldFilename != "" )
  {
    using FieldReaderType = itk::ImageFileReader< TDisplacementField >;
    typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( m_InitialDisplacementFieldFilename.c_str() );
    try
    {
      fieldReader->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
      std::cerr << "Caught an ITK exception: " << std::endl;
      throw err;
    }
    if ( this->GetOutDebug() )
    {
      std::cout << "\nReading Displacement fields.\n";
    }
    m_InitialDisplacementField = fieldReader->GetOutput();
    //  typename ImageType::DirectionType DisplacementOrientation;
    //  DisplacementOrientation=displacementField->GetDirection();
  }
#ifdef __USE_BRAINS2_INTEGRATION
  else if ( m_InitialTransformFilename != "" )
  {
    //  REFACTOR continuing here:
    //  At this point, setting an AffineTransform to generate the initial
    // deformation field
    //  reads and converts a brains2 .xfrm file.

    // read brains2 transform file
    using B2AffineTransformType = B2AffineTransform< ImageType >;
    using AffineTransformType = typename B2AffineTransformType::TransformType;
    B2AffineTransformType transform;
    transform.Read( m_InitialTransformFilename );
    typename AffineTransformType::Pointer inputAffineTransform = transform.GetAffineTransformPointer();
    if ( inputAffineTransform.IsNull() )
    {
      std::cerr << "Can't read transform file" << m_InitialTransformFilename << std::endl;
    }

    typename TDisplacementField::RegionType::SizeType size = transform.GetFixedImageSize();
    typename TDisplacementField::SpacingType          spacing = transform.GetFixedImageSpacing();

    // convert brains2 transform, which is in index, to ITK transform, which is
    // in mm
    using ITKAffineTransformType = itk::AffineTransform< double, 3 >;
    using VectorType = itk::Vector< double, 3 >;
    VectorType const fixedImageScaleReciprocal( Reciprocal< double, 3 >( transform.GetFixedImageSpacing() ) );

    VectorType const movingImageScale( transform.GetMovingImageSpacing() );

    using CrossOverAffineSystemType = CrossOverAffineSystem< double, 3 >;
    CrossOverAffineSystemType::Pointer crossOverAffineSystem = CrossOverAffineSystemType::New();
    crossOverAffineSystem->EncloseInScaling( fixedImageScaleReciprocal, movingImageScale );

    const bool                      ApplyUpstream = false;
    ITKAffineTransformType::Pointer InitialITKAffineTransform = ITKAffineTransformType::New();
    InitialITKAffineTransform->SetIdentity();
    InitialITKAffineTransform->Compose( crossOverAffineSystem->GetInhaleEncodeConversion(), ApplyUpstream );
    InitialITKAffineTransform->Compose( inputAffineTransform, ApplyUpstream );
    InitialITKAffineTransform->Compose( crossOverAffineSystem->GetInhaleDecodeConversion(), ApplyUpstream );
#  ifdef USE_TRANSFORM_INVERSE_FOR_INIT_FROM_AFFINE_TRANSFORM
    ITKAffineTransformType::Pointer InitialITKAffineTransformInverse = ITKAffineTransformType::New();
    InitialITKAffineTransform->GetInverse( InitialITKAffineTransformInverse );
    m_InitialDisplacementField = TransformToDisplacementField( m_TheFixedImage, InitialITKAffineTransformInverse );
#  else
    m_InitialDisplacementField = TransformToDisplacementField( m_TheFixedImage, InitialITKAffineTransform );
#  endif
  }
#endif
#ifdef __USE_BRAINS2_INTEGRATION
  else if ( m_InitialCoefficientFilename != "" )
  {
    DisplacementFieldFFTType::Pointer mu; // mu1, mu2, mu3;
    std::string                       CoeffNameInput( m_InitialCoefficientFilename.c_str() );

    {
      if ( this->GetOutDebug() )
      {
        std::cout << "Reading: " << CoeffNameInput << std::endl;
      }
      HarmonicReadAll3D( mu, CoeffNameInput );
    }
    if ( this->GetOutDebug() )
    {
      std::cout << "\nCreating Displacement fields from Coefficient files\n";
    }
    m_InitialDisplacementField = CreateITKDisplacementFieldFromCoeffs( mu );
  }
#endif

  // Print out the parameters.
  if ( this->GetOutDebug() )
  {
    std::cout << "NumberOfHistogramLevels : " << m_NumberOfHistogramLevels << std::endl;
    std::cout << "NumberOfMatchPoints : " << m_NumberOfMatchPoints << std::endl;
    std::cout << "NumberOfLevels : " << m_NumberOfLevels << std::endl;
    std::cout << "NumberOfIterations : " << m_NumberOfIterations << std::endl;
    std::cout << "TheMovingImageShrinkFactors : " << m_TheMovingImageShrinkFactors << std::endl;
    std::cout << "TheFixedImageShrinkFactors : " << m_TheFixedImageShrinkFactors << std::endl;
  }
}
} // namespace itk

#endif
