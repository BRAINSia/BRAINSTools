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
#ifndef __ValidationInputParser_hxx
#define __ValidationInputParser_hxx

#include "ValidationInputParser.h"
#include "itkMetaDataObject.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIO.h"
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkScaleVersor3DTransform.h"
#include "ConvertToRigidAffine.h"

#include <fstream>
// The following was just copied out of iccdefWarpImage.cc.  Sorry.
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkSpatialOrientation.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include <itkIO.h>
#include <metaCommand.h>

#include "TransformToDisplacementField.h"
#include <itkDisplacementFieldJacobianDeterminantFilter.h>

#ifdef __USE_BRAINS2_INTEGRATION
#include "iccdefdeformation/Deformation.h"
#include "iccdefdeformation/Utils.h"
#include "b2Affine_rw.h"
#include "iccdefdeformation/HarmonicArrayIO.h"
#endif

namespace itk
{
template <typename TImage>
ValidationInputParser<TImage>
::ValidationInputParser() :
  m_TheMovingImageFilename(""),
  m_TheFixedImageFilename(""),
  m_ParameterFilename(""),
  m_TheMovingImage(NULL),
  m_TheFixedImage(NULL),
  m_ForceCoronalZeroOrigin(false),
  m_OutDebug(false),
  m_NumberOfHistogramLevels(1024),
  m_NumberOfMatchPoints(7),
  m_NumberOfLevels(1),
  m_NumberOfIterations(IterationsArrayType(1) )
{
  m_TheMovingImageShrinkFactors.Fill(1);
  m_TheFixedImageShrinkFactors.Fill(1);
  m_NumberOfIterations.Fill(10);
}

template <typename TImage>
void
ValidationInputParser<TImage>
::Execute()
{
  /*************************
    * Read in the images
    */
  if( this->m_ForceCoronalZeroOrigin == true )
    {
    itkGenericExceptionMacro(<< "---Forcing Brains2 Orientation not yet implemented");
    }
  else
    {
    m_TheFixedImage = itkUtil::ReadImage<TImage>(m_TheFixedImageFilename);
    m_TheMovingImage = itkUtil::ReadImage<TImage>(m_TheMovingImageFilename);
    }
  // TODO:  Need to ensure that the fixed and moving images have the same
  // orientations.

  // TODO:  Need to figure out how to read in the initial deformation field.
  // std::cerr << "About to check for deformation field file " <<
  // m_InitialDisplacementFieldFilename << std::endl;
  // std::cerr << "About to check for transform file " <<
  // m_InitialTransformFilename << std::endl;
  // std::cerr << "About to check for Coefficient file" <<
  // m_InitialCoefficientFilename << std::endl;
  if( m_InitialDisplacementFieldFilename != "" )
    {
    typedef   itk::ImageFileReader<TDisplacementField> FieldReaderType;
    typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( m_InitialDisplacementFieldFilename.c_str() );
    try
      {
      fieldReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Caught an ITK exception: "
                << err.GetDescription()
                << std::endl;
      throw;
      }
    if( this->GetOutDebug() )
      {
      std::cout << "\nReading Displacement fields.\n";
      }
    m_InitialDisplacementField = fieldReader->GetOutput();
    //  typename ImageType::DirectionType DisplacementOrientation;
    //  DisplacementOrientation=displacementField->GetDirection();
    }
  else if( m_InitialTransformFilename != "" )
    {
#if 1
    // TODO: Zhao,  please implement the following, and test that it works
    // within
    //      the regression test suite.

    //  #######Now use TransformToDisplacementFieldSource
#if (ITK_VERSION_MAJOR < 4)
    typedef itk::TransformToDeformationFieldSource<TDisplacementField, double> DisplacementFieldGeneratorType;
#else
    typedef itk::TransformToDisplacementFieldSource<TDisplacementField, double> DisplacementFieldGeneratorType;
#endif
    typedef typename DisplacementFieldGeneratorType::TransformType TransformType;
    // Only a templated base class.

    typename TransformType::Pointer trsf = ReadTransformFromDisk(m_InitialTransformFilename);

    typename DisplacementFieldGeneratorType::Pointer defGenerator = DisplacementFieldGeneratorType::New();
    defGenerator->SetOutputSpacing( this->GetTheFixedImage()->GetSpacing() );
    defGenerator->SetOutputOrigin( this->GetTheFixedImage()->GetOrigin() );
    defGenerator->SetOutputDirection( this->GetTheFixedImage()->GetDirection() );
    defGenerator->SetOutputSize( this->GetTheFixedImage()->GetLargestPossibleRegion().GetSize() );
    defGenerator->SetOutputIndex( this->GetTheFixedImage()->GetLargestPossibleRegion().GetIndex() );
    defGenerator->SetTransform(trsf);
    defGenerator->Update();
    m_InitialDisplacementField = defGenerator->GetOutput();
    // itkUtil::WriteImage<TDisplacementField>(m_InitialDisplacementField,
    // "initialDeformationfield.nii.gz");
#endif
    }
  else if( m_InitialCoefficientFilename != "" )
    {
#ifdef __USE_BRAINS2_INTEGRATION
    DisplacementFieldFFTType::Pointer mu;        // mu1, mu2, mu3;
    std::string                       CoeffNameInput(
      m_InitialCoefficientFilename.c_str() );

      {
      if( this->GetOutDebug() )
        {
        std::cout << "Reading: " << CoeffNameInput << std::endl;
        }
      HarmonicReadAll3D(mu, CoeffNameInput);
      }
    if( this->GetOutDebug() )
      {
      std::cout << "\nCreating Displacement fields from Coefficient files\n";
      }
    m_InitialDisplacementField = CreateITKDisplacementFieldFromCoeffs(mu);
#else
    itkGenericExceptionMacro(<< "ERROR:  InitialCoefficientFilename not supported yet");
#endif
    }

  // Print out the parameters.
  if( this->GetOutDebug() )
    {
    std::cout << "NumberOfHistogramLevels : " << m_NumberOfHistogramLevels
              << std::endl;
    std::cout << "NumberOfMatchPoints : " << m_NumberOfMatchPoints << std::endl;
    std::cout << "NumberOfLevels : " << m_NumberOfLevels << std::endl;
    std::cout << "NumberOfIterations : " << m_NumberOfIterations << std::endl;
    std::cout << "TheMovingImageShrinkFactors : "
              << m_TheMovingImageShrinkFactors << std::endl;
    std::cout << "TheFixedImageShrinkFactors : "
              << m_TheFixedImageShrinkFactors << std::endl;
    }
}
} // namespace itk

#endif
