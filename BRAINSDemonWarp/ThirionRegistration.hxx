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
#ifndef __ThirionRegistration_hxx
#define __ThirionRegistration_hxx

#include "ThirionRegistration.h"

namespace itk
{
template <typename TImage, typename TRealImage, typename TOutputImage>
ThirionRegistration<TImage, TRealImage, TOutputImage>
::ThirionRegistration() :
  m_TheMovingImageFilename(""),
  m_TheFixedImageFilename(""),
  m_DisplacementBaseName("none"),
  m_WarpedImageName("none"),
  m_ForceCoronalZeroOrigin(false),
  m_UseHistogramMatching(false),
  m_OutNormalized("OFF"),
  m_OutputFilename(""),
  m_CheckerBoardFilename("none"),
  m_DisplacementFieldOutputName("none"),
  m_AppendOutputFile(true),
  m_BOBFTargetMask("none"),
  m_BOBFTemplateMask("none"),
  m_Lower(NumericTraits<PixelType>::NonpositiveMin() ),
  m_Upper(NumericTraits<PixelType>::max() ),
  m_DefaultPixelValue(NumericTraits<PixelType>::ZeroValue()),
  m_NumberOfHistogramLevels(256),
  m_NumberOfMatchPoints(2),
  m_NumberOfLevels(4)
{
  m_CheckerBoardPattern.Fill(4);
  m_Radius.Fill(1);
  m_NumberOfIterations = IterationsArrayType(m_NumberOfLevels);
  m_NumberOfIterations[0] = 2000;
  m_NumberOfIterations[1] = 500;
  m_NumberOfIterations[2] = 250;
  m_NumberOfIterations[3] = 100;
  for( unsigned i = 0; i < ImageType::ImageDimension; ++i )
    {
    m_TheMovingImageShrinkFactors[i] = 4;   // 16;
    m_TheFixedImageShrinkFactors[i] = 4;    // 16;
    m_Seed[i] = 0;
    m_MedianFilterSize[i] = 0;
    }
}

/*This method initializes the input parser which reads in the moving image,
  * fixed image and parameter file.*/

template <typename TImage, typename TRealImage, typename TOutputImage>
void
ThirionRegistration<TImage, TRealImage, TOutputImage>
::InitializeParser()
{
  this->m_Parser->SetTheMovingImageFilename(
    this->m_TheMovingImageFilename.c_str() );

  this->m_Parser->SetTheFixedImageFilename( this->m_TheFixedImageFilename.c_str() );
  this->m_Parser->SetForceCoronalZeroOrigin( this->GetForceCoronalZeroOrigin() );

  this->m_Parser->SetInitialDisplacementFieldFilename(
    this->m_InitialDisplacementFieldFilename.c_str() );
  this->m_Parser->SetInitialCoefficientFilename(
    this->m_InitialCoefficientFilename.c_str() );
  this->m_Parser->SetInitialTransformFilename(
    this->m_InitialTransformFilename.c_str() );
  //            this->m_Parser->SetParameterFilename(
  // this->m_ParameterFilename.c_str() );
  this->m_Parser->SetNumberOfHistogramLevels( this->GetNumberOfHistogramLevels() );
  this->m_Parser->SetNumberOfMatchPoints( this->GetNumberOfMatchPoints() );
  this->m_Parser->SetNumberOfLevels( this->GetNumberOfLevels() );
  this->m_Parser->SetTheMovingImageShrinkFactors(
    this->GetTheMovingImageShrinkFactors() );
  this->m_Parser->SetTheFixedImageShrinkFactors(
    this->GetTheFixedImageShrinkFactors() );
  this->m_Parser->SetNumberOfIterations( this->GetNumberOfIterations() );
  this->m_Parser->SetOutDebug( this->GetOutDebug() );
}

/*This method initializes the preprocessor which processes the moving and fixed
  * images before registration. The image files which are read in using the
  *    parser are given to the preprocessor.*/

template <typename TImage, typename TRealImage, typename TOutputImage>
void
ThirionRegistration<TImage, TRealImage, TOutputImage>
::InitializePreprocessor()
{
  this->m_Preprocessor->SetInputFixedImage( this->m_Parser->GetTheFixedImage() );
  this->m_Preprocessor->SetInputMovingImage( this->m_Parser->GetTheMovingImage() );
  this->m_Preprocessor->SetInitialDisplacementField( this->m_Parser->GetInitialDisplacementField() );
  this->m_Preprocessor->SetUseHistogramMatching( this->GetUseHistogramMatching() );
  this->m_Preprocessor->SetNumberOfHistogramLevels( this->m_Parser->GetNumberOfHistogramLevels() );
  this->m_Preprocessor->SetNumberOfMatchPoints( this->m_Parser->GetNumberOfMatchPoints() );
  this->m_Preprocessor->SetBOBFTargetMask( this->GetBOBFTargetMask() );
  this->m_Preprocessor->SetBOBFTemplateMask( this->GetBOBFTemplateMask() );
  this->m_Preprocessor->SetLower( this->GetLower() );
  this->m_Preprocessor->SetUpper( this->GetUpper() );
  this->m_Preprocessor->SetRadius( this->GetRadius() );
  this->m_Preprocessor->SetDefaultPixelValue( this->GetDefaultPixelValue() );
  this->m_Preprocessor->SetSeed( this->GetSeed() );
  this->m_Preprocessor->SetOutDebug( this->GetOutDebug() );
  this->m_Preprocessor->SetMedianFilterSize( this->GetMedianFilterSize() );
  this->m_Preprocessor->SetInitialDisplacementField( this->m_Parser->GetInitialDisplacementField() );
}

/*This method initializes the registration process. The preprocessed output
  * files are passed to the registrator.*/

template <typename TImage, typename TRealImage, typename TOutputImage>
void
ThirionRegistration<TImage, TRealImage, TOutputImage>
::InitializeRegistrator()
{
  this->m_Registrator->SetDisplacementBaseName( this->GetDisplacementBaseName() );
  this->m_Registrator->SetWarpedImageName( this->GetWarpedImageName() );
  this->m_Registrator->SetCheckerBoardFilename( this->GetCheckerBoardFilename() );
  this->m_Registrator->SetDisplacementFieldOutputName( this->GetDisplacementFieldOutputName() );
  this->m_Registrator->SetCheckerBoardPattern( this->GetCheckerBoardPattern() );
  this->m_Registrator->SetFixedImage( this->m_Preprocessor->GetOutputFixedImage() );
  this->m_Registrator->SetMovingImage( this->m_Preprocessor->GetOutputMovingImage() );
  this->m_Registrator->SetUnNormalizedMovingImage( this->m_Preprocessor->GetUnNormalizedMovingImage() );
  this->m_Registrator->SetUnNormalizedFixedImage( this->m_Preprocessor->GetUnNormalizedFixedImage() );
  this->m_Registrator->SetInitialDisplacementField( this->m_Parser->GetInitialDisplacementField() );
  this->m_Registrator->SetDefaultPixelValue( this->m_Preprocessor->GetDefaultPixelValue() );
  this->m_Registrator->SetUseHistogramMatching( this->GetUseHistogramMatching() );
  this->m_Registrator->SetNumberOfLevels( this->m_Parser->GetNumberOfLevels() );
  this->m_Registrator->SetNumberOfIterations( this->m_Parser->GetNumberOfIterations() );
  this->m_Registrator->SetInterpolationMode( this->GetInterpolationMode() );

  this->m_Registrator->SetFixedImageShrinkFactors( this->m_Parser->GetTheFixedImageShrinkFactors() );
  this->m_Registrator->SetMovingImageShrinkFactors( this->m_Parser->GetTheMovingImageShrinkFactors() );

  this->m_Registrator->SetOutNormalized( this->GetOutNormalized() );
  this->m_Registrator->SetOutDebug( this->GetOutDebug() );
  this->m_Registrator->SetDisplacementFieldOutputName( this->m_DisplacementFieldOutputName);
  this->m_Registrator->SetFixedLandmarkFilename(this->m_FixedLandmarkFilename);
  this->m_Registrator->SetMovingLandmarkFilename(this->m_MovingLandmarkFilename);
}
}   // namespace itk

#endif
