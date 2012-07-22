#ifndef _ICCDEFWarp_txx
#define _ICCDEFWarp_txx

#include "ICCDEFWarp.h"

namespace itk
{
template <typename TImage, typename TRealImage, typename TOutputImage>
ICCDEFWarp<TImage, TRealImage, TOutputImage>
::ICCDEFWarp()
{
  m_TheMovingImageFilename = "";
  m_TheFixedImageFilename = "";

  m_FixedBinaryVolume = "none";
  m_MovingBinaryVolume = "none";

  // m_ParameterFilename = "";
//  m_OutputFilename = "";
  m_AppendOutputFile = true;
//  m_WarpedImageName = "none";

  m_ForwardDisplacementFieldOutputName = "none";
  m_BackwardDisplacementFieldOutputName = "none";
//  m_DisplacementBaseName = "none";
  m_OutputPrefix = "";
  m_OutputDisplacement = false;
  m_OutputJacobianImage = false;
  m_OutputDisplacementField = false;

  m_DefaultPixelValue = NumericTraits<PixelType>::Zero;

  m_UseHistogramMatching = false;
  //    m_OutDebug = false;
  m_NumberOfHistogramLevels = 256;
  m_NumberOfMatchPoints = 2;
  m_NumberOfLevels = 4;   // if that fixes it, I'm going to do something else
  m_NumberOfIterations = IterationsArrayType(m_NumberOfLevels);
  m_NumberOfIterations[0] = 2000;
  m_NumberOfIterations[1] = 500;
  m_NumberOfIterations[2] = 250;
  m_NumberOfIterations[3] = 100;
  for( unsigned i = 0; i < ImageType::ImageDimension; i++ )
    {
    // m_Seed[i] = 0;
    m_MedianFilterSize[i] = 0;
    }
}

/*This method initializes the preprocessor which processes the moving and fixed
  images before registration. The image files which are read in using the parser
  are given to the preprocessor.*/

template <typename TImage, typename TRealImage, typename TOutputImage>
void
ICCDEFWarp<TImage, TRealImage, TOutputImage>
::InitializePreprocessor()
{
  this->m_Preprocessor->SetTheFixedImageFilename( this->GetTheFixedImageFilename() );
  this->m_Preprocessor->SetTheMovingImageFilename( this->GetTheMovingImageFilename() );
  this->m_Preprocessor->SetUseHistogramMatching( this->GetUseHistogramMatching() );
  this->m_Preprocessor->SetNumberOfHistogramLevels(this->GetNumberOfHistogramLevels() );
  this->m_Preprocessor->SetNumberOfMatchPoints(this->GetNumberOfMatchPoints() );
  this->m_Preprocessor->SetDefaultPixelValue( this->GetDefaultPixelValue() );
  this->m_Preprocessor->SetOutDebug( this->GetOutDebug() );
  this->m_Preprocessor->SetMedianFilterSize( this->GetMedianFilterSize() );
  this->m_Preprocessor->SetFixedBinaryVolume( this->GetFixedBinaryVolume() );
  this->m_Preprocessor->SetMovingBinaryVolume( this->GetMovingBinaryVolume() );
  this->m_Preprocessor->SetInitialTransformFilename(this->GetInitialTransformFilename() );
}

/*This method initializes the registration process. The preprocessed output
  files are passed to the registrator.*/

template <typename TImage, typename TRealImage, typename TOutputImage>
void
ICCDEFWarp<TImage, TRealImage, TOutputImage>
::InitializeRegistrator()
{
  this->m_Registrator->SetOutputDisplacement( this->GetOutputDisplacement() );
  this->m_Registrator->SetOutputPrefix( this->GetOutputPrefix() );
  this->m_Registrator->SetForwardDisplacementFieldOutputName(
    this->GetForwardDisplacementFieldOutputName() );
  this->m_Registrator->SetBackwardDisplacementFieldOutputName(
    this->GetBackwardDisplacementFieldOutputName() );
  this->m_Registrator->SetFixedImage( this->m_Preprocessor->GetOutputFixedImage() );
  this->m_Registrator->SetMovingImage(
    this->m_Preprocessor->GetOutputMovingImage() );
  this->m_Registrator->SetUnNormalizedMovingImage(
    this->m_Preprocessor->GetUnNormalizedMovingImage() );
  this->m_Registrator->SetUnNormalizedFixedImage(
    this->m_Preprocessor->GetUnNormalizedFixedImage() );
  this->m_Registrator->SetDefaultPixelValue(
    this->m_Preprocessor->GetDefaultPixelValue() );
  this->m_Registrator->SetUseHistogramMatching( this->GetUseHistogramMatching() );
  this->m_Registrator->SetNumberOfLevels( this->GetNumberOfLevels() );
  this->m_Registrator->SetNumberOfIterations(this->GetNumberOfIterations() );
  this->m_Registrator->SetOutDebug(this->GetOutDebug() );
  this->m_Registrator->SetInitialFixedDisplacementFieldFilename(this->m_InitialFixedDisplacementFieldFilename);
  this->m_Registrator->SetInitialMovingDisplacementFieldFilename(this->m_InitialMovingDisplacementFieldFilename);
  this->m_Registrator->SetOutputJacobianImage(this->GetOutputJacobianImage() );
  this->m_Registrator->SetOutputDisplacementField(this->GetOutputDisplacementField() );
}
}   // namespace itk

#endif
