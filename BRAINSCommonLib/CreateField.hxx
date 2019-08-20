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
#ifndef __CreateField_hxx
#define __CreateField_hxx
#include "CreateField.h"
#include "itkIOCommon.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkIO.h"
#include <cerrno>

namespace itk
{
template < typename TImage, typename T2Image >
CreateField< TImage, T2Image >::CreateField()
  : m_Image1Filename( "" )
  , m_Image2Filename( "" )
  , m_ParameterFilename( "" )
  , m_ImageOne( NULL )
  , m_ImageTwo( NULL )
  , m_NumberOfHistogramLevels( 1024 )
  , m_NumberOfMatchPoints( 7 )
  , m_NumberOfLevels( 1 )
  , m_FixedImage( NULL )
  , m_MovingImage( NULL )
  , m_NumberOfHistogramLevels( 256 )
  , m_NumberOfMatchPoints( 1 )
  , m_FixedImageMinimum( 0 )
  , m_MovingImageMinimum( 0 )
  , m_DisplacementField( NULL )
{
  m_Image1ShrinkFactors.Fill( 1 ), m_Image2ShrinkFactors.Fill( 1 )

                                     m_NumberOfIterations = IterationsArrayType( 1 );
  m_NumberOfIterations.Fill( 10 );

  m_FixedImagePyramid = FixedImagePyramidType::New();
  m_FixedImagePyramid->UseShrinkImageFilterOff();
  m_MovingImagePyramid = MovingImagePyramidType::New();
  m_MovingImagePyramid->UseShrinkImageFilterOff();
  m_Registration = RegistrationType::New();
  m_Registration->InPlaceOn();

  m_Registration->SetFixedImagePyramid( m_FixedImagePyramid );
  m_Registration->SetMovingImagePyramid( m_MovingImagePyramid );

  using CommandType = SimpleMemberCommand< Self >;
  typename CommandType::Pointer command = CommandType::New();
  command->SetCallbackFunction( this, &Self::StartNewLevel );

  m_Tag = m_Registration->AddObserver( IterationEvent(), command );
}

template < typename TImage, typename T2Image >
CreateField< TImage, T2Image >::~CreateField()
{
  m_Registration->RemoveObserver( m_Tag );
}

template < typename TImage, typename T2Image >
void
CreateField< TImage, T2Image >::Execute()
{
  try
  {
    std::cout << "  Reading Input Images and Parameters" << std::endl;
    m_ImageOne = itkUtil::ReadImage< TImage >( m_Image1Filename );
    m_ImageTwo = itkUtil::ReadImage< TImage >( m_Image2Filename );

    FILE * paramFile = fopen( m_ParameterFilename.c_str(), "r" );
    if ( !paramFile )
    {
      itkExceptionMacro( << "  Could not open parameter file. " );
    }

    unsigned int uNumber;

    if ( fscanf( paramFile, "%d", &uNumber ) != 1 )
    {
      itkExceptionMacro( << "  Could not find the number of histogram levels." );
    }
    m_NumberOfMatchPoints = uNumber;

    if ( fscanf( paramFile, "%d", &uNumber ) != 1 )
    {
      itkExceptionMacro( << "Could not find the number of match points." );
    }
    m_NumberOfHistogramLevels = uNumber;

    if ( fscanf( paramFile, "%d", &uNumber ) != 1 )
    {
      itkExceptionMacro( << "Could not find the number of levels." );
    }
    m_NumberOfLevels = uNumber;

    {
      itk::Array< unsigned int > temp( m_NumberOfLevels );
      temp.Fill( 0 );
      m_NumberOfIterations = temp;
    }
    for ( unsigned int j = 0; j < m_NumberOfLevels; ++j )
    {
      if ( fscanf( paramFile, "%d", &uNumber ) != 1 )
      {
        itkExceptionMacro( << "Could not find number of iterations per level. " );
      }
      m_NumberOfIterations[j] = uNumber;
    }
    for ( unsigned int j = 0; j < ImageDimension; ++j )
    {
      if ( fscanf( paramFile, "%d", &uNumber ) != 1 )
      {
        itkExceptionMacro( << "Could not find atlas starting shrink factor. " );
      }
      m_Image1ShrinkFactors[j] = uNumber;
    }
    for ( unsigned int j = 0; j < ImageDimension; ++j )
    {
      if ( fscanf( paramFile, "%d", &uNumber ) != 1 )
      {
        itkExceptionMacro( << "  Could not find subject starting shrink factor. " );
      }
      m_Image2ShrinkFactors[j] = uNumber;
    }
  }
  catch ( itk::ExceptionObject & err )
  {
    std::cout << "  Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
  }
  catch ( ... )
  {
    std::cout << "  Error occurred during input parsing." << std::endl;
    throw;
  }

  std::cout << "  Preprocessing the images" << std::endl;
  try
  {
    this->NormalizeImage( m_ImageTwo, m_FixedImage, m_FixedImageMinimum );
    this->NormalizeImage( m_ImageOne, m_MovingImage, m_MovingImageMinimum );

    using FilterType = HistogramMatchingImageFilter< OutputImageType, OutputImageType >;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput( m_MovingImage );
    filter->SetReferenceImage( m_FixedImage );
    filter->SetNumberOfHistogramLevels( m_NumberOfHistogramLevels );
    filter->SetNumberOfMatchPoints( m_NumberOfMatchPoints );
    filter->ThresholdAtMeanIntensityOn();
    filter->Update();

    m_MovingImage = filter->GetOutput();
  }
  catch ( itk::ExceptionObject & err )
  {
    std::cout << "  Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
  }
  catch ( ... )
  {
    std::cout << "  Error occured during preprocessing." << std::endl;
    throw;
  }

  std::cout << "  Registering the images" << std::endl;

  try
  {
    m_FixedImage->SetMetaDataDictionary( m_ImageTwo->GetMetaDataDictionary() );
    m_MovingImage->SetMetaDataDictionary( m_ImageOne->GetMetaDataDictionary() );
    if ( ( m_FixedImage->GetDirection() != m_MovingImage->GetDirection() )
         // INFO:  Remove dependance on RIP from
         ( itk::SpatialOrientationAdapter().FromDirectionCosines( m_FixedImage->GetDirection() ) !=
           itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP ) )
    {
      std::cout << "Image Directions are not the same or are not in RIP orientation " << std::endl
                << m_FixedImage->GetDirection() << "=============" << std::endl
                << m_MovingImage->GetDirection() << std::endl;
    }
    m_ImageOne = NULL;
    m_ImageTwo = NULL;

    m_FixedImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
    m_FixedImagePyramid->SetStartingShrinkFactors( m_Image2ShrinkFactors.GetDataPointer() );

    m_MovingImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
    m_MovingImagePyramid->SetStartingShrinkFactors( m_Image1ShrinkFactors.GetDataPointer() );

    m_Registration->SetFixedImage( m_FixedImage );
    m_Registration->SetMovingImage( m_MovingImage );
    m_Registration->SetNumberOfLevels( m_NumberOfLevels );
    m_Registration->SetNumberOfIterations( m_NumberOfIterations.data_block() );
    try
    {
      m_Registration->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
      std::cout << "  Caught an exception: " << std::endl;
      std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
    }
    catch ( ... )
    {
      std::cout << errno << " " << __FILE__ << " " << __LINE__ << std::endl;
      std::cout << "  Caught a non-ITK exception " << __FILE__ << " " << __LINE__ << std::endl;
    }

    try
    {
      m_Registration->ReleaseDataFlagOn();
      m_DisplacementField = m_Registration->GetOutput();
      // m_DisplacementField->DisconnectPipeline();
    }
    catch ( itk::ExceptionObject & err )
    {
      std::cout << "  Caught an exception: " << std::endl;
      std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
    }
    catch ( ... )
    {
      std::cout << "  Caught a non-ITK exception " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
    }
  }
  catch ( itk::ExceptionObject & err )
  {
    std::cout << "  Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
  }
  catch ( ... )
  {
    std::cout << "  Error occured during registration" << std::endl;
    throw;
  }
}

template < typename TImage, typename T2Image >
void
CreateField< TImage, T2Image >::ReleaseDataFlagOn()
{
  m_FixedImage = NULL;
  m_MovingImage = NULL;
}

template < typename TImage, typename T2Image >
void
CreateField< TImage, T2Image >::NormalizeImage( InputImageType * input, OutputImagePointer & output,
                                                InputPixelType & min )
{
  using MinMaxFilterType = MinimumMaximumImageFilter< InputImageType >;
  typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();

  minMaxFilter->SetInput( input );
  minMaxFilter->Update();

  min = minMaxFilter->GetMinimum();
  double shift = -1.0 * static_cast< double >( min );
  double scale = static_cast< double >( minMaxFilter->GetMaximum() );
  scale += shift;
  scale = 1.0 / scale;

  using FilterType = ShiftScaleImageFilter< InputImageType, OutputImageType >;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  filter->SetShift( shift );
  filter->SetScale( scale );
  filter->Update();

  output = filter->GetOutput();
}

template < typename OutputImageType, typename InputImageType >
void
FFCreateNewImageFromTemplate( typename OutputImageType::Pointer &      PointerToOutputImage,
                              const typename InputImageType::Pointer & PreInitializedImage )
{
  PointerToOutputImage = OutputImageType::New();
  PointerToOutputImage->SetRegions( PreInitializedImage->GetLargestPossibleRegion() );
  PointerToOutputImage->CopyInformation( PreInitializedImage );
  PointerToOutputImage->Allocate();
  PointerToOutputImage->FillBuffer( 0 );
  CHECK_CORONAL( PointerToOutputImage->GetDirection() );
}

template < typename TImage, typename T2Image >
void
CreateField< TImage, T2Image >::StartNewLevel()
{
  std::cout << "  Starting level " << m_Registration->GetCurrentLevel() << std::endl;
}
} // namespace itk
#endif
