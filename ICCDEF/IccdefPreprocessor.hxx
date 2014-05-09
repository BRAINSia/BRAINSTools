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
#ifndef _IccdefPreprocessor_txx
#define _IccdefPreprocessor_txx

#include "IccdefPreprocessor.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkBOBFFilter.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIO.h"
#include "itkMedianImageFilter.h"

#include "itkTransformToDisplacementFieldSource.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransform.h"
#include "itkWarpImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
IccdefPreprocessor<TInputImage, TOutputImage>
::IccdefPreprocessor()
{
  m_UseHistogramMatching = 0;
  m_NumberOfHistogramLevels = 256;
  m_NumberOfMatchPoints = 2;

  m_FixedImageMinimum =  NumericTraits<InputPixelType>::NonpositiveMin();
  m_MovingImageMinimum = NumericTraits<InputPixelType>::NonpositiveMin();

  m_FixedBinaryVolume = "none";
  m_MovingBinaryVolume = "none";
  for( unsigned i = 0; i < TInputImage::ImageDimension; ++i )
    {
    m_MedianFilterSize[i] = 0;
    }

  m_DefaultPixelValue = NumericTraits<PixelType>::One;
  m_OutDebug = false;
}

template <typename TInputImage, typename TOutputImage>
void
IccdefPreprocessor<TInputImage, TOutputImage>
::Execute()
{
  this->m_InputFixedImage = itkUtil::ReadImage<InputImageType>(this->GetTheFixedImageFilename() );
  this->m_InputMovingImage = itkUtil::ReadImage<InputImageType>(this->GetTheMovingImageFilename() );

  // resacle the image intensity to [0,255.0]
#if 1
  typedef itk::MinimumMaximumImageFilter<InputImageType> MinimumMaximumImageType;
  typename MinimumMaximumImageType::Pointer fixedMMImage = MinimumMaximumImageType::New();
  fixedMMImage->SetInput(this->m_InputFixedImage);
  fixedMMImage->Update();

  typename MinimumMaximumImageType::Pointer movingMMImage = MinimumMaximumImageType::New();
  movingMMImage->SetInput(this->m_InputMovingImage);
  movingMMImage->Update();

  // itk::NumericTraits< typename OutputImageType::PixelType >::max()
  // InputPixelType fixedMaximum = fixedMMImage->GetMaximum();
  typedef itk::RescaleIntensityImageFilter<InputImageType, InputImageType> RescaleIntensityImageType;
  if( fixedMMImage->GetMaximum() > (InputPixelType) itk::NumericTraits<unsigned char>::max()
      || fixedMMImage->GetMinimum() < (InputPixelType) itk::NumericTraits<unsigned char>::min() )
    {
    typename RescaleIntensityImageType::Pointer rescale_FixedImage = RescaleIntensityImageType::New();
    rescale_FixedImage->SetInput(this->m_InputFixedImage);
    rescale_FixedImage->SetOutputMaximum( (InputPixelType)255.0);
    rescale_FixedImage->SetOutputMinimum( (InputPixelType)0.0);
    rescale_FixedImage->Update();
    this->m_InputFixedImage = rescale_FixedImage->GetOutput();
    // rescale_FixedImage->DisconnectPipeline();
    if( this->GetOutDebug() )
      {
      std::cout << "Rescale the input fixed image from [" << fixedMMImage->GetMinimum() << ", "
                << fixedMMImage->GetMaximum() << "] to: [0,255]" << std::endl;
      }
    }
  if( movingMMImage->GetMaximum() > (InputPixelType) itk::NumericTraits<unsigned char>::max()
      || movingMMImage->GetMinimum() < (InputPixelType) itk::NumericTraits<unsigned char>::min() )
    {
    typename RescaleIntensityImageType::Pointer rescale_MovingImage = RescaleIntensityImageType::New();
    rescale_MovingImage->SetInput(this->m_InputMovingImage);
    rescale_MovingImage->SetOutputMaximum( (InputPixelType)255.0);
    rescale_MovingImage->SetOutputMinimum( (InputPixelType)0.0);
    rescale_MovingImage->Update();
    this->m_InputMovingImage = rescale_MovingImage->GetOutput();
    // rescale_MovingImage->DisconnectPipeline();
    if( this->GetOutDebug() )
      {
      std::cout << "Rescale the input moving image from [" << movingMMImage->GetMinimum() << ", "
                << movingMMImage->GetMaximum() << "] to: [0,255]" << std::endl;
      }
    }
#endif

  // Make BOBF Images if specified
  if( this->m_FixedBinaryVolume != std::string("none") )
    {
    if( this->GetOutDebug() )
      {
      std::cout << "Making BOBF \n";
      std::cout << "PRE Fixed Origin" << m_InputFixedImage->GetOrigin()
                << std::endl;
      }
    m_InputFixedImage = this->MakeBOBFImage( m_InputFixedImage,
                                             m_FixedBinaryVolume );
    if( this->GetOutDebug() )
      {
      std::cout << "Fixed Origin" << m_InputFixedImage->GetOrigin()
                << std::endl;
      std::cout << "PRE Moving Origin" << m_InputMovingImage->GetOrigin()
                << std::endl;
      }
    m_InputMovingImage = this->MakeBOBFImage( m_InputMovingImage,
                                              m_MovingBinaryVolume);
    if( this->GetOutDebug() )
      {
      std::cout << "Moving Origin" << m_InputMovingImage->GetOrigin()
                << std::endl;
      std::cout << "Writing Brain Only Background Filled Moving image"
                << std::endl;
      itkUtil::WriteImage<TInputImage>(m_InputMovingImage, "BOBF_Moving.nii.gz");
      itkUtil::WriteImage<TInputImage>(m_InputFixedImage, "BOBF_Fixed.nii.gz");
      }
    }

  if( m_InitialTransformFilename != "" )
    {
    //
    // read in the initial ITKTransform
    //
    typedef itk::TransformFileReader TransformReaderType;
    typename TransformReaderType::Pointer affineReader =  TransformReaderType::New();
    typedef typename TransformReaderType::TransformType BaseTransformType;
    BaseTransformType * baseTransform(0);

    std::cout << "Read ITK transform from text file: " << m_InitialTransformFilename << std::endl;

    affineReader->SetFileName( m_InitialTransformFilename.c_str() );
    try
      {
      affineReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      exit(-1);
      }
    const typename itk::TransformFileReader::TransformListType *  transforms = affineReader->GetTransformList();
    if( transforms->size() == 0 )
      {
      std::cout << "Number of transforms = " << transforms->size() << std::endl;
      std::cout << "FATAL ERROR: Failed to read transform"
                << m_InitialTransformFilename << std::endl;
      exit(-1);
      }

    //  #######Now use TransformToDisplacementFieldSource
    typedef itk::TransformToDisplacementFieldSource<TDisplacementField, double> DisplacementFieldGeneratorType;
    typedef typename DisplacementFieldGeneratorType::TransformType              TransformType; // Only a templated base
                                                                                               // class.

    std::cout << "Number of transforms = " << transforms->size() << std::endl;

    itk::TransformFileReader::TransformListType::const_iterator it
      = transforms->begin();

    TransformType* trsf = NULL;
    if( transforms->size() == 1 ) // There is no bulk transform.
      {
      // baseTransform = transforms->front();
      baseTransform = ( *it ).GetPointer();
      trsf = dynamic_cast<TransformType *>(baseTransform);

      const std::string firstNameOfClass = ( *it )->GetNameOfClass();
      std::cout << "FIRST (and only) NameOfClass = " << firstNameOfClass << std::endl;
      }
    else // Pick up what we presume was the bulk transform followed by a BSpline.
      {
      typename AffineTransformType::Pointer
      BulkTransform = static_cast<AffineTransformType *>( ( *it ).GetPointer() );
      const std::string firstNameOfClass = ( *it )->GetNameOfClass();
      std::cout << "First (Bulk) NameOfClass = " << firstNameOfClass << std::endl;
      ++it;
      const std::string secondNameOfClass = ( *it )->GetNameOfClass();
      std::cout << "SECOND NameOfClass = " << secondNameOfClass << std::endl;
      baseTransform = ( *it ).GetPointer();
      trsf = dynamic_cast<TransformType *>(baseTransform);
      if( secondNameOfClass == "BSplineTransform" )
        {
        typename BSplineTransformType::Pointer
        ITKTransform = static_cast<BSplineTransformType *>(baseTransform);
        ITKTransform->SetBulkTransform( BulkTransform );
        }
      else
        {
        std::cout << "Number of transforms in transform file " << m_InitialTransformFilename
                  <<
          " > 1, but ValidationInputParser (for BRAINSDemonWarp) only handles a transform list when the second transform is in fact a BSplineTransform, not this "
                  << secondNameOfClass << "." << std::endl;
        exit(-1);
        }
      }

    typename DisplacementFieldGeneratorType::Pointer defGenerator = DisplacementFieldGeneratorType::New();
    defGenerator->SetOutputSize( m_InputFixedImage->GetLargestPossibleRegion().GetSize() );
    defGenerator->SetOutputSpacing( m_InputFixedImage->GetSpacing() );
    defGenerator->SetOutputOrigin( m_InputFixedImage->GetOrigin() );
    defGenerator->SetOutputIndex( m_InputFixedImage->GetLargestPossibleRegion().GetIndex() );
    defGenerator->SetOutputDirection( m_InputFixedImage->GetDirection() );
    defGenerator->SetTransform(trsf);
    try
      {
      defGenerator->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      exit(-1);
      }
    m_InitialDisplacementField = defGenerator->GetOutput();
    m_InitialDisplacementField = itkUtil::OrientImage<TDisplacementField>(m_InitialDisplacementField,
                                                                          m_InputFixedImage->GetDirection() );

    // Warp the moving image with the initial deformation field
    typedef itk::WarpImageFilter<TInputImage, TInputImage, TDisplacementField> WarpImageType;
    typename WarpImageType::Pointer warper = WarpImageType::New();
    warper->SetInput(m_InputMovingImage);
    warper->SetOutputOrigin(m_InputFixedImage->GetOrigin() );
    warper->SetOutputSpacing(m_InputFixedImage->GetSpacing() );
    warper->SetOutputDirection(m_InputFixedImage->GetDirection() );
    warper->SetDisplacementField(m_InitialDisplacementField);
    warper->GetOutput()->SetRequestedRegion(m_InitialDisplacementField->GetRequestedRegion() );
    warper->Update();
    m_InputMovingImage = warper->GetOutput();
    }

  // Orientation the image
  m_InputFixedImage = itkUtil::OrientImage<InputImageType>(m_InputFixedImage,
                                                           itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
  m_InputMovingImage = itkUtil::OrientImage<InputImageType>(m_InputMovingImage,
                                                            itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);

//  std::cout<<"fixed direction:"<<this->m_InputFixedImage->GetDirection()<<std::endl;
//  std::cout<<"moving direction:"<<this->m_InputMovingImage->GetDirection()<<std::endl;

  if( m_MedianFilterSize[0] > 0  ||  m_MedianFilterSize[1] > 0 ||  m_MedianFilterSize[2] > 0 )
    {
    std::cout << "Using Medican fileter..." << std::endl;
    typedef typename itk::MedianImageFilter<TInputImage,
                                            TInputImage> MedianImageFilterType;
    typename MedianImageFilterType::Pointer medianFilter
      = MedianImageFilterType::New();
    medianFilter->SetRadius(m_MedianFilterSize);
    medianFilter->SetInput(m_InputFixedImage);
    medianFilter->Update();
    m_InputFixedImage = medianFilter->GetOutput();
    //
    // reinitialize
    medianFilter = MedianImageFilterType::New();
    medianFilter->SetRadius(m_MedianFilterSize);
    medianFilter->SetInput(m_InputMovingImage);
    medianFilter->Update();
    m_InputMovingImage = medianFilter->GetOutput();
    }

  // Create UnNormalized...Images
    {
    this->m_UnNormalizedFixedImage
      = itkUtil::PreserveCast<TInputImage, TOutputImage>(
          this->m_InputFixedImage);
    }
    {
    m_UnNormalizedMovingImage
      = itkUtil::PreserveCast<TInputImage, TOutputImage>(
          this->m_InputMovingImage);
    }

  m_OutputMovingImage = itkUtil::CopyImage<TOutputImage>(
      m_UnNormalizedMovingImage);

  if( this->GetUseHistogramMatching() )
    {
    typedef HistogramMatchingImageFilter<OutputImageType,
                                         OutputImageType> HistogramMatchingFilterType;
    typename HistogramMatchingFilterType::Pointer histogramfilter
      = HistogramMatchingFilterType::New();
    if( this->GetOutDebug() )
      {
      std::cout << "Performing Histogram Matching \n";
      }
    if( ( vcl_numeric_limits<typename OutputImageType::PixelType>::max()
          - vcl_numeric_limits<typename OutputImageType::PixelType>::min() ) <
        m_NumberOfHistogramLevels )
      {
      std::cout << "The intensity of range is less than Histogram levels!!"
                << std::endl;
      }
    histogramfilter->SetInput( m_UnNormalizedMovingImage  );
    histogramfilter->SetReferenceImage( m_UnNormalizedFixedImage);

    histogramfilter->SetNumberOfHistogramLevels( m_NumberOfHistogramLevels );
    histogramfilter->SetNumberOfMatchPoints( m_NumberOfMatchPoints );
    histogramfilter->ThresholdAtMeanIntensityOn();
    histogramfilter->Update();
    m_OutputMovingImage  = histogramfilter->GetOutput();
    }

  m_OutputFixedImage = itkUtil::CopyImage<TOutputImage>(
      m_UnNormalizedFixedImage);

//  itkUtil::Normalize<TOutputImage>(m_OutputFixedImage);
//  itkUtil::Normalize<TOutputImage>(m_OutputMovingImage);

  if( this->GetOutDebug() )
    {
    std::cout << "Writing Histogram equalized image" << std::endl;
    itkUtil::WriteImage<TOutputImage>(m_OutputFixedImage,
                                      "HistogramReferenceFixedImage.nii.gz");
    itkUtil::WriteImage<TOutputImage>(m_OutputMovingImage,
                                      "HistogramModifiedMovingImage.nii.gz");
    }

  m_InputMovingImage = NULL;
  m_InputFixedImage = NULL;
}

/*This function takes in a brain image and a whole brain mask and strips the
  skull of the image. It uses the BOBF filter to perform the skull stripping.*/

template <typename TInputImage, typename TOutputImage>
typename IccdefPreprocessor<TInputImage,
                            TOutputImage>::InputImagePointer IccdefPreprocessor<TInputImage,
                                                                                TOutputImage>
::MakeBOBFImage( InputImagePointer input, std::string MaskName )
{
  InputImagePointer Mask = itkUtil::ReadImage<InputImageType>(MaskName);

  typedef BOBFFilter<InputImageType, InputImageType> BOBFFilterType;
  typename BOBFFilterType::Pointer BOBFfilter = BOBFFilterType::New();
  if( this->GetOutDebug() )
    {
    std::cout
      <<
      "Making Brain only Background filled image with the following parameters. "
      << std::endl;
    std::cout << "Background fill Value:  " << m_DefaultPixelValue << std::endl;
    }

  BOBFfilter->SetReplaceValue( (InputPixelType)m_DefaultPixelValue );
  BOBFfilter->SetInputImage( input );
  BOBFfilter->SetInputMask( Mask );
  try
    {
    BOBFfilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    exit(-1);
    }

  InputImagePointer output = BOBFfilter->GetOutput();
  return output;
}
}   // namespace itk

#endif // _IccdefPreprocessor_txx
