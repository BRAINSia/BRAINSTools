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
#ifndef __DemonsRegistrator_hxx
#define __DemonsRegistrator_hxx

#include <string>
#include <sstream>
#include "DemonsRegistrator.h"
#include "itkCommand.h"
#include "ApplyField.h"
#include "itkStatisticsImageFilter.h"
#include "itkMetaImageIO.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "vector"
#include "itkCheckerBoardImageFilter.h"
#include "itkIO.h"

#include "itkMultiResolutionPDEDeformableRegistration.h"
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkWarpImageFilter.h"
#include "itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h"
#include "GenericTransformImage.h"
#include "DebugImageWrite.h"

namespace itk
{
/*This function writes the displacement fields of the Displacement.*/
template <class TRealImage, class TOutputImage, class TFieldValue>
void DemonsRegistrator<TRealImage, TOutputImage,
                       TFieldValue>::WriteDisplacementComponents()
{
  m_DefaultPixelValue = NumericTraits<PixelType>::OneValue();

  // we use the vector index selection filter to break the deformation field
  // into x,y,z components.
  typedef itk::Image<FieldValueType,
                     3>                  ComponentImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<TDisplacementField,
                                                   ComponentImageType> ComponentFilterType;

  std::string CurrentComponentFilename;
  try
    {
    char ext[3][14] = { "_xdisp.nii.gz", "_ydisp.nii.gz", "_zdisp.nii.gz" };

    typename ComponentFilterType::Pointer myComponentFilter =
      ComponentFilterType::New();
    myComponentFilter->SetInput(m_DisplacementField);
    for( unsigned int extiter = 0; extiter < 3; ++extiter )
      {
      CurrentComponentFilename = m_DisplacementBaseName + ext[extiter];
      if( this->GetOutDebug() )
        {
        std::cout << "Writing Transform Image: "
                  << CurrentComponentFilename << std::endl;
        }

      myComponentFilter->SetIndex(extiter);

      typename ComponentImageType::Pointer DisplacementComponentImagePtr =
        myComponentFilter->GetOutput();

      itkUtil::WriteImage<ComponentImageType>(DisplacementComponentImagePtr,
                                              CurrentComponentFilename);
      }
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    throw;
    }
}

/*Constructor to initialize the parameters.*/
template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue>
DemonsRegistrator<TRealImage, TOutputImage, TFieldValue>::DemonsRegistrator() :
  m_InitialDisplacementField(NULL),
  m_FixedImage(NULL),
  m_MovingImage(NULL),
  m_FixedImagePyramid(FixedImagePyramidType::New() ),
  m_MovingImagePyramid(MovingImagePyramidType::New() ),
  m_Registration(RegistrationType::New() ),
  m_DefaultPixelValue( NumericTraits<typename RealImageType::PixelType>::Zero),
  m_NumberOfLevels(1),
  m_NumberOfIterations(UnsignedIntArray(1) ),
  m_DisplacementField(NULL),
  m_DisplacementBaseName("none"),
  m_WarpedImageName("none"),
  m_CheckerBoardFilename("none"),
  m_DisplacementFieldOutputName("none"),
  m_OutNormalized("OFF"),
  m_OutDebug(false),
  m_UseHistogramMatching(false),
  m_InterpolationMode("Linear")
{
  //TODO: Needed for ITKv4 registration m_Registration->InPlaceOn();
  // Set up internal registrator with default components
  m_FixedImagePyramid->UseShrinkImageFilterOff();
  m_MovingImagePyramid->UseShrinkImageFilterOff();

  m_FixedImageShrinkFactors.Fill(1);
  m_MovingImageShrinkFactors.Fill(1);

  m_NumberOfIterations.Fill(10);
  m_CheckerBoardPattern.Fill(4);
}

template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue>
DemonsRegistrator<TRealImage, TOutputImage, TFieldValue>::~DemonsRegistrator()
{
  if( m_Tag )
    {
    m_Registration->RemoveObserver(m_Tag);
    }
}

/*Perform the registration of preprocessed images.*/
template <
  typename TRealImage,
  class TOutputImage,
  class TFieldValue>
void DemonsRegistrator<TRealImage, TOutputImage, TFieldValue>::Execute()
{
  if( this->m_FixedLandmarkFilename != ""
      && this->m_MovingLandmarkFilename != "" )
    {
    itkGenericExceptionMacro(<< "Registering Landmarks as an initializer is not yet implemented");
    }
#if 1
  // Setup the image pyramids
  m_FixedImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
  m_FixedImagePyramid->SetStartingShrinkFactors(m_FixedImageShrinkFactors.GetDataPointer() );

  m_MovingImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
  m_MovingImagePyramid->SetStartingShrinkFactors( m_MovingImageShrinkFactors.GetDataPointer() );
#endif
  // Setup the registrator

    {
    // Setup an registration observer
    typedef SimpleMemberCommand<Self> CommandType;
    typename CommandType::Pointer command = CommandType::New();
    command->SetCallbackFunction(this, &Self::StartNewLevel);

    m_Tag = m_Registration->AddObserver(IterationEvent(), command);

    typedef VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<
        TDisplacementField, double> FieldInterpolatorType;

    typename FieldInterpolatorType::Pointer VectorInterpolator =
      FieldInterpolatorType::New();

    m_Registration->GetModifiableFieldExpander()->SetInterpolator(VectorInterpolator);

    m_Registration->SetFixedImagePyramid(m_FixedImagePyramid);
    m_Registration->SetMovingImagePyramid(m_MovingImagePyramid);
    m_Registration->SetFixedImage(m_FixedImage);
    m_Registration->SetMovingImage(m_MovingImage);
    m_Registration->SetNumberOfLevels(m_NumberOfLevels);
    m_Registration->SetNumberOfIterations( m_NumberOfIterations.data_block() );

    DebugOutput(RealImageType, m_FixedImage);
    DebugOutput(RealImageType, m_MovingImage);

    m_MovingImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
    m_MovingImagePyramid->SetInput(m_MovingImage);
    m_MovingImagePyramid->UpdateLargestPossibleRegion();

    m_FixedImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
    m_FixedImagePyramid->SetInput(m_FixedImage);
    m_FixedImagePyramid->UpdateLargestPossibleRegion();
    for( unsigned int i = 0; i < m_NumberOfLevels; i++ )
      {
      DebugOutputN(RealImageType, this->m_FixedImagePyramid->GetOutput(i), i, m_FixedImagePyramid);
      DebugOutputN(RealImageType, this->m_MovingImagePyramid->GetOutput(i), i, m_FixedImagePyramid);
      }

    // Setup the initial deformation field
    if( this->m_InitialDisplacementField.IsNotNull() )
      {
      DebugOutput(TDisplacementField, this->m_InitialDisplacementField);
      m_Registration->SetInitialDisplacementField(this->m_InitialDisplacementField);
      }
    // Perform the registration.
    try
      {
      m_Registration->UpdateLargestPossibleRegion();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cout << "Caught an exception: " << std::endl;
      std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
      }
    catch( ... )
      {
      std::
      cout << "Caught a non-ITK exception " << __FILE__ << " " << __LINE__
           << std::endl;
      }
    if( this->GetOutDebug() )
      {
      std::cout
        << "Moving image shrink factors used in each level of MultiResolution Schedule\n"
        << m_MovingImagePyramid->GetSchedule() << std::endl;
      std::cout
        << "Fixed image shrink factors used in each level of MultiResolution Schedule\n"
        << m_FixedImagePyramid->GetSchedule() << std::endl;
      }
    try
      {
      m_DisplacementField = m_Registration->GetOutput();
      if( m_DisplacementField->GetDirection() != m_FixedImage->GetDirection() )
        {
        itkGenericExceptionMacro(<< "ERROR Directions don't match\n"
                                 << m_DisplacementField->GetDirection()
                                 << "\n"
                                 << m_FixedImage->GetDirection() );
        }
      if( m_Tag )
        {
        m_Registration->RemoveObserver(m_Tag);
        m_Tag = 0;
        }
      m_Registration = NULL;
      }
    catch( itk::ExceptionObject & err )
      {
      std::cout << "Caught an exception: " << std::endl;
      std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
      }
    catch( ... )
      {
      std::cout << "Caught a non-ITK exception " << __FILE__ << " " << __LINE__
                << std::endl;
      }
    }

  // Write the output deformation fields if specified by the user.
  if( this->m_DisplacementFieldOutputName != std::string("none")
      && this->m_DisplacementFieldOutputName != std::string("") )
    {
    itkUtil::WriteImage<TDisplacementField>(m_DisplacementField,
                                            this->m_DisplacementFieldOutputName);
    if( this->GetOutDebug() )
      {
      std::cout << "---Displacement field has been written "
                << this->m_DisplacementFieldOutputName << "--" << std::endl;
      }
    }

  //  Write out the displacement fields specified by the user.
  if( this->m_DisplacementBaseName != std::string("none") )
    {
    this->WriteDisplacementComponents();
    }

  if( this->m_WarpedImageName != std::string("none") || this->m_CheckerBoardFilename != std::string("none") )
    {
    /*Warp the image with the generated deformation field.*/
    typename RealImageType::Pointer DeformedMovingImagePtr(0);

      {
      typename RealImageType::Pointer sourceMovingImage = NULL;
      if( this->GetUseHistogramMatching() == true )
        {
        sourceMovingImage = m_MovingImage;
        }
      else
        {
        sourceMovingImage = m_UnNormalizedMovingImage;
        }
      DeformedMovingImagePtr = TransformWarp<RealImageType, RealImageType, TDisplacementField>(
          sourceMovingImage,
          m_FixedImage.GetPointer(),
          0,
          GetInterpolatorFromString<RealImageType>(this->m_InterpolationMode),
          m_DisplacementField);
      DebugOutput(RealImageType, sourceMovingImage);
      DebugOutput(RealImageType, m_FixedImage);
      DebugOutput(RealImageType, DeformedMovingImagePtr);
      DebugOutput(TDisplacementField, m_DisplacementField);
      }

    if( this->GetOutDebug() )
      {
      std::cout << "-----Direction of output warped image\n"
                << DeformedMovingImagePtr->GetDirection()
                << "\n-----Direction of deformation field\n"
                << this->m_DisplacementField->GetDirection() << std::endl;
      }

    /*Write the output image.*/

    if( this->m_WarpedImageName != std::string("") && this->m_WarpedImageName != std::string("none") )
      {
      typename TOutputImage::Pointer CastImageSptr =
        itkUtil::PreserveCast<RealImageType, TOutputImage>(
          DeformedMovingImagePtr);
      itkUtil::WriteImage<TOutputImage>(CastImageSptr,
                                        this->m_WarpedImageName);

      if( this->GetOutDebug() )
        {
        std::cout << "---Deformed Image has been written" <<  this->m_WarpedImageName << std::endl;
        }
      }
    /*Write the checkerboard image of the fixed image and the output image.*/
    if( this->m_CheckerBoardFilename != std::string("") && this->m_CheckerBoardFilename != std::string("none") )
      {
      typedef itk::CheckerBoardImageFilter<RealImageType> Checkerfilter;
      typename Checkerfilter::Pointer checker = Checkerfilter::New();
      if( this->GetUseHistogramMatching() == true )
        {
        checker->SetInput1(m_FixedImage);
        }
      else
        {
        checker->SetInput1(m_UnNormalizedFixedImage);
        }
      checker->SetInput2(DeformedMovingImagePtr);
      checker->SetCheckerPattern( this->GetCheckerBoardPattern() );
      try
        {
        checker->UpdateLargestPossibleRegion();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "Caught an ITK exception: " << std::endl;
        std::cout << err << " " << __FILE__ << " " << __LINE__ << std::
          endl;
        throw;
        }
      typename RealImageType::Pointer CheckerImagePtr = checker->GetOutput();
      itkUtil::WriteImage<RealImageType>(CheckerImagePtr,
                                         this->m_CheckerBoardFilename);
      if( this->GetOutDebug() )
        {
        std::cout << "---Checker Board Image has been written" << std::endl;
        }
      }
    }
}

// Print out the present registration level.
template <
  class TRealImage,
  class TOutputImage,
  class TFieldValue>
void DemonsRegistrator<TRealImage, TOutputImage, TFieldValue>::StartNewLevel()
{
  if( this->GetOutDebug() )
    {
    std::cout << "-- Starting level " << m_Registration->GetCurrentLevel() << std::endl;
    }
}
} // namespace itk
#endif
