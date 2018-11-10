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
#ifndef _IccdefRegistrator_txx
#define _IccdefRegistrator_txx

#include <sstream>
#include <string>
#include "IccdefRegistrator.h"
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

#include "itkImageFileWriter.h"

#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"

#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itksys/SystemTools.hxx"

namespace itk
{
/*This function writes the displacement fields of the Deformation.*/
template <
  typename TRealImage,
  typename TOutputImage,
  typename TFieldValue>
void IccdefRegistrator<TRealImage, TOutputImage,
                       TFieldValue>::WriteDisplacementComponents(TDisplacementField * df, std::string prefix)
{
//  using OutputImageType = itk::Image<float, 3>;
  this->m_DefaultPixelValue = NumericTraits<PixelType>::OneValue();

  // we use the vector index selection filter to break the deformation field
  // into x,y,z components.
  using ComponentImageType = itk::Image<FieldValueType,
                     3>;
  using ComponentFilterType = itk::VectorIndexSelectionCastImageFilter<TDisplacementField,
                                                   ComponentImageType>;

  std::string CurrentComponentFilename;
  try
    {
    char ext[3][14] = { "/xdisp.nii.gz", "/ydisp.nii.gz", "/zdisp.nii.gz"};

    typename ComponentFilterType::Pointer myComponentFilter
      = ComponentFilterType::New();
    myComponentFilter->SetInput(df);
    for( unsigned int extiter = 0; extiter < 3; extiter++ )
      {
      CurrentComponentFilename = prefix + ext[extiter];
      if( this->GetOutDebug() )
        {
        std::cout << "Writing Transform Image: "
                  << CurrentComponentFilename << std::endl;
        }

      myComponentFilter->SetIndex(extiter);

      typename ComponentImageType::Pointer DisplacementComponentImagePtr
        = myComponentFilter->GetOutput();

      itkUtil::WriteImage<ComponentImageType>(DisplacementComponentImagePtr,
                                              CurrentComponentFilename);
      }
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "exception in file Displacement File Writer("
              << CurrentComponentFilename << ")" << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
}

/*Constructor to initialize the parameters.*/
template <
  typename TRealImage,
  typename TOutputImage,
  typename TFieldValue>
IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::IccdefRegistrator()
{
  // Images need to be set from the outside
  this->m_FixedImage = nullptr;
  this->m_MovingImage = nullptr;
  this->m_DisplacementField = nullptr;
  this->m_BackwardDisplacementField = nullptr;

  // Set up internal registrator with default components
  this->m_FixedImagePyramid = FixedImagePyramidType::New();
  this->m_FixedImagePyramid->UseShrinkImageFilterOff();
  this->m_MovingImagePyramid = MovingImagePyramidType::New();
  this->m_MovingImagePyramid->UseShrinkImageFilterOff();
  this->m_Registration = RegistrationType::New();
  // m_Registration->InPlaceOn(); -- not in ITKv4

//  this->m_Registration->SetFixedImagePyramid (this->m_FixedImagePyramid);
//  this->m_Registration->SetMovingImagePyramid (this->m_MovingImagePyramid);

  this->m_DefaultPixelValue =  NumericTraits<typename RealImageType::PixelType>::ZeroValue();
  // Setup an registration observer
  using CommandType = SimpleMemberCommand<Self>;
  typename CommandType::Pointer command = CommandType::New();
  command->SetCallbackFunction(this, &Self::StartNewLevel);

  this->m_Tag = this->m_Registration->AddObserver(IterationEvent(), command);

  using FieldInterpolatorType = VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<
      TDisplacementField, double>;

  typename FieldInterpolatorType::Pointer VectorInterpolator12
    = FieldInterpolatorType::New();

  typename FieldInterpolatorType::Pointer VectorInterpolator21
    = FieldInterpolatorType::New();

  this->m_Registration->GetModifiableFieldExpander12()->SetInterpolator(VectorInterpolator12);
  this->m_Registration->GetModifiableFieldExpander21()->SetInterpolator(VectorInterpolator21);

  // Default parameters
  this->m_NumberOfLevels = 1;

  this->m_NumberOfIterations = UnsignedIntArray(1);
  this->m_NumberOfIterations.Fill(10);
  this->m_OutputPrefix = "none";

  this->m_ForwardDisplacementFieldOutputName = "none";
  this->m_BackwardDisplacementFieldOutputName = "none";
  this->m_InitialFixedDisplacementFieldFilename = "none";
  this->m_InitialMovingDisplacementFieldFilename = "none";
  this->m_OutputJacobianImage = false;
  this->m_OutputDisplacement = false;
  this->m_OutputDisplacementField = false;

  this->m_UseHistogramMatching = false;
  this->m_OutDebug = false;

  this->m_InitialDisplacementField = nullptr;
  this->m_ForwardDir = "forward";
  this->m_BackwardDir = "backward";
}

template <
  typename TRealImage,
  typename TOutputImage,
  typename TFieldValue>
IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::~IccdefRegistrator()
{
  if( this->m_Tag )
    {
    this->m_Registration->RemoveObserver(this->m_Tag);
    }
}

/*Perform the registration of preprocessed images.*/
template <
  typename TRealImage,
  typename TOutputImage,
  typename TFieldValue>
void IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::Execute()
{

  this->m_Registration->SetFixedImage(this->m_FixedImage);
  this->m_Registration->SetMovingImage(this->m_MovingImage);

  this->m_Registration->SetNumberOfLevels(this->m_NumberOfLevels);
  this->m_Registration->SetNumberOfIterations( this->m_NumberOfIterations.
                                               data_block() );

  this->m_ForwardDir = "forward";
  this->m_BackwardDir = "backward";

  if( this->m_OutputDisplacementField )
    {
    this->m_Registration->SetDisplacementFieldOutputNamePrefix(this->m_OutputPrefix);
    }

  // Setup the initial deformation field
  if( this->m_InitialFixedDisplacementFieldFilename != std::string("none")
      && this->m_InitialFixedDisplacementFieldFilename != std::string("") )
    {
    using FieldReaderType = itk::ImageFileReader<TDisplacementField>;
    typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( this->m_InitialFixedDisplacementFieldFilename.c_str() );
    try
      {
      fieldReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Caught an ITK exception: " << std::endl;
      throw;
      }

    this->m_Registration->SetInitialFixedDisplacementField(fieldReader->GetOutput() );
    }

  if( this->m_InitialMovingDisplacementFieldFilename != std::string("none")
      && this->m_InitialMovingDisplacementFieldFilename != std::string("") )
    {
    using FieldReaderType = itk::ImageFileReader<TDisplacementField>;
    typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( this->m_InitialMovingDisplacementFieldFilename.c_str() );
    try
      {
      fieldReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Caught an ITK exception: " << std::endl;
      throw;
      }

    this->m_Registration->SetInitialMovingDisplacementField(fieldReader->GetOutput() );
    }

  // Perform the registration.
  try
    {
    this->m_Registration->Update();
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

  try
    {
    this->m_DisplacementField = this->m_Registration->GetOutput(0);
    this->m_BackwardDisplacementField = this->m_Registration->GetOutput(1);
    if( this->m_DisplacementField->GetDirection() != this->m_FixedImage->GetDirection() )
      {
      std::cout << "ERROR Directions don't match\n"
                << this->m_DisplacementField->GetDirection()
                << "\n"
                << this->m_FixedImage->GetDirection()
                << std::endl;
      exit(-1);
      }
    if( this->m_Tag )
      {
      this->m_Registration->RemoveObserver(this->m_Tag);
      this->m_Tag = 0;
      }
    this->m_Registration = nullptr;
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

#if 1
  // Write the output deformation fields if specified by the user.
  if( this->m_ForwardDisplacementFieldOutputName != std::string("none")
      && this->m_ForwardDisplacementFieldOutputName != std::string("") )
    {
    itkUtil::WriteImage<TDisplacementField>(this->m_DisplacementField,
                                            this->m_ForwardDisplacementFieldOutputName);
    if( this->GetOutDebug() )
      {
      std::cout << "---Forward deformation field has been written "
                << this->m_ForwardDisplacementFieldOutputName << "--" << std::endl;
      }
    }
  if( this->m_BackwardDisplacementFieldOutputName != std::string("none")
      && this->m_BackwardDisplacementFieldOutputName != std::string("") )
    {
    itkUtil::WriteImage<TDisplacementField>(this->m_BackwardDisplacementField,
                                            this->m_ForwardDisplacementFieldOutputName);
    if( this->GetOutDebug() )
      {
      std::cout << "---Backward deformation field has been written "
                << this->m_BackwardDisplacementFieldOutputName << "--" << std::endl;
      }
    }
#endif

  //  Write out the displacement fields specified by the user.
  if( this->GetOutputDisplacement() )
    {
    this->WriteDisplacementComponents(this->m_DisplacementField, this->m_ForwardDir);
    this->WriteDisplacementComponents(this->m_BackwardDisplacementField, this->m_BackwardDir);
    }

  if( this->GetOutputJacobianImage() )
    {
    using JacobianFilterType = itk::DisplacementFieldJacobianDeterminantFilter<
        TDisplacementField, FieldValueType, TRealImage>;

    typename JacobianFilterType::Pointer localForwardJacobianFilter = JacobianFilterType::New();
    localForwardJacobianFilter->SetUseImageSpacing( true );
    localForwardJacobianFilter->SetInput(this->m_DisplacementField);
    localForwardJacobianFilter->UpdateLargestPossibleRegion();

    typename TRealImage::Pointer fjPtr = localForwardJacobianFilter->GetOutput();
    itkUtil::WriteImage<TRealImage>(fjPtr, this->m_ForwardDir + "/" + "Jacobian_forward.nii.gz");

    typename JacobianFilterType::Pointer localBackwardJacobianFilter = JacobianFilterType::New();
    localBackwardJacobianFilter->SetUseImageSpacing( true );
    localBackwardJacobianFilter->SetInput(this->m_BackwardDisplacementField);
    localBackwardJacobianFilter->UpdateLargestPossibleRegion();
    typename TRealImage::Pointer bjPtr = localBackwardJacobianFilter->GetOutput();
    itkUtil::WriteImage<TRealImage>(bjPtr, this->m_BackwardDir + "/" + "Jacobian_backward.nii.gz");
    }

  if( this->m_OutputPrefix != std::string("none") )
    {
    /*Warp the image with the generated deformation field.*/
    typename RealImageType::Pointer DeformedMovingImagePtr(nullptr);
    typename RealImageType::Pointer DeformedFixedImagePtr(nullptr);

    using WarperType = WarpImageFilter<RealImageType, RealImageType,
                            TDisplacementField>;
    typename WarperType::Pointer warper = WarperType::New();
    if( this->GetUseHistogramMatching() == true )
      {
      warper->SetInput( this->m_MovingImage );
      }
    else
      {
      warper->SetInput( this->m_UnNormalizedMovingImage );
      }

    warper->SetOutputSpacing( this->m_FixedImage->GetSpacing() );
    warper->SetOutputOrigin( this->m_FixedImage->GetOrigin() );
    warper->SetOutputDirection( this->m_FixedImage->GetDirection() );
    warper->SetDisplacementField( this->m_DisplacementField);
    warper->Update();
    DeformedMovingImagePtr = warper->GetOutput();

    typename WarperType::Pointer warperb = WarperType::New();
    warperb->SetInput(this->m_FixedImage);
    warperb->SetOutputSpacing( this->m_MovingImage->GetSpacing() );
    warperb->SetOutputOrigin( this->m_MovingImage->GetOrigin() );
    warperb->SetOutputDirection( this->m_MovingImage->GetDirection() );
    warperb->SetDisplacementField(this->m_BackwardDisplacementField);
    warperb->Update();
    DeformedFixedImagePtr = warperb->GetOutput();

    if( this->GetOutDebug() )
      {
      std::cout << "-----Direction of output warped image\n"
                << DeformedMovingImagePtr->GetDirection()
                << "\n-----Direction of deformation field\n"
                << this->m_DisplacementField->GetDirection() << std::endl;
      }
    /*Write the output image.*/
    if( this->m_OutputPrefix != std::string("none") )
      {
      typename OutputImageType::Pointer CastImageSptr
        = itkUtil::PreserveCast<RealImageType, OutputImageType>( DeformedMovingImagePtr);

      const std::string deformedMovingFilename(this->m_ForwardDir + "/" + "deformedMovingImage.nii.gz" );
      itkUtil::WriteImage<OutputImageType>(CastImageSptr, deformedMovingFilename);

      if( this->GetOutDebug() )
        {
        std::cout << "---Deformed Moving Image has been written: " << deformedMovingFilename << std::endl;
        }

      typename OutputImageType::Pointer CastImageSptr1
        = itkUtil::PreserveCast<RealImageType, OutputImageType>(DeformedFixedImagePtr);
      const std::string deformedFixedFilename(this->m_BackwardDir + "/" + "deformedFixedImage.nii.gz");
      itkUtil::WriteImage<OutputImageType>(CastImageSptr1, deformedFixedFilename);

      if( this->GetOutDebug() )
        {
        std::cout << "---Deformed Fixed Image has been written: " << deformedFixedFilename << std::endl;
        }
      }
    }
}

// Print out the present registration level.
template <
  typename TRealImage,
  typename TOutputImage,
  typename TFieldValue>
void IccdefRegistrator<TRealImage, TOutputImage, TFieldValue>::StartNewLevel()
{
  if( this->GetOutDebug() )
    {
    std::cout << "--- Starting level " << this->m_Registration->GetCurrentLevel()
              << std::endl;
    }
}
} // namespace itk
#endif
