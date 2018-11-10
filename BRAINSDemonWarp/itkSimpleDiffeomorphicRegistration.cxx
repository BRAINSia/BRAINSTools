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
/*=========================================================================
 *
 *  Program:   BRAINSDemonsWarp
 *  Module:    $RCSfile: ITKHeader.h,v $
 *  Language:  C++
 *  Date:      $Date: 2006-04-25 23:20:16 $
 *  Version:   $Revision: 1.1 $
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/

#include "itkSimpleDiffeomorphicRegistration.h"
#include "itkVector.h"
#include "itkArray.h"

#include "itkDiffeomorphicDemonsRegistrationWithMaskFilter.h"
#include "itkPDEDeformableRegistrationFilter.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"
#include "itkImageMaskSpatialObject.h"

itkSimpleDiffeomorphicRegistration::itkSimpleDiffeomorphicRegistration()
{
  m_DemonsPreprocessor = DemonsPreprocessorType::New();
  m_DemonsRegistrator = DemonsRegistratorType::New();
  m_DisplacementField = NULL;
}

/*
  * void itkSimpleDiffeomorphicRegistration::SetDeformedImageName(std::string
  * name) {
  * m_DeformedImageName = name;
  * }
  */
/*
  * std::string itkSimpleDiffeomorphicRegistration::GetDeformedImageName(void) {
  * return m_DeformedImageName;
  * }
  */

void itkSimpleDiffeomorphicRegistration::InitializePreprocessor()
{
  m_DemonsPreprocessor->SetInputFixedImage(m_FixedImage);
  m_DemonsPreprocessor->SetInputMovingImage(m_MovingImage);
  m_DemonsPreprocessor->SetUseHistogramMatching(true);
  m_DemonsPreprocessor->SetNumberOfHistogramLevels(NumberOfHistogramLevels);
  m_DemonsPreprocessor->SetNumberOfMatchPoints(NumberOfMatchPoints);
}

void itkSimpleDiffeomorphicRegistration::Initialization()
{
  using RegistrationFilterType = itk::DiffeomorphicDemonsRegistrationWithMaskFilter<TRealImage,
                                                             TRealImage,
                                                             TDisplacementField>;
  RegistrationFilterType::Pointer filter = RegistrationFilterType::New();
  using BaseRegistrationFilterType = itk::PDEDeformableRegistrationFilter<TRealImage, TRealImage,
                                               TDisplacementField>;
  BaseRegistrationFilterType::Pointer actualfilter;
  // using TGradientType = RegistrationFilterType::GradientType;
  TRealImage::Pointer movingBinaryVolumeImage;
  TRealImage::Pointer fixedBinaryVolumeImage;
  constexpr double otsuPercentileThreshold = 0.01;
  constexpr int closingSize = 7;
  //   fixedBinaryVolumeImage = FindLargestForgroundFilledMask<TRealImage>(
  //     m_FixedImage,
  //     otsuPercentileThreshold,
  //     closingSize);
  //   movingBinaryVolumeImage = FindLargestForgroundFilledMask<TRealImage>(
  //     m_MovingImage,
  //     otsuPercentileThreshold,
  //     closingSize);
  using LFFMaskFilterType = itk::LargestForegroundFilledMaskImageFilter<TRealImage>;
  LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
  LFF->SetInput(m_FixedImage);
  LFF->SetOtsuPercentileThreshold(otsuPercentileThreshold);
  LFF->SetClosingSize(closingSize);
  LFF->UpdateLargestPossibleRegion();
  fixedBinaryVolumeImage = LFF->GetOutput();
  //  LFF = LFFMaskFilterType::New();
  LFF->SetInput(m_MovingImage);
  LFF->SetOtsuPercentileThreshold(otsuPercentileThreshold);
  LFF->SetClosingSize(closingSize);
  LFF->UpdateLargestPossibleRegion();
  movingBinaryVolumeImage = LFF->GetOutput();

  using MaskPixelType = unsigned char;
  using MaskImageType = itk::Image<MaskPixelType, DIM>;
  using CastImageFilter = itk::CastImageFilter<TRealImage, MaskImageType>;

  using ImageMaskType = itk::SpatialObject<DIM>;
  using ImageMaskPointer = ImageMaskType::Pointer;

  CastImageFilter::Pointer castFixedMaskImage = CastImageFilter::New();
  castFixedMaskImage->SetInput(fixedBinaryVolumeImage);
  castFixedMaskImage->Update();

  // convert mask image to mask
  using ImageMaskSpatialObjectType = itk::ImageMaskSpatialObject<DIM>;
  ImageMaskSpatialObjectType::Pointer fixedMask = ImageMaskSpatialObjectType::New();
  fixedMask->SetImage( castFixedMaskImage->GetOutput() );
  fixedMask->ComputeObjectToWorldTransform();

  CastImageFilter::Pointer castMovingMaskImage = CastImageFilter::New();
  castMovingMaskImage->SetInput(movingBinaryVolumeImage);
  castMovingMaskImage->Update();

  // convert mask image to mask

  ImageMaskSpatialObjectType::Pointer movingMask = ImageMaskSpatialObjectType::New();
  movingMask->SetImage( castMovingMaskImage->GetOutput() );
  movingMask->ComputeObjectToWorldTransform();

  filter->SetFixedImageMask( dynamic_cast<ImageMaskType *>( fixedMask.GetPointer() ) );
  filter->SetMovingImageMask( dynamic_cast<ImageMaskType *>( fixedMask.GetPointer() ) );

  filter->SetMaximumUpdateStepLength(MaxStepLength);
  //  filter->SetUseGradientType(static_cast<TGradientType> (0));
  filter->SmoothDisplacementFieldOn();
  filter->SetStandardDeviations(SmoothDisplacementFieldSigma);
  filter->SmoothUpdateFieldOff();
  actualfilter = filter;
  m_DemonsRegistrator->SetRegistrationFilter(actualfilter);

  using IterationsArrayType = itk::Array<unsigned int>;
  IterationsArrayType numberOfIterations;
  numberOfIterations.SetSize(NumberOfLevels);
  numberOfIterations.SetElement(0, NumberOfIteration0);
  numberOfIterations.SetElement(1, NumberOfIteration1);
  numberOfIterations.SetElement(2, NumberOfIteration2);
  numberOfIterations.SetElement(3, NumberOfIteration3);
  numberOfIterations.SetElement(4, NumberOfIteration4);

  using ShrinkFactorsType = itk::FixedArray<unsigned int, 3>;
  ShrinkFactorsType theMovingImageShrinkFactors;
  ShrinkFactorsType theFixedImageShrinkFactors;
  for( int i = 0; i < 3; i++ )
    {
    theMovingImageShrinkFactors[i] = FixedPyramid;
    theFixedImageShrinkFactors[i] = FixedPyramid;
    }

  m_DemonsRegistrator->SetFixedImage( m_DemonsPreprocessor->GetOutputFixedImage() );
  m_DemonsRegistrator->SetMovingImage(
    m_DemonsPreprocessor->GetOutputMovingImage() );
  m_DemonsRegistrator->SetMovingImageShrinkFactors(theMovingImageShrinkFactors);
  m_DemonsRegistrator->SetFixedImageShrinkFactors(theFixedImageShrinkFactors);
  m_DemonsRegistrator->SetNumberOfLevels(NumberOfLevels);
  m_DemonsRegistrator->SetNumberOfIterations(numberOfIterations);
  m_DemonsRegistrator->SetWarpedImageName( this->GetDeformedImageName() );
  // m_DemonsRegistrator->SetDisplacementFieldOutputName(m_DisplacementFieldName);
  m_DemonsRegistrator->SetUseHistogramMatching(true);
}

// EXPORT
void itkSimpleDiffeomorphicRegistration::Update()
{
  InitializePreprocessor();
  try
    {
    m_DemonsPreprocessor->Execute();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
    }
  catch( ... )
    {
    std::cout << "Error occured during preprocessing." << std::endl;
    throw;
    }
  Initialization();
  try
    {
    m_DemonsRegistrator->Execute();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
    }
  catch( ... )
    {
    std::
    cout << "Caught a non-ITK exception " << __FILE__ << " " << __LINE__
         << std::endl;
    }

  m_DisplacementField = m_DemonsRegistrator->GetDisplacementField();
}
