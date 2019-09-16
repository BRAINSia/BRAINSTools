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
/*
 *  itkMultiResolutionICCDeformableRegistration.txx
 *  iccdefRegistrationNew
 *
 *  Created by Yongqiang Zhao on 5/8/09.
 *  Copyright 2009 __UI__. All rights reserved.
 *
 */

#ifndef __itkMultiResolutionICCDeformableRegistration_hxx
#define __itkMultiResolutionICCDeformableRegistration_hxx
#include "itkMultiResolutionICCDeformableRegistration.h"

#include "itkRecursiveGaussianImageFilter.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkMath.h"

#include "itkWarpImageFilter.h"

namespace itk
{
/**
 * Default constructor
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::
  MultiResolutionICCDeformableRegistration()
{
  this->SetNumberOfRequiredInputs(2);
  // Primary input is optional in this filter
  this->RemoveRequiredInputName("Primary");

  typename RegistrationType::Pointer registrator = RegistrationType::New();
  registrator->InPlaceOn();
  m_RegistrationFilter = dynamic_cast<RegistrationType *>(registrator.GetPointer());

  m_MovingImagePyramid = MovingImagePyramidType::New();
  m_FixedImagePyramid = FixedImagePyramidType::New();
  m_FieldExpander12 = FieldExpanderType::New();
  m_FieldExpander21 = FieldExpanderType::New();
  m_InitialDisplacementField = nullptr;

  this->SetNumberOfRequiredOutputs(2);
  this->SetNthOutput(0, this->MakeOutput(0));
  this->SetNthOutput(1, this->MakeOutput(1));

  m_NumberOfLevels = 3;
  m_NumberOfIterations.SetSize(m_NumberOfLevels);
  m_FixedImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
  m_MovingImagePyramid->SetNumberOfLevels(m_NumberOfLevels);

  unsigned int ilevel;
  for (ilevel = 0; ilevel < m_NumberOfLevels; ilevel++)
  {
    m_NumberOfIterations[ilevel] = 10;
  }
  m_CurrentLevel = 0;

  m_StopRegistrationFlag = false;
  m_DisplacementFieldOutputNamePrefix = "none";
}

/*
 * Set the moving image image.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::SetMovingImage(
  const MovingImageType * ptr)
{
  this->ProcessObject::SetNthInput(2, const_cast<MovingImageType *>(ptr));
}

/*
 * Get the moving image image.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
const typename MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::
  MovingImageType *
  MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::GetMovingImage(
    void) const
{
  return dynamic_cast<const MovingImageType *>(this->ProcessObject::GetInput(2));
}

/*
 * Set the fixed image.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::SetFixedImage(
  const FixedImageType * ptr)
{
  this->ProcessObject::SetNthInput(1, const_cast<FixedImageType *>(ptr));
}

/*
 * Get the fixed image.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
const typename MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::
  FixedImageType *
  MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::GetFixedImage(
    void) const
{
  return dynamic_cast<const FixedImageType *>(this->ProcessObject::GetInput(1));
}

/*
 *
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
std::vector<SmartPointer<DataObject>>::size_type
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::
  GetNumberOfValidRequiredInputs() const
{
  typename std::vector<SmartPointer<DataObject>>::size_type num = 0;

  if (this->GetFixedImage())
  {
    num++;
  }

  if (this->GetMovingImage())
  {
    num++;
  }

  return num;
}

/**
 * Set the number of multi-resolution levels
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::SetNumberOfLevels(
  unsigned int num)
{
  if (m_NumberOfLevels != num)
  {
    this->Modified();
    m_NumberOfLevels = num;
    m_NumberOfIterations.SetSize(m_NumberOfLevels);
  }

  if (m_MovingImagePyramid && m_MovingImagePyramid->GetNumberOfLevels() != num)
  {
    m_MovingImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
  }
  if (m_FixedImagePyramid && m_FixedImagePyramid->GetNumberOfLevels() != num)
  {
    m_FixedImagePyramid->SetNumberOfLevels(m_NumberOfLevels);
  }
}

/**
 * Standard PrintSelf method.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfLevels: " << m_NumberOfLevels << std::endl;
  os << indent << "CurrentLevel: " << m_CurrentLevel << std::endl;

  os << indent << "NumberOfIterations: [";
  unsigned int ilevel;
  for (ilevel = 0; ilevel < m_NumberOfLevels - 1; ilevel++)
  {
    os << m_NumberOfIterations[ilevel] << ", ";
  }
  os << m_NumberOfIterations[ilevel] << "]" << std::endl;

  os << indent << "RegistrationFilter: ";
  os << m_RegistrationFilter.GetPointer() << std::endl;
  os << indent << "MovingImagePyramid: ";
  os << m_MovingImagePyramid.GetPointer() << std::endl;
  os << indent << "FixedImagePyramid: ";
  os << m_FixedImagePyramid.GetPointer() << std::endl;

  os << indent << "StopRegistrationFlag: ";
  os << m_StopRegistrationFlag << std::endl;
}

/*
 * Perform a the deformable registration using a multiresolution scheme
 * using an internal mini-pipeline
 *
 *  ref_pyramid ->  registrator  ->  field_expander --|| tempField
 * test_pyramid ->           |                              |
 *                           |                              |
 *                           --------------------------------
 *
 * A tempField image is used to break the cycle between the
 * registrator and field_expander.
 *
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::GenerateData()
{
  // Check for NULL images and pointers
  MovingImageConstPointer movingImage = this->GetMovingImage();
  FixedImageConstPointer  fixedImage = this->GetFixedImage();

  if (!movingImage || !fixedImage)
  {
    itkExceptionMacro(<< "Fixed and/or moving image not set");
  }

  if (!m_MovingImagePyramid || !m_FixedImagePyramid)
  {
    itkExceptionMacro(<< "Fixed and/or moving pyramid not set");
  }

  if (!m_RegistrationFilter)
  {
    itkExceptionMacro(<< "Registration filter not set");
  }

  if (this->m_InitialDisplacementField && this->GetInput(0))
  {
    itkExceptionMacro(<< "Only one initial deformation can be given. "
                      << "SetInitialDisplacementField should not be used in "
                      << "cunjunction with SetArbitraryInitialDisplacementField "
                      << "or SetInput.");
  }

  // as per suggestion in this bug report:
  // http://public.kitware.com/Bug/view.php?id=3590
  // this should allow input images to be released, since
  // they are no longer needed after generating the image pyramid
  this->RestoreInputReleaseDataFlags();

  // Create the image pyramids.
  m_MovingImagePyramid->SetInput(movingImage);
  m_MovingImagePyramid->UpdateLargestPossibleRegion();

  m_FixedImagePyramid->SetInput(fixedImage);
  m_FixedImagePyramid->UpdateLargestPossibleRegion();

  // Initializations
  m_CurrentLevel = 0;
  m_StopRegistrationFlag = false;

  unsigned int movingLevel = std::min((int)m_CurrentLevel, (int)m_MovingImagePyramid->GetNumberOfLevels());

  unsigned int fixedLevel = std::min((int)m_CurrentLevel, (int)m_FixedImagePyramid->GetNumberOfLevels());

  DisplacementFieldPointer tempField12 = this->m_InitialMovingDisplacementField;
  DisplacementFieldPointer tempField21 = this->m_InitialFixedDisplacementField;
  bool                     lastShrinkFactorsAllOnes = false;

  int iteration = 0;

  while (!this->Halt())
  {
    if (tempField12.IsNull())
    {
      m_RegistrationFilter->SetInitialDisplacementField(nullptr);
    }
    else
    {
      // Resample the field to be the same size as the fixed image
      // at the current level
      m_FieldExpander12->SetInput(tempField12);

      typename FloatImageType::Pointer fi = m_FixedImagePyramid->GetOutput(fixedLevel);
      m_FieldExpander12->SetSize(fi->GetLargestPossibleRegion().GetSize());
      m_FieldExpander12->SetOutputStartIndex(fi->GetLargestPossibleRegion().GetIndex());
      m_FieldExpander12->SetOutputOrigin(fi->GetOrigin());
      m_FieldExpander12->SetOutputSpacing(fi->GetSpacing());

      m_FieldExpander12->UpdateLargestPossibleRegion();
      // m_FieldExpander12->SetInput( NULL );
      tempField12 = m_FieldExpander12->GetOutput();
      tempField12->DisconnectPipeline();

      m_RegistrationFilter->SetInitialForwardDisplacementField(tempField12);

      m_FieldExpander21->SetInput(tempField21);
      typename FloatImageType::Pointer mi = m_MovingImagePyramid->GetOutput(movingLevel);
      m_FieldExpander21->SetSize(mi->GetLargestPossibleRegion().GetSize());
      m_FieldExpander21->SetOutputStartIndex(mi->GetLargestPossibleRegion().GetIndex());
      m_FieldExpander21->SetOutputOrigin(mi->GetOrigin());
      m_FieldExpander21->SetOutputSpacing(mi->GetSpacing());

      m_FieldExpander21->UpdateLargestPossibleRegion();
      // m_FieldExpander21->SetInput( NULL );
      tempField21 = m_FieldExpander21->GetOutput();
      tempField21->DisconnectPipeline();

      m_RegistrationFilter->SetInitialBackwardDisplacementField(tempField21);
    }

    // setup registration filter and pyramids
    m_RegistrationFilter->SetMovingImage(m_MovingImagePyramid->GetOutput(movingLevel));
    m_RegistrationFilter->SetFixedImage(m_FixedImagePyramid->GetOutput(fixedLevel));

    m_RegistrationFilter->SetNumberOfIterations(m_NumberOfIterations[m_CurrentLevel]);

    // cache shrink factors for computing the next expand factors.
    lastShrinkFactorsAllOnes = true;
    for (unsigned int idim = 0; idim < ImageDimension; idim++)
    {
      if (m_FixedImagePyramid->GetSchedule()[fixedLevel][idim] > 1)
      {
        lastShrinkFactorsAllOnes = false;
        break;
      }
    }

    // compute new deformation field
    m_RegistrationFilter->UpdateLargestPossibleRegion();
    tempField12 = m_RegistrationFilter->GetOutput(0);
    tempField12->DisconnectPipeline();
    tempField21 = m_RegistrationFilter->GetOutput(1);
    tempField21->DisconnectPipeline();
    //  tempField12->Print(std::cout,6);

    if (m_DisplacementFieldOutputNamePrefix != "" && m_DisplacementFieldOutputNamePrefix != "none")
    {
      std::string       name = m_DisplacementFieldOutputNamePrefix + "_iteration"; // resolution";
      std::stringstream ss1;
      std::string       str1;
      iteration += m_NumberOfIterations[m_CurrentLevel];
      ss1 << iteration;
      ss1 >> str1;
      name = name + str1;

      std::stringstream ss;
      std::string       str;
      int               level = m_NumberOfLevels - m_CurrentLevel - 1;
      //  int resolution = 10000*(std::pow(0.5, level ));
      //  std::cout<<level<<std::endl;
      //  std::cout<<"r:"<<resolution<<std::endl;
      ss << 10000 * (std::pow(0.5, level));
      ss >> str;
      name = name + "resolution" + str;

      using WriteImageType = ImageFileWriter<TDisplacementField>;
      typename WriteImageType::Pointer writer12 = WriteImageType::New();
      writer12->SetInput(tempField12);
      writer12->SetFileName("forward/" + name + "_forward.nii.gz");
      writer12->Update();

      typename WriteImageType::Pointer writer21 = WriteImageType::New();
      writer21->SetInput(tempField21);
      writer21->SetFileName("backward/" + name + "_backward.nii.gz");
      writer21->Update();
    }

    // Increment level counter.
    std::cout << "Level: " << m_CurrentLevel << std::endl;
    m_CurrentLevel++;
    movingLevel = std::min((int)m_CurrentLevel, (int)m_MovingImagePyramid->GetNumberOfLevels());
    fixedLevel = std::min((int)m_CurrentLevel, (int)m_FixedImagePyramid->GetNumberOfLevels());

    // Invoke an iteration event.
    this->InvokeEvent(IterationEvent());

    // We can release data from pyramid which are no longer required.
    if (movingLevel > 0)
    {
      m_MovingImagePyramid->GetOutput(movingLevel - 1)->ReleaseData();
    }
    if (fixedLevel > 0)
    {
      m_FixedImagePyramid->GetOutput(fixedLevel - 1)->ReleaseData();
    }
  } // while not Halt()

  if (!lastShrinkFactorsAllOnes)
  {
    // Some of the last shrink factors are not one
    // graft the output of the expander filter to
    // to output of this filter

    // resample the field to the same size as the fixed image
    m_FieldExpander12->SetInput(tempField12);
    m_FieldExpander12->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
    m_FieldExpander12->SetOutputStartIndex(fixedImage->GetLargestPossibleRegion().GetIndex());
    m_FieldExpander12->SetOutputOrigin(fixedImage->GetOrigin());
    m_FieldExpander12->SetOutputSpacing(fixedImage->GetSpacing());
    m_FieldExpander12->SetOutputDirection(fixedImage->GetDirection());
    m_FieldExpander12->UpdateLargestPossibleRegion();
    this->GraftNthOutput(0, m_FieldExpander12->GetOutput());
#if 1
    m_FieldExpander21->SetInput(tempField21);
    m_FieldExpander21->SetSize(movingImage->GetLargestPossibleRegion().GetSize());
    m_FieldExpander21->SetOutputStartIndex(movingImage->GetLargestPossibleRegion().GetIndex());
    m_FieldExpander21->SetOutputOrigin(movingImage->GetOrigin());
    m_FieldExpander21->SetOutputSpacing(movingImage->GetSpacing());
    m_FieldExpander21->SetOutputDirection(movingImage->GetDirection());
    m_FieldExpander21->UpdateLargestPossibleRegion();
    //    tempField21 = m_FieldExpander21->GetOutput();
    //    tempField21->DisconnectPipeline();
    //    tempField21->Print(std::cout,6);
    this->GraftNthOutput(1, m_FieldExpander21->GetOutput());
#endif
  }
  else
  {
    // all the last shrink factors are all ones
    // graft the output of registration filter to
    // to output of this filter
    this->GraftNthOutput(0, tempField12);
    this->GraftNthOutput(1, tempField21);
  }

  // Release memory
  m_FieldExpander12->SetInput(nullptr);
  m_FieldExpander12->GetOutput()->ReleaseData();
  m_FieldExpander21->SetInput(nullptr);
  m_FieldExpander21->GetOutput()->ReleaseData();
  m_RegistrationFilter->SetInput(nullptr);
  m_RegistrationFilter->GetOutput(0)->ReleaseData();
  m_RegistrationFilter->GetOutput(1)->ReleaseData();
}

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::StopRegistration()
{
  m_RegistrationFilter->StopRegistration();
  m_StopRegistrationFlag = true;
}

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
bool
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::Halt()
{
  // Halt the registration after the user-specified number of levels
  if (m_NumberOfLevels != 0)
  {
    this->UpdateProgress(static_cast<float>(m_CurrentLevel) / static_cast<float>(m_NumberOfLevels));
  }

  if (m_CurrentLevel >= m_NumberOfLevels)
  {
    return true;
  }
  if (m_StopRegistrationFlag)
  {
    return true;
  }
  else if (m_NumberOfIterations[m_CurrentLevel] == 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::
  GenerateOutputInformation()
{
  typename DataObject::Pointer output;

  if (this->GetInput(0))
  {
    // Initial deformation field is set.
    // Copy information from initial field.
    this->Superclass::GenerateOutputInformation();
  }
  else if (this->GetFixedImage())
  {
    // Initial deforamtion field is not set.
    // Copy information from the fixed image.
    for (unsigned int idx = 0; idx < this->GetNumberOfIndexedOutputs(); ++idx)
    {
      output = this->GetOutput(idx);
      if (output)
      {
        output->CopyInformation(this->GetFixedImage());
      }
    }
  }
}

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::
  GenerateInputRequestedRegion()
{
  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the moving image
  MovingImagePointer movingPtr = const_cast<MovingImageType *>(this->GetMovingImage());

  if (movingPtr)
  {
    movingPtr->SetRequestedRegionToLargestPossibleRegion();
  }

  // just propagate up the output requested region for
  // the fixed image and initial deformation field.
  DisplacementFieldPointer inputPtr = const_cast<DisplacementFieldType *>(this->GetInput());
  DisplacementFieldPointer outputPtr = this->GetOutput();
  FixedImagePointer        fixedPtr = const_cast<FixedImageType *>(this->GetFixedImage());

  if (inputPtr)
  {
    inputPtr->SetRequestedRegion(outputPtr->GetRequestedRegion());
  }

  if (fixedPtr)
  {
    fixedPtr->SetRequestedRegion(outputPtr->GetRequestedRegion());
  }
}

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
void
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::
  EnlargeOutputRequestedRegion(DataObject * ptr)
{
  // call the superclass's implementation
  Superclass::EnlargeOutputRequestedRegion(ptr);

  // set the output requested region to largest possible.
  DisplacementFieldType * outputPtr;

  outputPtr = dynamic_cast<DisplacementFieldType *>(ptr);

  if (outputPtr)
  {
    outputPtr->SetRequestedRegionToLargestPossibleRegion();
  }
}

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TRealType>
ProcessObject::DataObjectPointer
MultiResolutionICCDeformableRegistration<TFixedImage, TMovingImage, TDisplacementField, TRealType>::MakeOutput(
  ProcessObject::DataObjectPointerArraySizeType idx)
{
  switch (idx)
  {
    case 0:
    {
      return static_cast<DataObject *>(TDisplacementField::New().GetPointer());
    }
    break;
    case 1:
    {
      return static_cast<DataObject *>(TDisplacementField::New().GetPointer());
    }
    break;
    default:
      return nullptr;
  }
}
} // end namespace itk

#endif
