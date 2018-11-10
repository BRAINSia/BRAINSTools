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
#ifndef CleanBrainLabelMap_h
#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkVotingBinaryHoleFillingImageFilter.h"

template <typename TInputImage, typename TOutputImage>
typename TOutputImage::Pointer
CleanBrainLabelMap(const TInputImage *inputImage)
{
  using BinaryThresholdFilterType = typename
    itk::BinaryThresholdImageFilter<TInputImage, TInputImage>;

  typename BinaryThresholdFilterType::Pointer binaryThresholdFilter =
    BinaryThresholdFilterType::New();
  binaryThresholdFilter->SetLowerThreshold(1);
  binaryThresholdFilter->SetUpperThreshold(255);
  binaryThresholdFilter->SetInput(inputImage);
  binaryThresholdFilter->Update();
  typename TInputImage::Pointer emsBrainMask(binaryThresholdFilter->GetOutput() );

  using KernelType = typename itk::FlatStructuringElement<TInputImage::ImageDimension>;
  typename KernelType::RadiusType erodeRadius = { { 2, 2, 2 } };
  KernelType erodeKernel = KernelType::Ball(erodeRadius);

  using BinaryErodeFilterType = typename itk::BinaryErodeImageFilter<TInputImage, TInputImage, KernelType>;
  typename BinaryErodeFilterType::Pointer erodeFilter = BinaryErodeFilterType::New();

  erodeFilter->SetInput(emsBrainMask);
  erodeFilter->SetKernel(erodeKernel);

  using RelabelComponentFilterType = typename itk::RelabelComponentImageFilter<TInputImage, TInputImage>;
  typename RelabelComponentFilterType::Pointer relabelFilter =
    RelabelComponentFilterType::New();
  relabelFilter->SetInput(erodeFilter->GetOutput() );
  relabelFilter->SetMinimumObjectSize(30000);

  typename KernelType::RadiusType dilateRadius = { { 4, 4, 4 } };
  KernelType dilateKernel = KernelType::Ball(dilateRadius);

  using BinaryDilateFilterType = typename itk::BinaryDilateImageFilter<TInputImage, TInputImage, KernelType>;
  typename BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();

  dilateFilter->SetKernel(dilateKernel);
  dilateFilter->SetInput(relabelFilter->GetOutput() );

  using AndFilterType = typename itk::AndImageFilter<TInputImage, TInputImage, TInputImage>;
  typename AndFilterType::Pointer andFilter = AndFilterType::New();
  andFilter->SetInput1(emsBrainMask);
  andFilter->SetInput2(dilateFilter->GetOutput() );

  typename TInputImage::SizeType holeFillingRadius = { { 3, 3, 3 } };

  using HoleFillingFilterType = typename itk::VotingBinaryHoleFillingImageFilter<TInputImage, TOutputImage>;
  typename HoleFillingFilterType::Pointer holeFillingFilter =
    HoleFillingFilterType::New();
  holeFillingFilter->SetInput(andFilter->GetOutput() );
  holeFillingFilter->SetRadius(holeFillingRadius);
  holeFillingFilter->SetForegroundValue(1);
  holeFillingFilter->SetBackgroundValue(0);
  holeFillingFilter->Update();
  return holeFillingFilter->GetOutput();
}

#endif // CleanBrainLabelMap_h
