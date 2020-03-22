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
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *    This software is distributed WITHOUT ANY WARRANTY; without even
 *    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *    PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __itkBRAINSROIAutoImageFilter_hxx
#define __itkBRAINSROIAutoImageFilter_hxx
#include "itkBRAINSROIAutoImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkProgressAccumulator.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
BRAINSROIAutoImageFilter<TInputImage, TOutputImage>::BRAINSROIAutoImageFilter()
  : m_ResultMaskPointer(nullptr)
{
  // this filter requires two input images
  this->SetNumberOfRequiredInputs(1);
}

template <typename TInputImage, typename TOutputImage>
void
BRAINSROIAutoImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  m_ResultMaskPointer = nullptr; // Need to make this null during every re-run of
                                 // the data.
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  using LFFMaskFilterType = itk::LargestForegroundFilledMaskImageFilter<TInputImage, TOutputImage>;
  typename LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
  // Register the filter with the with progress accumulator using
  // equal weight proportion
  progress->RegisterInternalFilter(LFF, 1.0f);
  LFF->SetInput(this->GetInput());
  LFF->SetOtsuPercentileThreshold(m_OtsuPercentileThreshold);
  LFF->SetClosingSize(m_ClosingSize);
  LFF->SetDilateSize(m_DilateSize);
  LFF->SetThresholdCorrectionFactor(m_ThresholdCorrectionFactor);
  LFF->Update();
  this->GraftOutput(LFF->GetOutput());
}

template <typename TInputImage, typename TOutputImage>
void
BRAINSROIAutoImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "OtsuPercentileThreshold: " << m_OtsuPercentileThreshold << std::endl;
  os << indent << "ThresholdCorrectionFactor: " << m_ThresholdCorrectionFactor << std::endl;
  os << indent << "ClosingSize: " << m_ClosingSize << std::endl;
  os << indent << "DilateSize: " << m_DilateSize << std::endl;
}
} // end namespace itk
#endif
