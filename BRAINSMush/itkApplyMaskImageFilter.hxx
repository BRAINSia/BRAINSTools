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
#ifndef _itkApplyMaskImageFilter_hxx
#define _itkApplyMaskImageFilter_hxx

#include "itkApplyMaskImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
ApplyMaskImageFilter<TInputImage, TOutputImage>::ApplyMaskImageFilter()
{
  this->SetNumberOfRequiredInputs(2);

  m_InvertMask = false;
}

template <typename TInputImage, typename TOutputImage>
void
ApplyMaskImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "InvertMask: ";
  os << m_InvertMask << std::endl;
}

template <typename TInputImage, typename TOutputImage>
void
ApplyMaskImageFilter<TInputImage, TOutputImage>::SetMaskImage(const InputImageType * reference)
{
  this->ProcessObject::SetNthInput(1, const_cast<InputImageType *>(reference));
}

template <typename TInputImage, typename TOutputImage>
const typename ApplyMaskImageFilter<TInputImage, TOutputImage>::InputImageType *
ApplyMaskImageFilter<TInputImage, TOutputImage>::GetMaskImage()
{
  if (this->GetNumberOfInputs() < 2)
  {
    return NULL;
  }

  return dynamic_cast<TInputImage *>(this->ProcessObject::GetInput(1));
}

template <typename TInputImage, typename TOutputImage>
void
ApplyMaskImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  this->AllocateOutputs();

  // get the input and output pointers
  InputImageConstPointer inputImage = this->GetInput(0);
  InputImageConstPointer maskImage = this->GetInput(1);
  OutputImagePointer     outputVolume = this->GetOutput();

  // mask the image
  using InputConstIterator = ImageRegionConstIterator<InputImageType>;
  using OutputIterator = ImageRegionIterator<OutputImageType>;

  InputConstIterator itImage(inputImage, inputImage->GetRequestedRegion());

  InputConstIterator itMask(maskImage, maskImage->GetRequestedRegion());

  OutputIterator itOut(outputVolume, outputVolume->GetRequestedRegion());

  // support progress methods/callbacks
  unsigned long updateVisits = 0;
  unsigned long totalPixels = 0;

  totalPixels = inputImage->GetRequestedRegion().GetNumberOfPixels();
  updateVisits = totalPixels / 10;
  if (updateVisits < 1)
  {
    updateVisits = 1;
  }

  itImage.GoToBegin();
  itMask.GoToBegin();
  unsigned long i = 0;
  for (itOut.GoToBegin(); !itOut.IsAtEnd(); ++itOut)
  {
    if (!(i % updateVisits))
    {
      this->UpdateProgress(static_cast<float>(i )/ static_cast<float>(totalPixels));
    }

    if (!m_InvertMask)
    {
      if (itMask.Get() != 0)
      {
        itOut.Set(static_cast<OutputPixelType>(itImage.Get()));
      }
      else
      {
        itOut.Set(static_cast<OutputPixelType>(0));
      }
    }
    else
    {
      if (itMask.Get() != 0)
      {
        itOut.Set(static_cast<OutputPixelType>(0));
      }
      else
      {
        itOut.Set(static_cast<OutputPixelType>(itImage.Get()));
      }
    }
    ++i;
    ++itMask;
    ++itImage;
  }
}
} // end namespace itk

#endif
