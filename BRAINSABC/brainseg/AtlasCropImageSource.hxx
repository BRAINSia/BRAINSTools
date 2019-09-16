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
#ifndef __AtlasCropImageSource_hxx
#define __AtlasCropImageSource_hxx

#include "itkImageRegionIteratorWithIndex.h"

#include "AtlasCropImageSource.h"

template <typename TInputImage, typename TProbabilityImage>
AtlasCropImageSource<TInputImage, TProbabilityImage>::AtlasCropImageSource()
{
  m_Padding = 8.0;

  m_LowerBound.Fill(0);
  m_UpperBound.Fill(0);

  m_OriginalSize.Fill(0);

  m_SlabMode = false;
}

template <typename TInputImage, typename TProbabilityImage>
bool
AtlasCropImageSource<TInputImage, TProbabilityImage>::CheckBounds()
{
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (m_LowerBound[i] > m_UpperBound[i])
    {
      return false;
    }
  }

  InputImageSizeType croppedSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    croppedSize[i] = m_UpperBound[i] - m_LowerBound[i] + 1;
  }
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (croppedSize[i] > m_OriginalSize[i])
    {
      return false;
    }
  }

  return true;
}

template <typename TInputImage, typename TProbabilityImage>
void
AtlasCropImageSource<TInputImage, TProbabilityImage>::UseProbabilities(ProbabilityImageList probs)
{
  if (probs.size() == 0)
  {
    itkExceptionMacro(<< "Need at least one class probability image");
  }

  ProbabilityImageSizeType size = probs[0]->GetLargestPossibleRegion().size();

  ProbabilityImageSpacingType spacing = probs[0]->GetSpacing();
  // Make sure all class probabilities have the same space
  for (unsigned int i = 1; i < probs.size(); i++)
  {
    ProbabilityImageSizeType othersize = probs[i]->GetLargestPossibleRegion().size();
    if (size != othersize)
    {
      itkExceptionMacro(<< "Probability list size mismatch");
    }
  }

  // Transform padding to voxel counts
  InputImageOffsetType padding;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    padding[i] = (unsigned int)std::floor(m_Padding / spacing[i] + 0.5);
  }
  // Make sure padding is sensible
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if ((long)size[i] <= padding[i])
    {
      itkExceptionMacro(<< "Bounding box padding larger than or equal to image size");
    }
  }

  // Save size info
  m_OriginalSize = size;
  // Initial bounds: whole image
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    m_LowerBound[i] = size[i] - 1;
    m_UpperBound[i] = 0;
  }

  // Go through image and update bounds
  using IteratorType = itk::ImageRegionIteratorWithIndex<ProbabilityImageType>;

  IteratorType it(probs[0], probs[0]->GetLargestPossibleRegion());

  it.GoToBegin();
  while (!it.IsAtEnd())
  {
    ProbabilityImageIndexType ind = it.GetIndex();

    double sumProb = 0;
    for (unsigned int i = 0; i < probs.size(); i++)
    {
      sumProb += probs[i]->GetPixel(ind);
    }

    if (sumProb > 0)
    {
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        if (ind[i] < m_LowerBound[i])
        {
          m_LowerBound[i] = ind[i];
        }
        if (ind[i] > m_UpperBound[i])
        {
          m_UpperBound[i] = ind[i];
        }
      }
    }

    ++it;
  }

  // Only generate a slab?
  if (m_SlabMode)
  {
    long len = m_UpperBound[ImageDimension - 1] - m_LowerBound[ImageDimension - 1] + 1;
    long offt = (long)(0.2 * len);
    m_LowerBound[ImageDimension - 1] += offt;
    m_UpperBound[ImageDimension - 1] -= offt;
  }
  // Enlarge/pad bounding box
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (m_LowerBound[i] > padding[i])
    {
      m_LowerBound[i] -= padding[i];
    }
    else
    {
      m_LowerBound[i] = 0;
    }

    if (m_UpperBound[i] < ((long)size[i] - 1 - padding[i]))
    {
      m_UpperBound[i] += padding[i];
    }
    else
    {
      m_UpperBound[i] = size[i] - 1;
    }
  }

  // Store origin information
  m_InputOrigin = probs[0]->GetOrigin();
  probs[0]->TransformIndexToPhysicalPoint(m_LowerBound, m_CropOrigin);

  // m_CropInfo.offset =
  // m_CropInfo.size =
}

template <typename TInputImage, typename TProbabilityImage>
typename AtlasCropImageSource<TInputImage, TProbabilityImage>::InputImagePointer
AtlasCropImageSource<TInputImage, TProbabilityImage>::Restore(InputImagePointer img)
{
  if (!this->CheckBounds())
  {
    itkExceptionMacro(<< "Invalid bounds");
  }

  InputImageSizeType size = img->GetLargestPossibleRegion().size();

  InputImageSizeType croppedSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    croppedSize[i] = m_UpperBound[i] - m_LowerBound[i] + 1;
  }

  if (size != croppedSize)
  {
    itkExceptionMacro(<< "Input size does not match size of cropped image");
  }

  InputImageRegionType region;
  region.SetSize(m_OriginalSize);

  InputImagePointer output = InputImageType::New();
  output->CopyInformation(img);
  // Adjust origin
  output->SetOrigin(m_InputOrigin);
  output->SetRegions(region);
  output->Allocate();
  output->FillBuffer(0);

  InputImageOffsetType offt;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    offt[i] = m_LowerBound[i];
  }

  using IteratorType = itk::ImageRegionIteratorWithIndex<InputImageType>;

  IteratorType inIt(img, img->GetLargestPossibleRegion());

  inIt.GoToBegin();
  while (!inIt.IsAtEnd())
  {
    InputImageIndexType ind = inIt.GetIndex();

    output->SetPixel(ind + offt, inIt.Get());

    ++inIt;
  }

  return output;
}

template <typename TInputImage, typename TProbabilityImage>
typename AtlasCropImageSource<TInputImage, TProbabilityImage>::InputImagePointer
AtlasCropImageSource<TInputImage, TProbabilityImage>::Crop(InputImagePointer img)
{
  if (!this->CheckBounds())
  {
    itkExceptionMacro(<< "Invalid bounds");
  }

  InputImageSizeType size = img->GetLargestPossibleRegion().size();

  if (size != m_OriginalSize)
  {
    itkExceptionMacro(<< "Input size does not match size of probability images");
  }

  InputImageSizeType croppedSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    croppedSize[i] = m_UpperBound[i] - m_LowerBound[i] + 1;
  }

  InputImageRegionType cropRegion;
  cropRegion.SetSize(croppedSize);

  InputImagePointer output = InputImageType::New();
  output->CopyInformation(img);
  // Adjust origin
  output->SetOrigin(m_CropOrigin);
  output->SetRegions(cropRegion);
  output->Allocate();

  InputImageOffsetType offt;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    offt[i] = m_LowerBound[i];
  }

  using IteratorType = itk::ImageRegionIteratorWithIndex<InputImageType>;

  IteratorType outIt(output, cropRegion);

  outIt.GoToBegin();
  while (!outIt.IsAtEnd())
  {
    InputImageIndexType ind = outIt.GetIndex();

    outIt.Set(img->GetPixel(ind + offt));

    ++outIt;
  }

  return output;
}

#endif
