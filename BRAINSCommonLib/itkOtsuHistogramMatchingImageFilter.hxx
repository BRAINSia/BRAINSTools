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
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *    This software is distributed WITHOUT ANY WARRANTY; without even
 *    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *    PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __itkOtsuHistogramMatchingImageFilter_hxx
#define __itkOtsuHistogramMatchingImageFilter_hxx

#include "itkOtsuHistogramMatchingImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include <vector>

#include "BRAINSFitUtils.h"

namespace itk
{
/**
 *
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::OtsuHistogramMatchingImageFilter()
  : m_SourceMinValue(0)
  , m_ReferenceMinValue(0)
  , m_SourceHistogram(HistogramType::New())
  , m_ReferenceHistogram(HistogramType::New())
  , m_OutputHistogram(HistogramType::New())
  , m_SourceMask(nullptr)
  , m_ReferenceMask(nullptr)
{
  this->DynamicMultiThreadingOff(); // NEEDED FOR ITKv5 backwards compatibility
  this->SetNumberOfRequiredInputs(2);

  m_QuantileTable.set_size(3, m_NumberOfMatchPoints + 2);
  m_QuantileTable.fill(0);
  m_Gradients.set_size(m_NumberOfMatchPoints + 1);
  m_Gradients.fill(0);
}

/*
 *
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
void
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::PrintSelf(std::ostream & os,
                                                                                              Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfHistogramLevels: ";
  os << m_NumberOfHistogramLevels << std::endl;
  os << indent << "NumberOfMatchPoints: ";
  os << m_NumberOfMatchPoints << std::endl;

  os << indent << "m_SourceMinValue: ";
  os << m_SourceMinValue << std::endl;
  os << indent << "m_ReferenceMinValue: ";
  os << m_ReferenceMinValue << std::endl;
  os << indent << "m_ReferenceMinValue: ";
  os << m_ReferenceMinValue << std::endl;
  os << indent << "Source histogram: ";
  os << m_SourceHistogram.GetPointer() << std::endl;
  os << indent << "Reference histogram: ";
  os << m_ReferenceHistogram.GetPointer() << std::endl;
  os << indent << "Output histogram: ";
  os << m_OutputHistogram.GetPointer() << std::endl;
  os << indent << "QuantileTable: " << std::endl;
  os << m_QuantileTable << std::endl;
  os << indent << "Gradients: " << std::endl;
  os << m_Gradients << std::endl;
  os << indent << "LowerGradient: ";
  os << m_LowerGradient << std::endl;
  os << indent << "UpperGradient: ";
  os << m_UpperGradient << std::endl;
}

/*
 *
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
void
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::SetReferenceImage(
  const InputImageType * reference)
{
  this->ProcessObject::SetNthInput(1, const_cast<InputImageType *>(reference));
}

/*
 *
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
const typename OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::InputImageType *
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::GetReferenceImage()
{
  if (this->GetNumberOfInputs() < 2)
  {
    return nullptr;
  }

  TInputImage const * const temp = dynamic_cast<TInputImage *>(this->ProcessObject::GetInput(1));
  if (temp == nullptr)
  {
    itkGenericExceptionMacro(<< "Invalid mask converstion attempted.");
  }
  return temp;
}

/*
 * This filter requires all of the input images to be
 * in the buffer.
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
void
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::GenerateInputRequestedRegion()
{
  this->Superclass::GenerateInputRequestedRegion();
  for (unsigned int idx = 0; idx < this->GetNumberOfInputs(); ++idx)
  {
    if (this->GetInput(idx))
    {
      InputImagePointer image = const_cast<InputImageType *>(this->GetInput(idx));
      image->SetRequestedRegionToLargestPossibleRegion();
    }
  }
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
void
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::BeforeThreadedGenerateData()
{
  InputImageConstPointer source = this->GetSourceImage();
  InputImageConstPointer reference = this->GetReferenceImage();

  std::cout << "source stats" << std::endl;

  ComputeRobustMinMaxMean<TInputImage, SpatialObjectType>(
    0.005, source.GetPointer(), m_SourceMask.GetPointer(), m_SourceMinValue, m_SourceMaxValue, m_SourceMeanValue);
  std::cout << "reference stats" << std::endl;
  ComputeRobustMinMaxMean<TInputImage, SpatialObjectType>(0.005,
                                                          reference.GetPointer(),
                                                          m_ReferenceMask.GetPointer(),
                                                          m_ReferenceMinValue,
                                                          m_ReferenceMaxValue,
                                                          m_ReferenceMeanValue);

  this->ConstructHistogram(source, m_SourceMask, m_SourceHistogram, m_SourceMinValue, m_SourceMaxValue);
  this->ConstructHistogram(reference, m_ReferenceMask, m_ReferenceHistogram, m_ReferenceMinValue, m_ReferenceMaxValue);

  // Fill in the quantile table.
  m_QuantileTable.set_size(3, m_NumberOfMatchPoints + 2);
  m_QuantileTable[0][0] = m_SourceMinValue;
  m_QuantileTable[1][0] = m_ReferenceMinValue;

  m_QuantileTable[0][m_NumberOfMatchPoints + 1] = m_SourceMaxValue;
  m_QuantileTable[1][m_NumberOfMatchPoints + 1] = m_ReferenceMaxValue;

  const double delta = 1.0 / (double(m_NumberOfMatchPoints) + 1.0);
  for (unsigned int j = 1; j < m_NumberOfMatchPoints + 1; ++j)
  {
    m_QuantileTable[0][j] = m_SourceHistogram->Quantile(0, double(j) * delta);
    m_QuantileTable[1][j] = m_ReferenceHistogram->Quantile(0, double(j) * delta);
  }

  // Fill in the gradient array.
  m_Gradients.set_size(m_NumberOfMatchPoints + 1);
  for (unsigned int j = 0; j < m_NumberOfMatchPoints + 1; ++j)
  {
    const double denominator = m_QuantileTable[0][j + 1] - m_QuantileTable[0][j];
    if (denominator != 0)
    {
      m_Gradients[j] = m_QuantileTable[1][j + 1] - m_QuantileTable[1][j];
      m_Gradients[j] /= denominator;
    }
    else
    {
      m_Gradients[j] = 0.0;
    }
  }

  {
    const double denominator = m_QuantileTable[0][0] - m_SourceMinValue;
    if (denominator != 0)
    {
      m_LowerGradient = m_QuantileTable[1][0] - m_ReferenceMinValue;
      m_LowerGradient /= denominator;
    }
    else
    {
      m_LowerGradient = 0.0;
    }
  }

  {
    const double denominator = m_QuantileTable[0][m_NumberOfMatchPoints + 1] - m_SourceMaxValue;
    if (denominator != 0)
    {
      m_UpperGradient = m_QuantileTable[1][m_NumberOfMatchPoints + 1] - m_ReferenceMaxValue;
      m_UpperGradient /= denominator;
    }
    else
    {
      m_UpperGradient = 0.0;
    }
  }
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
void
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::AfterThreadedGenerateData()
{
  OutputImagePointer output = this->GetOutput();

  std::cout << "output stats" << std::endl;

  ComputeRobustMinMaxMean<TOutputImage, SpatialObjectType>(
    0.005, output.GetPointer(), this->m_SourceMask.GetPointer(), m_OutputMinValue, m_OutputMaxValue, m_OutputMeanValue);

  this->m_OutputIntensityThreshold = static_cast<OutputPixelType>(m_OutputMinValue);

  this->ConstructHistogram(output, m_SourceMask, m_OutputHistogram, m_OutputIntensityThreshold, m_OutputMaxValue);

  // Fill in the quantile table.
  m_QuantileTable[2][0] = m_OutputIntensityThreshold;

  m_QuantileTable[2][m_NumberOfMatchPoints + 1] = m_OutputMaxValue;

  const double delta = 1.0 / (double(m_NumberOfMatchPoints) + 1.0);
  for (unsigned int j = 1; j < m_NumberOfMatchPoints + 1; ++j)
  {
    m_QuantileTable[2][j] = m_OutputHistogram->Quantile(0, double(j) * delta);
  }
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
void
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  // Get the input and output pointers;
  InputImageConstPointer input = this->GetInput();
  OutputImagePointer     output = this->GetOutput();

  // Transform the source image and write to output.
  using InputConstIterator = ImageRegionConstIterator<InputImageType>;
  using OutputIterator = ImageRegionIterator<OutputImageType>;

  InputConstIterator inIter(input, outputRegionForThread);

  OutputIterator outIter(output, outputRegionForThread);

  // support progress methods/callbacks
  unsigned long updateVisits = 0;
  unsigned long totalPixels = 0;
  if (threadId == 0)
  {
    totalPixels = outputRegionForThread.GetNumberOfPixels();
    updateVisits = totalPixels / 10;
    if (updateVisits < 1)
    {
      updateVisits = 1;
    }
  }
  for (int i = 0; !outIter.IsAtEnd(); ++inIter, ++outIter, ++i)
  {
    if (threadId == 0 && !(i % updateVisits))
    {
      this->UpdateProgress((float)i / (float)totalPixels);
    }

    const double srcValue = static_cast<double>(inIter.Get());

    {
      unsigned int j = 0;
      for (; j < m_NumberOfMatchPoints + 2; ++j)
      {
        if (srcValue < m_QuantileTable[0][j])
        {
          break;
        }
      }

      double mappedValue;
      if (j == 0)
      {
        // Linear interpolate from min to point[0]
        mappedValue = m_ReferenceMinValue + (srcValue - m_SourceMinValue) * m_LowerGradient;
      }
      else if (j == m_NumberOfMatchPoints + 2)
      {
        // Linear interpolate from point[m_NumberOfMatchPoints+1] to max
        mappedValue = m_ReferenceMaxValue + (srcValue - m_SourceMaxValue) * m_UpperGradient;
      }
      else
      {
        // Linear interpolate from point[j] and point[j+1].
        mappedValue = m_QuantileTable[1][j - 1] + (srcValue - m_QuantileTable[0][j - 1]) * m_Gradients[j - 1];
      }

      // Clamp values to the min/max of the source histogram range (the values
      // are supposed to be aligned anyway)
      //      mappedValue=std::max(mappedValue,(double)m_SourceMinValue);
      //      mappedValue=std::min(mappedValue,(double)m_SourceMaxValue);
      outIter.Set(static_cast<OutputPixelType>(mappedValue));
    }
  }
}

/**
 * Construct a histogram from an image.
 */
template <typename TInputImage, typename TOutputImage, typename THistogramMeasurement>
void
OtsuHistogramMatchingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>::ConstructHistogram(
  const InputImageType *           image,
  const SpatialObjectType::Pointer mask,
  HistogramType *                  histogram,
  const THistogramMeasurement      minValue,
  const THistogramMeasurement      maxValue)
{
  {
    // allocate memory for the histogram
    typename HistogramType::SizeType              size;
    typename HistogramType::MeasurementVectorType lowerBound;
    typename HistogramType::MeasurementVectorType upperBound;

    size.SetSize(1);
    lowerBound.SetSize(1);
    upperBound.SetSize(1);
    histogram->SetMeasurementVectorSize(1);
    size[0] = m_NumberOfHistogramLevels;
    lowerBound.Fill(minValue);
    upperBound.Fill(maxValue);

    // Initialize with equally spaced bins.
    histogram->Initialize(size, lowerBound, upperBound);
    histogram->SetToZero();
  }

  typename HistogramType::MeasurementVectorType measurement;
  measurement.SetSize(1);

  using MeasurementType = typename HistogramType::MeasurementType;
  measurement[0] = NumericTraits<MeasurementType>::ZeroValue();

  {
    // put each image pixel into the histogram
    using ConstIterator = ImageRegionConstIteratorWithIndex<InputImageType>;
    ConstIterator iter(image, image->GetBufferedRegion());

    iter.GoToBegin();
    while (!iter.IsAtEnd())
    {
      const InputPixelType value = iter.Get();
      bool                 inMeasurementRegion = false;
      // First check if value is in histogram range.
      if (static_cast<double>(value) >= minValue && static_cast<double>(value) <= maxValue)
      {
        if (mask.IsNull()) // Assume entire area is valid
        {
          inMeasurementRegion = true;
        }
        else
        {
          typename TInputImage::PointType physicalPoint;
          image->TransformIndexToPhysicalPoint(iter.GetIndex(), physicalPoint);
          if (mask.IsNotNull() && mask->IsInsideInWorldSpace(physicalPoint))
          {
            inMeasurementRegion = true;
          }
        }
      }
      if (inMeasurementRegion == true)
      {
        // add sample to histogram
        measurement[0] = value;
        histogram->IncreaseFrequencyOfMeasurement(measurement, 1.0F);
      }
      ++iter;
    }
  }
}
} // end namespace itk

#endif
