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
 *  Module:    $RCSfile: itkScalarImagePortionToHistogramGenerator.h,v $
 *  Language:  C++
 *  Date:      $Date: 2009-08-08 14:18:12 $
 *  Version:   $Revision: 1.2 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __itkScalarImagePortionToHistogramGenerator_h
#define __itkScalarImagePortionToHistogramGenerator_h

#include "itkImageToListSampleFilter.h"
#include "itkSampleToHistogramFilter.h"

#include "itkHistogram.h"
#include "itkObject.h"

namespace itk
{
namespace Statistics
{
/** \class ScalarImagePortionToHistogramGenerator
  *
  * \brief TODO
  */
template <typename TImageType, typename TMaskType>
class ScalarImagePortionToHistogramGenerator : public Object
{
public:
  /** Standard type alias */
  using Self = ScalarImagePortionToHistogramGenerator;
  using Superclass = Object;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScalarImagePortionToHistogramGenerator, Object);

  /** standard New() method support */
  itkNewMacro(Self);

  using ImageType = TImageType;
  using PixelType = typename TImageType::PixelType;
  using RealPixelType = typename NumericTraits<PixelType>::RealType;

  using HistogramType = itk::Statistics::Histogram<double>;
  using HistogramPointer = typename HistogramType::Pointer;
  using HistogramConstPointer = typename HistogramType::ConstPointer;

  using ListSampleGeneratorType = itk::Statistics::ImageToListSampleFilter<ImageType, TMaskType>;
  using ListSampleGeneratorPointer = typename ListSampleGeneratorType::Pointer;
  using ListSampleType = typename ListSampleGeneratorType::ListSampleType;

  using GeneratorType = itk::Statistics::SampleToHistogramFilter<ListSampleType, HistogramType>;
  using GeneratorPointer = typename GeneratorType::Pointer;
public:

  /** Triggers the Computation of the histogram */
  void Compute(void);

  /** Connects the input image for which the histogram is going to be computed
    */
  void SetInput(const TImageType *);

  /** Connects the input image for which the histogram is going to be computed
    */
  void SetBinaryPortionImage(const TMaskType *);

  /** Return the histogram.
    * \warning This output is only valid after the Compute() method has been
    *    invoked
    * \sa Compute */
  const HistogramType * GetOutput() const;

  /** Set number of histogram bins */
  void SetNumberOfBins(unsigned int numberOfBins);

  /** Set marginal scale value to be passed to the histogram generator */
  void SetMarginalScale(double marginalScale);

  /** Set the minimum value from which the bins will be computed */
  void SetHistogramMin(RealPixelType minimumValue);

  /** Set the maximum value from which the bins will be computed */
  void SetHistogramMax(RealPixelType maximumValue);

protected:
  ScalarImagePortionToHistogramGenerator();
  ~ScalarImagePortionToHistogramGenerator() override
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const override;

private:

  ListSampleGeneratorPointer m_ImageToListSampleGenerator;

  HistogramPointer m_Histogram;
  GeneratorPointer m_HistogramGenerator;

  ScalarImagePortionToHistogramGenerator(const Self &); // purposely not
                                                        // implemented
  void operator=(const Self &);                         // purposely not

  // implemented
};
}   // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScalarImagePortionToHistogramGenerator.hxx"
#endif

#endif
