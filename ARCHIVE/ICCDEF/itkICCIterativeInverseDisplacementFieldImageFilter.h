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

#ifndef __itkICCIterativeInverseDisplacementFieldImageFilter_h
#define __itkICCIterativeInverseDisplacementFieldImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkWarpVectorImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkTimeProbe.h"

namespace itk
{
/** \class ICCIterativeInverseDisplacementFieldImageFilter
 * \brief Computes the inverse of a deformation field.
 *
 * ICCIterativeInverseDisplacementFieldImageFilter takes a deformation field as input and
 * computes the deformation field that is its inverse. If the input deformation
 * field was mapping coordinates from a space A into a space B, the output of
 * this filter will map coordinates from the space B into the space A.
 *
 * The algorithm implemented in this filter uses an iterative method for
 * progresively refining the values of the inverse field. Starting from the
 * direct field, at every pixel the direct mapping of this point is found, and
 * a the nevative of the current deformation is stored in the inverse field at
 * the nearest pixel. Then, subsequent iterations verify if any of the neigbor pixels
 * provide a better return to the current pixel, in which case its value is taken for
 * updating the vector in the inverse field.
 *
 * This method was discussed in the users-list during February 2004.
 *
 * \author  Corinne Mattmann
 *
 */

template <typename TInputImage, typename TOutputImage>
class ICCIterativeInverseDisplacementFieldImageFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(ICCIterativeInverseDisplacementFieldImageFilter);

  /** Standard class type alias. */
  using Self = ICCIterativeInverseDisplacementFieldImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ICCIterativeInverseDisplacementFieldImageFilter, ImageToImageFilter);

  /** Some type alias. */
  using InputImageType = TInputImage;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImagePointType = typename InputImageType::PointType;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputImageSizeType = typename InputImageType::SizeType;
  using InputImageSpacingType = typename InputImageType::SpacingType;
  using InputImageIndexType = typename InputImageType::IndexType;
  using PixelType = typename InputImageType::PixelType;
  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputImagePixelType = typename OutputImageType::PixelType;
  using OutputImagePointType = typename OutputImageType::PointType;
  using OutputImageIndexType = typename OutputImageType::IndexType;
  using OutputImageValueType = typename OutputImagePixelType::ValueType;

  using TimeType = TimeProbe;

  using InputConstIterator = ImageRegionConstIterator<InputImageType>;
  using InputIterator = ImageRegionIterator<InputImageType>;
  using OutputIterator = ImageRegionIterator<OutputImageType>;

  using VectorWarperType = WarpVectorImageFilter<TOutputImage, TInputImage, TOutputImage>;

  using FieldInterpolatorType = VectorLinearInterpolateImageFunction<TInputImage, double>;
  using FieldInterpolatorPointer = typename FieldInterpolatorType::Pointer;
  using FieldInterpolatorOutputType = typename FieldInterpolatorType::OutputType;

  itkSetMacro(NumberOfIterations, unsigned int);
  itkGetConstMacro(NumberOfIterations, unsigned int);
  using ThreadRegionType = typename OutputImageType::RegionType;

  // If the error (in mm) between forward and backward mapping is smaller than the StopValue,
  // the algorithm stops.
  // This value can be used to speed up the calculation.
  itkSetMacro(StopValue, double);
  itkGetConstMacro(StopValue, double);

  char * GetReport()
  {
    return this->m_Report;
  }

  OutputImagePixelType TrilinearInterpolationFast(float& fDesiredX, float& fDesiredY, float& fDesiredZ,
                                                  InputImageSizeType size);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<OutputImageValueType> ) );
  /** End concept checking */
#endif
protected:
  ICCIterativeInverseDisplacementFieldImageFilter();
  ~ICCIterativeInverseDisplacementFieldImageFilter() override
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const override;

  void MakeReport();

  void GenerateData() override;

  void ComputeInverse(InputImageConstPointer &, OutputImagePointer &);

  void ThreadedComputeInverse(InputImageConstPointer &, OutputImagePointer &, const ThreadRegionType & regionToProcess,
                              int);

  InputImageIndexType  BoundaryIndexing(int, int, int, const int, const int, const int);

  unsigned int m_NumberOfIterations;
  double       m_StopValue;
  double       m_Time;

  struct ThreadStruct
    {
    ICCIterativeInverseDisplacementFieldImageFilter *Filter;
    InputImageConstPointer inputPtr;
    OutputImagePointer outputPtr;
    };
private:
  static ITK_THREAD_RETURN_TYPE ComputeInverseThreaderCallback(void * arg);
};
} // end namespace itk

#include "itkICCIterativeInverseDisplacementFieldImageFilter.hxx"

#endif
