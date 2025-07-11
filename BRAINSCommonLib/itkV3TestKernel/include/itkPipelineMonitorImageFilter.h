/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef __itkPipelineMonitorImageFilter_h
#define __itkPipelineMonitorImageFilter_h

#include "itkInPlaceImageFilter.h"

namespace itk
{
/** \class PipelineMonitorImageSource
 * \brief Enables monitoring, recording and debuging of the pipeline
 * execution and information exchange.
 *
 * This filter is useful for testing, debuging, and understanding the
 * pipeline. When DebugOn is enabled and compiled in Debug mode, many
 * itkDebug messages are printed. This filter also features, several
 * Verify methods which check the recorded information, for certain
 * conditions, which should occour when well behaved filters are
 * executed.
 *
 * There are two meta verify methods that should primarily be used
 * depending on the expected capabilities of the pipeline:
 * - \b VerifyAllInputCanStream
 * - \b VerifyAllInputCanNotStream
 *
 * During the pipeline execution this filter records a variety of
 * information to aid if verifying correct pipeline behavior:
 * - \b NumberOfUpdate the number of times GenerateData was executed
 * - \b OutputRequestedRegions an array of the output of this filter's
 *      requested region
 * - \b InputRequestedRegions an array of the input image's requested
 *      region after PropagateRequestedRegion
 * - \b UpdatedBufferedRegions an array of the input image's buffered
 *      region after upstream GenerateData
 * - \b UpdatedRequestedRegions an array of the input image's
 *      requested region after upstream GenerateData
 *
 * The following are recorded from the input image after the input's
 * output information is generated:
 * - \b UpdatedOutputOrigin
 * - \b UpdatedOutputDirection
 * - \b UpdatedOutputSpacing
 * - \b UpdatedOutputLargestPossibleRegion
 *
 * This filter always runs in-place so it has no per-pixel overhead.
 *
 * \ingroup ITKTestKernel
 */
template <typename TImageType>
class PipelineMonitorImageFilter : public InPlaceImageFilter<TImageType, TImageType>
{
public:
  using Self = PipelineMonitorImageFilter;
  using Superclass = InPlaceImageFilter<TImageType, TImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using PointType = typename TImageType::PointType;
  using DirectionType = typename TImageType::DirectionType;
  using SpacingType = typename TImageType::SpacingType;
  using InputImagePointer = typename TImageType::Pointer;
  using InputImageConstPointer = typename TImageType::ConstPointer;
  using ImageRegionType = typename Superclass::InputImageRegionType;

  using RegionVectorType = std::vector<typename TImageType::RegionType>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(PipelineMonitorImageFilter);

  /** Enable/Disable clearing all saved pipeline information when
   * GenerateOutputInformation is called.
   *
   * The NumberOfClearPipelines is incremented, to aid in detection
   * of certain pipeline errors caused but excessive execution of
   * GenerateOutputInformation.
   *
   * Defaults to On
   */
  itkSetMacro(ClearPipelineOnGenerateOutputInformation, bool);
  itkGetMacro(ClearPipelineOnGenerateOutputInformation, bool);
  itkBooleanMacro(ClearPipelineOnGenerateOutputInformation);

  /** This a meta verify method to check expected pipeline execution
   * when the pipeline is capable of streaming. See
   * VerifyInputFilterExecutedStreaming for information on the
   * expectedNumber parameter.
   */
  bool
  VerifyAllInputCanStream(int expectedNumber);

  /** Checks that the input filter didn't stream, and just updated
   * the largest possible region along with other correct behaviors.
   */
  bool
  VerifyAllInputCanNotStream(void);

  /** This method verifies that propagation was executed yet no
   * updating was needed.
   */
  bool
  VerifyAllNoUpdate(void);

  bool
  VerifyDownStreamFilterExecutedPropagation(void);

  /** Verifies the the GenerateData executed the expected number of
   * times.
   *
   * If expecetedNumber is positive then the number of updates must
   * match. If expectedNumber is negative then the number of updates
   * must at least be |expectedNumber|. If expectedNumber is zero,
   * then this method always returns true, and no verification is
   * performed.
   */
  bool
  VerifyInputFilterExecutedStreaming(int expectedNumber);

  /** Verifies that the output information didn't change between the
   * GenerateOutputInformation and the UpdateData phases of the
   * pipeline.
   */
  bool
  VerifyInputFilterMatchedUpdateOutputInformation(void);

  /** Verifies that the input filter buffered the requested region */
  bool
  VerifyInputFilterBufferedRequestedRegions(void);

  bool
  VerifyInputFilterMatchedRequestedRegions(void);

  bool
  VerifyInputFilterRequestedLargestRegion(void);

  unsigned int
  GetNumberOfUpdates(void) const
  {
    return m_NumberOfUpdates;
  }

  RegionVectorType
  GetOutputRequestedRegions(void) const
  {
    return m_OutputRequestedRegions;
  }

  RegionVectorType
  GetInputRequestedRegions(void) const
  {
    return m_InputRequestedRegions;
  }

  RegionVectorType
  GetUpdatedBufferedRegions(void) const
  {
    return m_UpdatedBufferedRegions;
  }

  RegionVectorType
  GetUpdatedRequestedRegions(void) const
  {
    return m_UpdatedRequestedRegions;
  }

  /** Clears all saved pipeline information, but increments
   * NumberOfClearPipeline. */
  void
  ClearPipelineSavedInformation(void);

  /** Standard pipeline methods are overloaded to call superclass's
   * implementation and record information.
   */
  virtual void
  GenerateOutputInformation(void);

  virtual void
  PropagateRequestedRegion(DataObject * output);

  virtual void
  EnlargeOutputRequestedRegion(DataObject * output);

  virtual void
  GenerateInputRequestedRegion(void);

  virtual void
  GenerateData(void);

protected:
  PipelineMonitorImageFilter(void);

  // ~PipelineMonitorImageFilter() { } default implementation OK

  void
  PrintSelf(std::ostream & os, Indent indent) const;

private:
  PipelineMonitorImageFilter(const PipelineMonitorImageFilter &); // not implemented
  void
  operator=(const PipelineMonitorImageFilter &); // not implemented

  bool m_ClearPipelineOnGenerateOutputInformation;

  unsigned int m_NumberOfUpdates;

  unsigned int m_NumberOfClearPipeline;

  RegionVectorType m_OutputRequestedRegions;
  RegionVectorType m_InputRequestedRegions;
  RegionVectorType m_UpdatedBufferedRegions;
  RegionVectorType m_UpdatedRequestedRegions;

  PointType       m_UpdatedOutputOrigin;
  DirectionType   m_UpdatedOutputDirection;
  SpacingType     m_UpdatedOutputSpacing;
  ImageRegionType m_UpdatedOutputLargestPossibleRegion;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkPipelineMonitorImageFilter.hxx"
#endif

#endif // __itkPipelineMonitorImageFilter_hxx
