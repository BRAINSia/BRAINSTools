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
#ifndef __itkFindeCenterOfBrainFilter_h
#define __itkFindeCenterOfBrainFilter_h
#include <itkImageToImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkLargestForegroundFilledMaskImageFilter.h"

namespace itk
{
/**
 * \class FindCenterOfBrainFilter
 */
template <typename TInputImage, typename TMaskImage = itk::Image<unsigned char, 3>>
class FindCenterOfBrainFilter : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  using Self = FindCenterOfBrainFilter;
  using Superclass = ImageToImageFilter<TInputImage, TInputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;
  itkNewMacro(Self);
  itkTypeMacro(FindCenterOfBrainFilter, Superclass);

  using ImageType = TInputImage;
  using MaskImageType = TMaskImage;
  using MaskImagePointer = typename MaskImageType::Pointer;
  using InputImagePointer = typename ImageType::Pointer;
  using PixelType = typename ImageType::PixelType;
  using PointType = typename ImageType::PointType;
  using SizeType = typename ImageType::SizeType;
  using SpacingType = typename ImageType::SpacingType;
  using IndexType = typename ImageType::IndexType;
  using ImageIteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using ImageConstIteratorType = typename itk::ImageRegionConstIteratorWithIndex<ImageType>;
  using LFFMaskFilterType = LargestForegroundFilledMaskImageFilter<ImageType, MaskImageType>;
  using DistanceImageType = typename itk::Image<float, 3>;
  using DistanceImagePointer = typename DistanceImageType::Pointer;
  /** Image related type alias. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  itkSetMacro(Maximize, bool);
  itkGetConstMacro(Maximize, bool);
  itkSetMacro(Axis, unsigned int);
  itkGetConstMacro(Axis, unsigned int);
  itkSetMacro(OtsuPercentileThreshold, double);
  itkGetConstMacro(OtsuPercentileThreshold, double);
  itkSetMacro(ClosingSize, unsigned int);
  itkGetConstMacro(ClosingSize, unsigned int);
  itkSetMacro(HeadSizeLimit, double);
  itkGetConstMacro(HeadSizeLimit, double);
  itkSetMacro(HeadSizeEstimate, double);
  itkGetConstMacro(HeadSizeEstimate, double);
  itkSetMacro(BackgroundValue, PixelType);
  itkGetConstMacro(BackgroundValue, PixelType);

  itkGetConstMacro(CenterOfBrain, PointType);
  itkGetModifiableObjectMacro(TrimmedImage, TInputImage);

  itkSetConstObjectMacro(ImageMask, TMaskImage);
  itkGetConstObjectMacro(ImageMask, TMaskImage);

  // THIS IS OUTPUT ONLY  itkSetObjectMacro(ClippedImageMask, TMaskImage);
  itkGetConstObjectMacro(ClippedImageMask, TMaskImage);

  // DEBUGGING STUFF
  itkSetMacro(GenerateDebugImages, bool);
  itkGetMacro(GenerateDebugImages, bool);
  DistanceImagePointer
  GetDebugDistanceImage() const
  {
    return m_DebugDistanceImage;
  }

  InputImagePointer
  GetDebugGridImage() const
  {
    return m_DebugGridImage;
  }

  MaskImagePointer
  GetDebugAfterGridComputationsForegroundImage() const
  {
    return m_DebugAfterGridComputationsForegroundImage;
  }

  MaskImagePointer
  GetDebugClippedImageMask() const
  {
    return m_DebugClippedImageMask;
  }

  InputImagePointer
  GetDebugTrimmedImage() const
  {
    return m_DebugTrimmedImage;
  }

protected:
  FindCenterOfBrainFilter();
  ~FindCenterOfBrainFilter() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  AllocateOutputs() override;

  void
  GenerateData() override;

private:
  bool         m_Maximize{ true };
  unsigned int m_Axis{ 2 };
  double       m_OtsuPercentileThreshold{ 0.001 };
  unsigned int m_ClosingSize{ 7 };
  double       m_HeadSizeLimit{ 1000 };
  double       m_HeadSizeEstimate{ 0 };
  PixelType    m_BackgroundValue;
  PointType    m_CenterOfBrain;
  //
  // DEBUGGING
  bool m_GenerateDebugImages{ false };
  /** The foreground mask, computed automatically if not specified
   * on the command line. **/
  typename TMaskImage::ConstPointer m_ImageMask;
  /** The foreground mask, computed automatically if
   * not specified on the command line. **/
  typename TMaskImage::Pointer  m_ClippedImageMask;
  typename TInputImage::Pointer m_TrimmedImage;
  DistanceImagePointer          m_DebugDistanceImage;
  InputImagePointer             m_DebugGridImage;
  MaskImagePointer              m_DebugAfterGridComputationsForegroundImage;
  MaskImagePointer              m_DebugClippedImageMask;
  InputImagePointer             m_DebugTrimmedImage;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkFindCenterOfBrainFilter.hxx"
#endif

#endif // itkFindeCenterOfBrainFilter_hxx
