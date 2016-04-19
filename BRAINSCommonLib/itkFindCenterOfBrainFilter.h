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
template <class TInputImage, class TMaskImage = itk::Image<unsigned char, 3> >
class FindCenterOfBrainFilter :
  public         ImageToImageFilter<TInputImage, TInputImage>
{
public:
  typedef FindCenterOfBrainFilter                      Self;
  typedef ImageToImageFilter<TInputImage, TInputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
  itkNewMacro(Self);
  itkTypeMacro(FindCenterOfBrain, Superclass);

  typedef TInputImage                     ImageType;
  typedef TMaskImage                      MaskImageType;
  typedef typename MaskImageType::Pointer MaskImagePointer;
  typedef typename ImageType::Pointer     InputImagePointer;
  typedef typename ImageType::PixelType   PixelType;
  typedef typename ImageType::PointType   PointType;
  typedef typename ImageType::SizeType    SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::IndexType   IndexType;
  typedef typename itk::ImageRegionIteratorWithIndex<ImageType>
    ImageIteratorType;
  typedef typename itk::ImageRegionConstIteratorWithIndex<ImageType>
    ImageConstIteratorType;
  typedef LargestForegroundFilledMaskImageFilter<ImageType, MaskImageType>
    LFFMaskFilterType;
  typedef typename itk::Image<float, 3>       DistanceImageType;
  typedef typename DistanceImageType::Pointer DistanceImagePointer;
  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

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
  DistanceImagePointer GetDebugDistanceImage() const
  {
    return m_DebugDistanceImage;
  }

  InputImagePointer GetDebugGridImage() const
  {
    return m_DebugGridImage;
  }

  MaskImagePointer GetDebugAfterGridComputationsForegroundImage() const
  {
    return m_DebugAfterGridComputationsForegroundImage;
  }

  MaskImagePointer GetDebugClippedImageMask() const
  {
    return m_DebugClippedImageMask;
  }

  InputImagePointer GetDebugTrimmedImage() const
  {
    return m_DebugTrimmedImage;
  }

protected:
  FindCenterOfBrainFilter();
  ~FindCenterOfBrainFilter();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void AllocateOutputs() ITK_OVERRIDE;

  virtual void GenerateData() ITK_OVERRIDE;

private:
  bool         m_Maximize;
  unsigned int m_Axis;
  double       m_OtsuPercentileThreshold;
  unsigned int m_ClosingSize;
  double       m_HeadSizeLimit;
  double       m_HeadSizeEstimate;
  PixelType    m_BackgroundValue;
  PointType    m_CenterOfBrain;
  //
  // DEBUGGING
  bool m_GenerateDebugImages;
  /** The foreground mask, computed automatically if not specified
   * on the command line. **/
  typename TMaskImage::ConstPointer m_ImageMask;
  /** The foreground mask, computed automatically if
   * not specified on the command line. **/
  typename TMaskImage::Pointer m_ClippedImageMask;
  typename TInputImage::Pointer m_TrimmedImage;
  DistanceImagePointer m_DebugDistanceImage;
  InputImagePointer    m_DebugGridImage;
  MaskImagePointer     m_DebugAfterGridComputationsForegroundImage;
  MaskImagePointer     m_DebugClippedImageMask;
  InputImagePointer    m_DebugTrimmedImage;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFindCenterOfBrainFilter.hxx"
#endif

#endif // itkFindeCenterOfBrainFilter_hxx
