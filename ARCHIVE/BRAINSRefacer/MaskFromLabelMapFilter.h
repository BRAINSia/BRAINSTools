#ifndef BRAINSTOOLS_MASKFROMLABELMAPFILTER_H
#define BRAINSTOOLS_MASKFROMLABELMAPFILTER_H


#include "itkImageToImageFilter.h"
#include "itkMacro.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


template <typename TImage, typename TAtlas, typename TOutputMask>
class MaskFromLabelMapFilter : public itk::ImageToImageFilter<TImage, TOutputMask>
{
public:
  /** Standard class type alias. */
  using Self = MaskFromLabelMapFilter;
  using Superclass = itk::ImageToImageFilter<TImage, TImage>;
  using Pointer = itk::SmartPointer<Self>;

  itkNewMacro(Self);
  itkTypeMacro(ImageFilterMultipleInputsDifferentType, ImageToImageFilter);

  itkSetMacro(ReferenceImage, typename TImage::Pointer);
  itkGetMacro(ReferenceImage, typename TImage::Pointer);

  itkSetMacro(InputAtlas, typename TAtlas::Pointer);
  itkGetMacro(InputAtlas, typename TAtlas::Pointer);

private:
  typename TImage::Pointer m_ReferenceImage;
  typename TAtlas::Pointer m_InputAtlas;

protected:
  MaskFromLabelMapFilter() { this->SetInput(TImage::New()); };
  ~MaskFromLabelMapFilter() override = default;
  ;

  void
  GenerateData() override
  {
    // resample LabelImage
    using NN_InterpolatorType = itk::NearestNeighborInterpolateImageFunction<TAtlas, double>;
    typename NN_InterpolatorType::Pointer NN_interpolator = NN_InterpolatorType::New();

    constexpr int Dimension = 3;
    using IdentityTransformType = itk::IdentityTransform<double, Dimension>;
    typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

    using maskResamplerType = itk::ResampleImageFilter<TAtlas, TAtlas>;
    typename maskResamplerType::Pointer maskResampler = maskResamplerType::New();

    std::cout << "Resampling atlas map:" << std::endl;
    maskResampler->SetInput(this->GetInputAtlas());
    maskResampler->SetInterpolator(NN_interpolator);
    maskResampler->SetTransform(identityTransform);
    maskResampler->SetReferenceImage(this->GetReferenceImage());
    maskResampler->UseReferenceImageOn();
    maskResampler->Update();

    using MaskFilterType = itk::BinaryThresholdImageFilter<TAtlas, TOutputMask>;
    typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();

    maskFilter->SetInput(maskResampler->GetOutput());
    maskFilter->SetOutsideValue(1);
    maskFilter->SetInsideValue(0);
    maskFilter->SetLowerThreshold(0);
    maskFilter->SetUpperThreshold(0);

    maskFilter->Update();
    typename TOutputMask::Pointer maskResult = maskFilter->GetOutput();
    typename TOutputMask::Pointer output = this->GetOutput();
    output->SetRegions(maskResult->GetLargestPossibleRegion());
    itk::Image<double, 3>::Pointer test = itk::Image<double, 3>::New();
    output->SetOrigin(maskResult->GetOrigin());
    output->SetSpacing(maskResult->GetSpacing());
    output->SetDirection(maskResult->GetDirection());
    output->Allocate();

    typename itk::ImageRegionConstIterator<TOutputMask> maskResultIterator(maskResult,
                                                                           maskResult->GetLargestPossibleRegion());

    typename itk::ImageRegionIterator<TOutputMask> outputIterator(output, output->GetLargestPossibleRegion());

    while (!outputIterator.IsAtEnd())
    {
      outputIterator.Set(maskResultIterator.Get());
      ++outputIterator;
      ++maskResultIterator;
    }
  };

private:
  MaskFromLabelMapFilter(const Self &) = delete; // purposely not implemented
  void
  operator=(const Self &) = delete; // purposely not implemented
};
#endif // BRAINSTOOLS_MASKFROMLABELMAPFILTER_H
