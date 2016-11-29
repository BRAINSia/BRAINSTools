#ifndef BRAINSTOOLS_MASKFROMLABELMAPFILTER_H
#define BRAINSTOOLS_MASKFROMLABELMAPFILTER_H


#include "itkImageToImageFilter.h"
#include "itkMacro.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


template< typename TImage, typename TAtlas, typename TOutputMask>
class MaskFromLabelMapFilter : public itk::ImageToImageFilter< TImage, TOutputMask >
  {
public:
  /** Standard class typedefs. */
  typedef MaskFromLabelMapFilter             Self;
  typedef itk::ImageToImageFilter< TImage, TImage > Superclass;
  typedef itk::SmartPointer< Self >        Pointer;

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
  MaskFromLabelMapFilter()
    {
    this->SetInput(TImage::New());
    };
  ~MaskFromLabelMapFilter(){};

  virtual void GenerateData() ITK_OVERRIDE
    {
    //resample LabelImage
    typedef itk::NearestNeighborInterpolateImageFunction<TAtlas, double> NN_InterpolatorType;
    typename NN_InterpolatorType::Pointer NN_interpolator = NN_InterpolatorType::New();

    const int Dimension = 3;
    typedef itk::IdentityTransform<double, Dimension> IdentityTransformType;
    typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

    typedef itk::ResampleImageFilter<TAtlas, TAtlas> maskResamplerType;
    typename maskResamplerType::Pointer maskResampler = maskResamplerType::New();

    std::cout << "Resampling atlas map:" << std::endl;
    maskResampler->SetInput(this->GetInputAtlas());
    maskResampler->SetInterpolator(NN_interpolator);
    maskResampler->SetTransform(identityTransform);
    maskResampler->SetReferenceImage(this->GetReferenceImage());
    maskResampler->UseReferenceImageOn();
    maskResampler->Update();

    typedef itk::BinaryThresholdImageFilter<TAtlas, TOutputMask> MaskFilterType;
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
    itk::Image<double,3>::Pointer test = itk::Image<double, 3>::New();
    output->SetOrigin(maskResult->GetOrigin());
    output->SetSpacing(maskResult->GetSpacing());
    output->SetDirection(maskResult->GetDirection());
    output->Allocate();

    typename itk::ImageRegionConstIterator<TOutputMask> maskResultIterator(maskResult, maskResult->GetLargestPossibleRegion());

    typename itk::ImageRegionIterator<TOutputMask> outputIterator(output, output->GetLargestPossibleRegion());

    while(!outputIterator.IsAtEnd())
      {
      outputIterator.Set(maskResultIterator.Get());
      ++outputIterator;
      ++maskResultIterator;
      }
    };

private:
  MaskFromLabelMapFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  };
#endif //BRAINSTOOLS_MASKFROMLABELMAPFILTER_H