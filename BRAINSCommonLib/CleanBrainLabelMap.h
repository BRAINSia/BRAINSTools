#ifndef CleanBrainLabelMap_h
#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkVotingBinaryHoleFillingImageFilter.h"

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer
CleanBrainLabelMap(const TInputImage *inputImage)
{
  typedef typename
    itk::BinaryThresholdImageFilter<TInputImage, TInputImage> BinaryThresholdFilterType;

  typename BinaryThresholdFilterType::Pointer binaryThresholdFilter =
    BinaryThresholdFilterType::New();
  binaryThresholdFilter->SetLowerThreshold(1);
  binaryThresholdFilter->SetUpperThreshold(255);
  binaryThresholdFilter->SetInput(inputImage);
  binaryThresholdFilter->Update();
  typename TInputImage::Pointer emsBrainMask(binaryThresholdFilter->GetOutput() );

  typedef typename itk::FlatStructuringElement<TInputImage::ImageDimension> KernelType;
  typename KernelType::RadiusType erodeRadius = { { 2, 2, 2 } };
  KernelType erodeKernel = KernelType::Ball(erodeRadius);

  typedef typename itk::BinaryErodeImageFilter<TInputImage, TInputImage, KernelType> BinaryErodeFilterType;
  typename BinaryErodeFilterType::Pointer erodeFilter = BinaryErodeFilterType::New();

  erodeFilter->SetInput(emsBrainMask);
  erodeFilter->SetKernel(erodeKernel);

  typedef typename itk::RelabelComponentImageFilter<TInputImage, TInputImage>
    RelabelComponentFilterType;
  typename RelabelComponentFilterType::Pointer relabelFilter =
    RelabelComponentFilterType::New();
  relabelFilter->SetInput(erodeFilter->GetOutput() );
  relabelFilter->SetMinimumObjectSize(30000);

  typename KernelType::RadiusType dilateRadius = { { 4, 4, 4 } };
  KernelType dilateKernel = KernelType::Ball(dilateRadius);

  typedef typename itk::BinaryDilateImageFilter<TInputImage, TInputImage, KernelType> BinaryDilateFilterType;
  typename BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();

  dilateFilter->SetKernel(dilateKernel);
  dilateFilter->SetInput(relabelFilter->GetOutput() );

  typedef typename itk::AndImageFilter<TInputImage, TInputImage, TInputImage> AndFilterType;
  typename AndFilterType::Pointer andFilter = AndFilterType::New();
  andFilter->SetInput1(emsBrainMask);
  andFilter->SetInput2(dilateFilter->GetOutput() );

  typename TInputImage::SizeType holeFillingRadius = { { 3, 3, 3 } };

  typedef typename itk::VotingBinaryHoleFillingImageFilter<TInputImage, TOutputImage> HoleFillingFilterType;
  typename HoleFillingFilterType::Pointer holeFillingFilter =
    HoleFillingFilterType::New();
  holeFillingFilter->SetInput(andFilter->GetOutput() );
  holeFillingFilter->SetRadius(holeFillingRadius);
  holeFillingFilter->SetForegroundValue(1);
  holeFillingFilter->SetBackgroundValue(0);
  holeFillingFilter->Update();
  return holeFillingFilter->GetOutput();
}

#endif // CleanBrainLabelMap_h
