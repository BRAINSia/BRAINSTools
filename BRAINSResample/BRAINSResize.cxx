#include <iostream>
#include "itkIO.h"
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkExceptionObject.h>
#include "BRAINSResizeCLP.h"

template <typename ImageType>
int Resample(const std::string & inputVolume, const std::string & outputVolume, const double scaleFactor)
{
  typename ImageType::Pointer inputImage;
  inputImage = itkUtil::ReadImage<ImageType>(inputVolume);

  typename ImageType::RegionType region = inputImage->GetLargestPossibleRegion();
  typename ImageType::SizeType    size(region.GetSize() );
  typename ImageType::SpacingType spacing(inputImage->GetSpacing() );

  typedef typename itk::ResampleImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter(FilterType::New() );

  typedef typename itk::IdentityTransform<double, ImageType::ImageDimension> TransformType;
  typename TransformType::Pointer transform(TransformType::New() );
  transform->SetIdentity();
  filter->SetTransform(transform);

  typedef typename itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator(InterpolatorType::New() );
  filter->SetInterpolator(interpolator);
  for( unsigned i = 0; i < 3; i++ )
    {
    spacing[i] *= scaleFactor;
    size[i] = static_cast<typename ImageType::SizeType::SizeValueType>
      (static_cast<double>(size[i]) / scaleFactor);
    }

  filter->SetOutputSpacing(spacing);
  filter->SetOutputOrigin(inputImage->GetOrigin() );
  filter->SetOutputDirection(inputImage->GetDirection() );
  filter->SetSize(size);
  filter->SetInput(inputImage);

  typename ImageType::Pointer outputImage;
  try
    {
    filter->Update();
    outputImage = filter->GetOutput();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  itkUtil::WriteImage<ImageType>(outputImage, outputVolume);
  return EXIT_SUCCESS;
}

int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  try
    {
    if( pixelType == "short" )
      {
      return Resample<itk::Image<short, 3> >(inputVolume, outputVolume, scaleFactor);
      }
    else if( pixelType == "uint" )
      {
      return Resample<itk::Image<unsigned int, 3> >(inputVolume, outputVolume, scaleFactor);
      }
    else if( pixelType == "float" )
      {
      return Resample<itk::Image<float, 3> >(inputVolume, outputVolume, scaleFactor);
      }
    else if( pixelType == "uchar" )
      {
      return Resample<itk::Image<unsigned char, 3> >(inputVolume, outputVolume, scaleFactor);
      }
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "BRAINSResize " << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_FAILURE;
}
