#include <iostream>
#include "itkIO.h"
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>

template <typename ImageType>
int Resample(const std::string & inputVolume,
             const std::string & outputVolume)
{
  typename ImageType::Pointer inputImage;
  try
    {
    inputImage = itkUtil::ReadImage<ImageType>(inputVolume);
    }
  catch( ... )
    {
    std::cout << "Caught an exception: " << std::endl;
    std::cout << " " << __FILE__ << " " << __LINE__ << std::endl;
    exit(1);
    }

  typename ImageType::RegionType region =
    inputImage->GetLargestPossibleRegion();
  typename ImageType::SizeType    size(region.GetSize() );
  typename ImageType::SpacingType spacing(inputImage->GetSpacing() );

  typedef typename
    itk::ResampleImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter(FilterType::New() );

  typedef typename
    itk::IdentityTransform<double, ImageType::ImageDimension>  TransformType;
  typename TransformType::Pointer transform(TransformType::New() );
  transform->SetIdentity();
  filter->SetTransform(transform);

  typedef typename itk::LinearInterpolateImageFunction<
      ImageType, double>  InterpolatorType;
  typename InterpolatorType::Pointer interpolator(InterpolatorType::New() );
  filter->SetInterpolator(interpolator);
  for( unsigned i = 0; i < 3; i++ )
    {
    spacing[i] *= 4.0;
    size[i] = static_cast<typename ImageType::SizeType::SizeValueType>
      (static_cast<double>(size[i]) / 4.0);
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
    }

  try
    {
    itkUtil::WriteImage<ImageType>(outputImage, outputVolume);
    }
  catch( ... )
    {
    std::cout << "Caught an exception: " << std::endl;
    std::cout << " " << __FILE__ << " " << __LINE__ << std::endl;
    exit(1);
    }
  exit(0);
}

int
main(int argc, char * *argv)
{
  if( argc < 4 )
    {
    std::cerr << "resamp: "
              << "Usage resamp pixType <inputImage> <outputImage"
              << std::endl;
    exit(1);
    }
  std::string pixType(argv[1]);
  std::string inputVolume(argv[2]);
  std::string outputVolume(argv[3]);
  if( pixType == "short" )
    {
    Resample<itk::Image<short, 3> >(inputVolume, outputVolume);
    }
  else if( pixType == "uint" )
    {
    Resample<itk::Image<unsigned int, 3> >(inputVolume, outputVolume);
    }
  else if( pixType == "float" )
    {
    Resample<itk::Image<float, 3> >(inputVolume, outputVolume);
    }
  else if( pixType == "uchar" )
    {
    Resample<itk::Image<unsigned char, 3> >(inputVolume, outputVolume);
    }
}
