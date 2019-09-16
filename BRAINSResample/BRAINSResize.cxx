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
#include <iostream>
#include "itkIO.h"
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include "BRAINSResizeCLP.h"
#include <BRAINSCommonLib.h>

template <typename ImageType>
int
Resample(const std::string & inputVolume, const std::string & outputVolume, const double scaleFactor)
{
  typename ImageType::Pointer inputImage;
  inputImage = itkUtil::ReadImage<ImageType>(inputVolume);

  typename ImageType::RegionType  region = inputImage->GetLargestPossibleRegion();
  typename ImageType::SizeType    size(region.GetSize());
  typename ImageType::SpacingType spacing(inputImage->GetSpacing());

  using FilterType = typename itk::ResampleImageFilter<ImageType, ImageType>;
  typename FilterType::Pointer filter(FilterType::New());

  using TransformType = typename itk::IdentityTransform<double, ImageType::ImageDimension>;
  typename TransformType::Pointer transform(TransformType::New());
  transform->SetIdentity();
  filter->SetTransform(transform);

  using InterpolatorType = typename itk::LinearInterpolateImageFunction<ImageType, double>;
  typename InterpolatorType::Pointer interpolator(InterpolatorType::New());
  filter->SetInterpolator(interpolator);
  for (unsigned i = 0; i < 3; i++)
  {
    spacing[i] *= scaleFactor;
    size[i] = static_cast<typename ImageType::SizeType::SizeValueType>(static_cast<double>(size[i]) / scaleFactor);
  }

  filter->SetOutputSpacing(spacing);
  filter->SetOutputOrigin(inputImage->GetOrigin());
  filter->SetOutputDirection(inputImage->GetDirection());
  filter->SetSize(size);
  filter->SetInput(inputImage);

  typename ImageType::Pointer outputImage;
  try
  {
    filter->Update();
    outputImage = filter->GetOutput();
  }
  catch (itk::ExceptionObject & excep)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  itkUtil::WriteImage<ImageType>(outputImage, outputVolume);
  return EXIT_SUCCESS;
}

int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  try
  {
    if (pixelType == "short")
    {
      return Resample<itk::Image<short, 3>>(inputVolume, outputVolume, scaleFactor);
    }
    else if (pixelType == "uint")
    {
      return Resample<itk::Image<unsigned int, 3>>(inputVolume, outputVolume, scaleFactor);
    }
    else if (pixelType == "float")
    {
      return Resample<itk::Image<float, 3>>(inputVolume, outputVolume, scaleFactor);
    }
    else if (pixelType == "uchar")
    {
      return Resample<itk::Image<unsigned char, 3>>(inputVolume, outputVolume, scaleFactor);
    }
  }
  catch (itk::ExceptionObject & excep)
  {
    std::cerr << "BRAINSResize " << excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_FAILURE;
}
