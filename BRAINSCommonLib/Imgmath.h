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
#ifndef __Imgmath_h
#define __Imgmath_h

/* This file contains the functions to perform pixel by pixel mathematical
 * operations on 2 images. All the functions are performed by using ITK
 * filters. */

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMinimumImageFilter.h"
#include "itkMaximumImageFilter.h"
// INFO:  add these correctly so we get multi-threading.
#include "itkInvertIntensityImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSqrtImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkExpNegativeImageFilter.h"
#include "itkLogImageFilter.h"
#include <vcl_compiler.h>
#include <iostream>

/* Iadd adds 2 images at every pixel location and outputs the resulting image.*/
template <typename ImageType>
typename ImageType::Pointer
Iadd(const typename ImageType::Pointer input1, typename ImageType::Pointer input2)
{
  using FilterType = itk::AddImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Isub subtracts 2 images at every pixel location and outputs the resulting
 * image.*/
template <typename ImageType>
typename ImageType::Pointer
Isub(const typename ImageType::Pointer input1, const typename ImageType::Pointer input2)
{
  using FilterType = itk::SubtractImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Imul multiplies 2 images at every pixel location and outputs the resulting
 * image.*/
template <typename ImageType>
typename ImageType::Pointer
Imul(const typename ImageType::Pointer input1, const typename ImageType::Pointer input2)
{
  using FilterType = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Idiv divides 2 images at every pixel location and outputs the resulting
 * image.*/
template <typename ImageType>
typename ImageType::Pointer
Idiv(const typename ImageType::Pointer input1, const typename ImageType::Pointer input2)
{
  using FilterType = itk::DivideImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Imax does the numerical generalization of OR on 2 (non-negative) images at
 * every pixel location
 * and outputs the resulting image.*/
template <typename ImageType>
typename ImageType::Pointer
Imax(const typename ImageType::Pointer input1, const typename ImageType::Pointer input2)
{
  using FilterType = itk::MaximumImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Imin does the numerical generalization of AND on 2 (non-negative) images at
 * every pixel location
 * and outputs the resulting image.*/
template <typename ImageType>
typename ImageType::Pointer
Imin(const typename ImageType::Pointer input1, const typename ImageType::Pointer input2)
{
  using FilterType = itk::MinimumImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Iavg takes an image and the number of images as inputs , divides each pixel
  location of the image by the number of images outputs the resulting image.*/

template <typename ImageType>
typename ImageType::Pointer
Iavg(typename ImageType::Pointer input1, int nimgs)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(input1->GetLargestPossibleRegion());
  image->CopyInformation(input1);
  image->Allocate();
  using lConstIteratorType = typename itk::ImageRegionConstIterator<ImageType>;
  lConstIteratorType in1(input1, input1->GetLargestPossibleRegion());
  using IteratorType = typename itk::ImageRegionIterator<ImageType>;
  IteratorType out(image, image->GetLargestPossibleRegion());
  for (in1.GoToBegin(), out.GoToBegin(); !in1.IsAtEnd(); ++in1, ++out)
  {
    out.Set(in1.Get() / nimgs);
  }

  return image;
}

template <typename ImageType>
typename ImageType::Pointer
IMask(typename ImageType::Pointer input1, typename ImageType::Pointer mask)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(input1->GetLargestPossibleRegion());
  image->CopyInformation(input1);
  image->Allocate();
  using RegionIteratorType = typename itk::ImageRegionIterator<ImageType>;
  RegionIteratorType in1(input1, input1->GetLargestPossibleRegion());
  RegionIteratorType in2(mask, mask->GetLargestPossibleRegion());
  RegionIteratorType out(image, image->GetLargestPossibleRegion());
  for (in1.GoToBegin(), out.GoToBegin(), in2.GoToBegin(); !in1.IsAtEnd(); ++in1, ++in2, ++out)
  {
    const typename ImageType::PixelType temp = in1.Get();
    out.Set((in2.Get() > 0) ? temp : 0);
  }
  return image;
}

template <typename ImageType>
typename ImageType::Pointer
ImageAddConstant(const typename ImageType::Pointer input, const double shiftvalue)
{
  // INFO:  This should be a UnaryImageFunctor operation to get multi-threading.
  // KENT: Replace the use of this filter with itk::AddImageFilter
  typename ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions(input->GetLargestPossibleRegion());
  outImage->CopyInformation(input);
  outImage->Allocate();
  using RegionIteratorType = typename itk::ImageRegionIterator<ImageType>;
  RegionIteratorType in(input, input->GetLargestPossibleRegion());
  RegionIteratorType out(outImage, outImage->GetLargestPossibleRegion());
  out.GoToBegin();
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    out.Set(static_cast<typename ImageType::PixelType>((in.Get() + shiftvalue)));
    ++out;
  }
  return outImage;
}

template <typename ImageType>
typename ImageType::Pointer
ImageMultiplyConstant(const typename ImageType::Pointer input, const double scalevalue)
{
  using MultFilterType =
    typename itk::MultiplyImageFilter<ImageType, itk::Image<double, ImageType::ImageDimension>, ImageType>;
  typename MultFilterType::Pointer filt = MultFilterType::New();
  filt->SetInput(input);
  filt->SetConstant(scalevalue);
  filt->Update();
  return filt->GetOutput();
}

/* ImageComplementConstant does the numerical generalization of NOT on one
 * (non-negative) image
 * at every pixel location and outputs the resulting image.  For finding the
 *binary mask image Not,
 * referencevalue should be 1;  however there is no defaulting to 1 here. */
template <typename ImageType>
typename ImageType::Pointer
ImageComplementConstant(const typename ImageType::Pointer input, const double referencevalue)
{
  // INFO:  This should be a UnaryImageFunctor operation to get multi-threading.
  typename ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions(input->GetLargestPossibleRegion());
  outImage->CopyInformation(input);
  outImage->Allocate();
  using RegionIteratorType = typename itk::ImageRegionIterator<ImageType>;
  RegionIteratorType in(input, input->GetLargestPossibleRegion());
  RegionIteratorType out(outImage, outImage->GetLargestPossibleRegion());
  out.GoToBegin();
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    out.Set(static_cast<typename ImageType::PixelType>((referencevalue - in.Get())));
    ++out;
  }
  return outImage;
}

template <typename ImageType>
typename ImageType::Pointer
ImageDivideConstant(typename ImageType::Pointer input, const double denominator)
{
  typename ImageType::Pointer DivImage = ImageMultiplyConstant<ImageType>(input, 1.0 / denominator);
  return DivImage;
}

template <typename ImageType>
void
ImageSqrtValue(typename ImageType::Pointer Output, const typename ImageType::Pointer Input)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(Input->GetLargestPossibleRegion());
  image->CopyInformation(Input);
  image->Allocate();
  using RegionIteratorType = typename itk::ImageRegionIterator<ImageType>;
  RegionIteratorType in(Input, Input->GetLargestPossibleRegion());
  RegionIteratorType out(image, image->GetLargestPossibleRegion());
  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd(); ++in, ++out)
  {
    out.Set(static_cast<typename ImageType::PixelType>(in.Get()));
  }

  Output = image;
}

template <typename ImageType>
typename ImageType::Pointer
ImageSqrtValue(const typename ImageType::Pointer input)
{
  typename ImageType::Pointer rval;
  ImageSqrtValue<ImageType>(rval, input);
  return rval;
}

#endif
