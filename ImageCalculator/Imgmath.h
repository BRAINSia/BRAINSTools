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
#ifndef _ITKMATH_H_
#define _ITKMATH_H_

/* This file contains the functions to perform pixel by pixel mathematical
 * operations on 2 images. All the functions are performed by using ITK
 * filters. */

#include "itkImage.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include <vcl_compiler.h>
#include <iostream>
#include <cmath>

/* Iadd adds 2 images at every pixel location and outputs the resulting image.*/

template <typename ImageType>
typename ImageType::Pointer Iadd(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  using FilterType = itk::AddImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Isub subtracts 2 images at every pixel location and outputs the resulting image.*/

template <typename ImageType>
typename ImageType::Pointer Isub(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  using FilterType = itk::SubtractImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Imul multiplies 2 images at every pixel location and outputs the resulting image.*/

template <typename ImageType>
typename ImageType::Pointer Imul(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  using FilterType = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Idiv divides 2 images at every pixel location and outputs the resulting image.*/

template <typename ImageType>
typename ImageType::Pointer Idiv(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  using FilterType = itk::DivideImageFilter<ImageType, ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Iavg takes an image and the number of images as inputs , divides each pixel location of the image by the number of
  images outputs the resulting image.*/

template <typename ImageType>
typename ImageType::Pointer Iavg(typename ImageType::Pointer input1, int nimgs)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(input1->GetLargestPossibleRegion() );
  image->CopyInformation(input1);
  image->Allocate();
  using ConstIteratorType = typename itk::ImageRegionIterator<ImageType>;
  ConstIteratorType in1(input1, input1->GetLargestPossibleRegion() );
  ConstIteratorType out(image, image->GetLargestPossibleRegion() );
  for( in1.GoToBegin(), out.GoToBegin(); !in1.IsAtEnd(); ++in1, ++out )
    {
    out.Set(in1.Get() / nimgs);
    }

  return image;
}

template <typename ImageType>
typename ImageType::Pointer IMask(typename ImageType::Pointer input1, typename ImageType::Pointer mask )

{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(input1->GetLargestPossibleRegion() );
  image->CopyInformation(input1);
  image->Allocate();
  using ConstIteratorType = typename itk::ImageRegionIterator<ImageType>;
  ConstIteratorType in1(input1, input1->GetLargestPossibleRegion() );
  ConstIteratorType in2(mask, mask->GetLargestPossibleRegion() );
  ConstIteratorType out(image, image->GetLargestPossibleRegion() );
  for( in1.GoToBegin(), out.GoToBegin(), in2.GoToBegin(); !in1.IsAtEnd(); ++in1, ++in2, ++out )
    {
    const typename ImageType::PixelType temp = in1.Get();
    out.Set( (in2.Get() > 0) ? temp : 0);
    }
  return image;
}

/*ImageMultiplyConstant multiplies the entire image with a constant value and outputs the resultant image*/

template <typename ImageType>
typename ImageType::Pointer ImageMultiplyConstant(typename ImageType::Pointer input1,
                                                  typename ImageType::PixelType constant)
{
  using ConstIteratorType = typename itk::ImageRegionIterator<ImageType>;
  ConstIteratorType in1(input1, input1->GetLargestPossibleRegion() );
  for( in1.GoToBegin(); !in1.IsAtEnd(); ++in1 )
    {
    in1.Set( (in1.Get() * constant) );
    }

  return input1;
}

template <typename ImageType>
typename ImageType::Pointer ImageDivideConstant(typename ImageType::Pointer input1,
                                                typename ImageType::PixelType constant)
{
  using ConstIteratorType = typename itk::ImageRegionIterator<ImageType>;
  ConstIteratorType in1(input1, input1->GetRequestedRegion() );
  for( in1.GoToBegin(); !in1.IsAtEnd(); ++in1 )
    {
    in1.Set( (in1.Get() / constant) );
    }

  return input1;
}

template <typename ImageType>
void ImageSqrtValue(typename ImageType::Pointer Output,
                    const typename ImageType::Pointer Input)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(Input->GetLargestPossibleRegion() );
  image->CopyInformation(Input);
  image->Allocate();
  using ConstIteratorType = typename itk::ImageRegionIterator<ImageType>;
  ConstIteratorType in(Input, Input->GetLargestPossibleRegion() );
  ConstIteratorType out(image, image->GetLargestPossibleRegion() );
  for( in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd(); ++in, ++out )
    {
    out.Set(static_cast<typename ImageType::PixelType>(in.Get() ) );
    }

  Output = image;
  return;
}

template <typename ImageType>
typename ImageType::Pointer
ImageSqrtValue(typename ImageType::Pointer input)
{
  typename ImageType::Pointer rval;
  ImageSqrtValue<ImageType>(rval, input);
  return rval;
}

#endif
