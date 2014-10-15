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
/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile$
 *  Language:  C++
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/

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
// TODO:  add these correctly so we get multi-threading.
#include "itkInvertIntensityImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSqrtImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkExpNegativeImageFilter.h"
#include "itkLogImageFilter.h"
#include <vcl_cmath.h>

/* Iadd adds 2 images at every pixel location and outputs the resulting image.*/
template <class ImageType>
typename ImageType::Pointer Iadd(const typename ImageType::Pointer input1,
                                 typename ImageType::Pointer input2)
{
  typedef itk::AddImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Isub subtracts 2 images at every pixel location and outputs the resulting
  * image.*/
template <class ImageType>
typename ImageType::Pointer Isub(const typename ImageType::Pointer input1,
                                 const typename ImageType::Pointer input2)
{
  typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Imul multiplies 2 images at every pixel location and outputs the resulting
  * image.*/
template <class ImageType>
typename ImageType::Pointer Imul(const typename ImageType::Pointer input1,
                                 const typename ImageType::Pointer input2)
{
  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Idiv divides 2 images at every pixel location and outputs the resulting
  * image.*/
template <class ImageType>
typename ImageType::Pointer Idiv(const typename ImageType::Pointer input1,
                                 const typename ImageType::Pointer input2)
{
  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> FilterType;
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
template <class ImageType>
typename ImageType::Pointer Imax(const typename ImageType::Pointer input1,
                                 const typename ImageType::Pointer input2)
{
  typedef itk::MaximumImageFilter<ImageType, ImageType, ImageType> FilterType;
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
template <class ImageType>
typename ImageType::Pointer Imin(const typename ImageType::Pointer input1,
                                 const typename ImageType::Pointer input2)
{
  typedef itk::MinimumImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(input1);
  filter->SetInput2(input2);
  filter->Update();
  typename ImageType::Pointer image = filter->GetOutput();
  return image;
}

/* Iavg takes an image and the number of images as inputs , divides each pixel
  location of the image by the number of images outputs the resulting image.*/

template <class ImageType>
typename ImageType::Pointer Iavg(typename ImageType::Pointer input1, int nimgs)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions( input1->GetLargestPossibleRegion() );
  image->CopyInformation(input1);
  image->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> ConstIteratorType;
  ConstIteratorType in1( input1, input1->GetLargestPossibleRegion() );
  ConstIteratorType out( image, image->GetLargestPossibleRegion() );
  for( in1.GoToBegin(), out.GoToBegin(); !in1.IsAtEnd(); ++in1, ++out )
    {
    out.Set(in1.Get() / nimgs);
    }

  return image;
}

template <class ImageType>
typename ImageType::Pointer IMask(typename ImageType::Pointer input1,
                                  typename ImageType::Pointer mask)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions( input1->GetLargestPossibleRegion() );
  image->CopyInformation(input1);
  image->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> RegionIteratorType;
  RegionIteratorType in1( input1, input1->GetLargestPossibleRegion() );
  RegionIteratorType in2( mask, mask->GetLargestPossibleRegion() );
  RegionIteratorType out( image, image->GetLargestPossibleRegion() );
  for( in1.GoToBegin(), out.GoToBegin(), in2.GoToBegin();
       !in1.IsAtEnd();
       ++in1, ++in2, ++out )
    {
    const typename ImageType::PixelType temp = in1.Get();
    out.Set( ( in2.Get() > 0 ) ? temp : 0 );
    }
  return image;
}

template <class ImageType>
typename ImageType::Pointer ImageAddConstant(
  const typename ImageType::Pointer input,
  const double shiftvalue)
{
  // TODO:  This should be a UnaryImageFunctor operation to get multi-threading.
  // KENT: Replace the use of this filter with itk::AddImageFilter
  typename ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( input->GetLargestPossibleRegion() );
  outImage->CopyInformation(input);
  outImage->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> RegionIteratorType;
  RegionIteratorType in( input, input->GetLargestPossibleRegion() );
  RegionIteratorType out( outImage, outImage->GetLargestPossibleRegion() );
  out.GoToBegin();
  for( in.GoToBegin(); !in.IsAtEnd(); ++in )
    {
    out.Set( static_cast<typename ImageType::PixelType>( ( in.Get()
                                                           + shiftvalue ) ) );
    ++out;
    }
  return outImage;
}

template <class ImageType>
typename ImageType::Pointer ImageMultiplyConstant(
  const typename ImageType::Pointer input,
  const double scalevalue)
{
  typedef typename itk::MultiplyImageFilter<ImageType, itk::Image<double, 3>, ImageType>
    MultFilterType;
  typedef typename MultFilterType::Pointer MultFilterPointer;
  MultFilterPointer filt = MultFilterType::New();
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
template <class ImageType>
typename ImageType::Pointer ImageComplementConstant(
  const typename ImageType::Pointer input,
  const double referencevalue)
{
  // TODO:  This should be a UnaryImageFunctor operation to get multi-threading.
  typename ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions( input->GetLargestPossibleRegion() );
  outImage->CopyInformation(input);
  outImage->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> RegionIteratorType;
  RegionIteratorType in( input, input->GetLargestPossibleRegion() );
  RegionIteratorType out( outImage, outImage->GetLargestPossibleRegion() );
  out.GoToBegin();
  for( in.GoToBegin(); !in.IsAtEnd(); ++in )
    {
    out.Set( static_cast<typename ImageType::PixelType>( ( referencevalue
                                                           - in.Get() ) ) );
    ++out;
    }
  return outImage;
}

template <class ImageType>
typename ImageType::Pointer ImageDivideConstant(
  typename ImageType::Pointer input,
  const double denominator)
{
  typename ImageType::Pointer DivImage = ImageMultiplyConstant<ImageType>(input,
                                                                          1.0 / denominator);
  return DivImage;
}

template <class ImageType>
void ImageSqrtValue(typename ImageType::Pointer Output,
                    const typename ImageType::Pointer Input)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions( Input->GetLargestPossibleRegion() );
  image->CopyInformation(Input);
  image->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> RegionIteratorType;
  RegionIteratorType in( Input, Input->GetLargestPossibleRegion() );
  RegionIteratorType out( image, image->GetLargestPossibleRegion() );
  for( in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd(); ++in, ++out )
    {
    out.Set( static_cast<typename ImageType::PixelType>( in.Get() ) );
    }

  Output = image;
  return;
}

template <class ImageType>
typename ImageType::Pointer
ImageSqrtValue(const typename ImageType::Pointer input)
{
  typename ImageType::Pointer rval;
  ImageSqrtValue<ImageType>(rval, input);
  return rval;
}

#endif
