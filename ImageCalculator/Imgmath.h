/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    Imgmath.h
  Language:  C++

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

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
#include <vcl_cmath.h>

/* Iadd adds 2 images at every pixel location and outputs the resulting image.*/

template <class ImageType>
typename ImageType::Pointer Iadd(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  typedef itk::AddImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    exit(-1);
    }
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Isub subtracts 2 images at every pixel location and outputs the resulting image.*/

template <class ImageType>
typename ImageType::Pointer Isub(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    exit(-1);
    }
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Imul multiplies 2 images at every pixel location and outputs the resulting image.*/

template <class ImageType>
typename ImageType::Pointer Imul(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    exit(-1);
    }
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Idiv divides 2 images at every pixel location and outputs the resulting image.*/

template <class ImageType>
typename ImageType::Pointer Idiv(typename ImageType::Pointer input1, typename ImageType::Pointer input2 )
{
  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input1  );
  filter->SetInput2( input2 );
  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    exit(-1);
    }
  typename ImageType::Pointer image = filter->GetOutput();

  return image;
}

/* Iavg takes an image and the number of images as inputs , divides each pixel location of the image by the number of
  images outputs the resulting image.*/

template <class ImageType>
typename ImageType::Pointer Iavg(typename ImageType::Pointer input1, int nimgs)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(input1->GetLargestPossibleRegion() );
  image->CopyInformation(input1);
  image->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> ConstIteratorType;
  ConstIteratorType in1(input1, input1->GetLargestPossibleRegion() );
  ConstIteratorType out(image, image->GetLargestPossibleRegion() );
  for( in1.GoToBegin(), out.GoToBegin(); !in1.IsAtEnd(); ++in1, ++out )
    {
    out.Set(in1.Get() / nimgs);
    }

  return image;
}

template <class ImageType>
typename ImageType::Pointer IMask(typename ImageType::Pointer input1, typename ImageType::Pointer mask )

{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(input1->GetLargestPossibleRegion() );
  image->CopyInformation(input1);
  image->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> ConstIteratorType;
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

template <class ImageType>
typename ImageType::Pointer ImageMultiplyConstant(typename ImageType::Pointer input1,
                                                  typename ImageType::PixelType constant)
{
  typedef typename itk::ImageRegionIterator<ImageType> ConstIteratorType;
  ConstIteratorType in1(input1, input1->GetLargestPossibleRegion() );
  for( in1.GoToBegin(); !in1.IsAtEnd(); ++in1 )
    {
    in1.Set( (in1.Get() * constant) );
    }

  return input1;
}

template <class ImageType>
typename ImageType::Pointer ImageDivideConstant(typename ImageType::Pointer input1,
                                                typename ImageType::PixelType constant)
{
  typedef typename itk::ImageRegionIterator<ImageType> ConstIteratorType;
  ConstIteratorType in1(input1, input1->GetRequestedRegion() );
  for( in1.GoToBegin(); !in1.IsAtEnd(); ++in1 )
    {
    in1.Set( (in1.Get() / constant) );
    }

  return input1;
}

template <class ImageType>
void ImageSqrtValue(typename ImageType::Pointer Output,
                    const typename ImageType::Pointer Input)
{
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(Input->GetLargestPossibleRegion() );
  image->CopyInformation(Input);
  image->Allocate();
  typedef typename itk::ImageRegionIterator<ImageType> ConstIteratorType;
  ConstIteratorType in(Input, Input->GetLargestPossibleRegion() );
  ConstIteratorType out(image, image->GetLargestPossibleRegion() );
  for( in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd(); ++in, ++out )
    {
    out.Set(static_cast<typename ImageType::PixelType>(in.Get() ) );
    }

  Output = image;
  return;
}

template <class ImageType>
typename ImageType::Pointer
ImageSqrtValue(typename ImageType::Pointer input)
{
  typename ImageType::Pointer rval;
  ImageSqrtValue<ImageType>(rval, input);
  return rval;
}

#endif
