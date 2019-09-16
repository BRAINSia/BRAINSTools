//
// Created by Hans Johnson on 7/13/16.
//
#ifndef SRTypes_h_
#define SRTypes_h_

#include "itkComplexToRealImageFilter.h"

#include "itkFFTWForwardFFTImageFilter.h"
#include "itkFFTWInverseFFTImageFilter.h"
#include "itkFFTWRealToHalfHermitianForwardFFTImageFilter.h"
#include "itkFFTWHalfHermitianToRealInverseFFTImageFilter.h"
#include "itkComplexToComplexFFTImageFilter.h"
#include "itkFFTWComplexToComplexFFTImageFilter.h"
#include <itkVectorImage.h>
//#include <itkGradientImageFilter.h>
#include <rtkForwardDifferenceGradientImageFilter.h>
#include <rtkBackwardDifferenceDivergenceImageFilter.h>


#ifdef USE_SIMPLEITK
#  include <SimpleITK.h>
#  include <sitkImageOperators.h>

namespace sitk = itk::simple;

// TODO, perhaps this should be part of SITK proper?
inline sitk::Image
ImageFill(sitk::Image & img1, double s)
{
  return sitk::Add(sitk::Multiply(img1, 0.0), s);
}
extern sitk::Image
SimpleOpWeightedL2(sitk::Image & norm01_lowres, sitk::Image & edgemask);

#endif


// using PrecisionType = double;
using PrecisionType = float;

#define CURR_IMG_DIM 3

using FloatImageType = itk::Image<PrecisionType, CURR_IMG_DIM>;
using ImageBaseType = itk::ImageBase<CURR_IMG_DIM>;

//#define USE_HALF_FFTW
#ifdef USE_HALF_FFTW
// using FloatFFTWFullFFTType = itk::FFTWRealToHalfHermitianForwardFFTImageFilter< FloatImageType >;
// using FloatFFTWFullIFFTType = itk::FFTWHalfHermitianToRealInverseFFTImageFilter<
// FloatFFTWFullFFTType::OutputImageType >; using HalfHermetianImageType = FloatFFTWFullFFTType::OutputImageType;
#else
using FloatFFTWFullFFTType = itk::FFTWForwardFFTImageFilter<FloatImageType>;
using FloatFFTWFullIFFTType = itk::FFTWInverseFFTImageFilter<FloatFFTWFullFFTType::OutputImageType>;
using HalfHermetianImageType = FloatFFTWFullFFTType::OutputImageType;

using ComplexFFTWFullIFFT = itk::FFTWComplexToComplexFFTImageFilter<FloatFFTWFullFFTType::OutputImageType>;
using C2FType = itk::ComplexToRealImageFilter<FloatFFTWFullFFTType::OutputImageType, FloatImageType>;
#endif

using CVType = itk::CovariantVector<PrecisionType, CURR_IMG_DIM>;
using CVImageType = itk::Image<CVType, CURR_IMG_DIM>;
using GradientType =
  rtk::ForwardDifferenceGradientImageFilter<FloatImageType, PrecisionType, PrecisionType, CVImageType>;
using DivergenceType = rtk::BackwardDifferenceDivergenceImageFilter<CVImageType, FloatImageType>;

using VarVecImageType = itk::VectorImage<PrecisionType, CURR_IMG_DIM>;
using VarVecType = VarVecImageType::PixelType;
using GradientVarVecType =
  rtk::ForwardDifferenceGradientImageFilter<FloatImageType, PrecisionType, PrecisionType, VarVecImageType>;
using DivergenceVarVecType = rtk::BackwardDifferenceDivergenceImageFilter<VarVecImageType, FloatImageType>;

template <typename ImageType>
typename ImageType::Pointer
CreateEmptyImage(typename itk::ImageBase<ImageType::ImageDimension> * in)
{
  typename ImageType::Pointer out = ImageType::New();
  out->CopyInformation(in);
  out->SetRegions(in->GetLargestPossibleRegion());
  out->Allocate();
  out->FillBuffer(itk::NumericTraits<typename ImageType::PixelType>::ZeroValue());
  return out;
}

template <typename ImageType>
typename ImageType::Pointer
CreateEmptyImage(typename ImageType::Pointer in)
{
  return CreateEmptyImage<ImageType>(in.GetPointer());
}


#endif // SRTypes_h_
