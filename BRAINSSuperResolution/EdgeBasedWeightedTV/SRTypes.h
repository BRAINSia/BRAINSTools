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
#include <SimpleITK.h>
#include <sitkImageOperators.h>

namespace sitk = itk::simple;

//TODO, perhaps this should be part of SITK proper?
inline sitk::Image ImageFill( sitk::Image &img1, double s) { return sitk::Add(sitk::Multiply(img1, 0.0),s); }
extern sitk::Image SimpleOpWeightedL2(sitk::Image & norm01_lowres, sitk::Image & edgemask);

#endif


//typedef double PrecisionType;
typedef float PrecisionType;

#define CURR_IMG_DIM 3

typedef itk::Image<PrecisionType, CURR_IMG_DIM> FloatImageType;
typedef itk::ImageBase<CURR_IMG_DIM> ImageBaseType;

//#define USE_HALF_FFTW
#ifdef USE_HALF_FFTW
//typedef itk::FFTWRealToHalfHermitianForwardFFTImageFilter< FloatImageType > FloatFFTWFullFFTType;
//typedef itk::FFTWHalfHermitianToRealInverseFFTImageFilter< FloatFFTWFullFFTType::OutputImageType > FloatFFTWFullIFFTType;
//typedef FloatFFTWFullFFTType::OutputImageType HalfHermetianImageType;
#else
typedef itk::FFTWForwardFFTImageFilter< FloatImageType > FloatFFTWFullFFTType;
typedef itk::FFTWInverseFFTImageFilter< FloatFFTWFullFFTType::OutputImageType > FloatFFTWFullIFFTType;
typedef FloatFFTWFullFFTType::OutputImageType HalfHermetianImageType;

typedef itk::FFTWComplexToComplexFFTImageFilter< FloatFFTWFullFFTType::OutputImageType > ComplexFFTWFullIFFT;
typedef itk::ComplexToRealImageFilter< FloatFFTWFullFFTType::OutputImageType, FloatImageType> C2FType;
#endif

typedef itk::CovariantVector < PrecisionType, CURR_IMG_DIM >  CVType;
typedef itk::Image< CVType, CURR_IMG_DIM > CVImageType;
typedef rtk::ForwardDifferenceGradientImageFilter<FloatImageType,PrecisionType ,PrecisionType, CVImageType> GradientType;
typedef rtk::BackwardDifferenceDivergenceImageFilter<CVImageType,FloatImageType > DivergenceType;

typedef itk::VectorImage< PrecisionType, CURR_IMG_DIM > VarVecImageType;
typedef VarVecImageType::PixelType VarVecType;
typedef rtk::ForwardDifferenceGradientImageFilter<FloatImageType,PrecisionType ,PrecisionType, VarVecImageType> GradientVarVecType;
typedef rtk::BackwardDifferenceDivergenceImageFilter<VarVecImageType,FloatImageType > DivergenceVarVecType;

template<typename ImageType>
typename ImageType::Pointer CreateEmptyImage(typename itk::ImageBase<ImageType::ImageDimension> * in)
{
  typename ImageType::Pointer out = ImageType::New();
  out->CopyInformation(in);
  out->SetRegions(in->GetLargestPossibleRegion());
  out->Allocate();
  out->FillBuffer(itk::NumericTraits<typename ImageType::PixelType>::ZeroValue());
  return out;
}

template<typename ImageType>
typename ImageType::Pointer CreateEmptyImage(typename ImageType::Pointer in)
{
  return CreateEmptyImage<ImageType>(in.GetPointer());
}



#endif //SRTypes_h_
