/*
 *  itkICCDeformableRegistrationFilter.txx
 *  ICCDeformationTools
 *
 *  Created by Yongqiang Zhao on 4/21/09.
 *  Copyright 2009 The University of Iowa. All rights reserved.
 *
 */

#include "itkICCDeformableRegistrationFilter.h"

#ifndef __itkICCDeformableRegistrationFilter_txx
#define __itkICCDeformableRegistrationFilter_txx

#include "itkICCDeformableRegistrationFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itkLabelStatisticsImageFilter.h"

namespace itk
{
/**
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ICCDeformableRegistrationFilter()
{
  typename ICCDeformableFunctionType::Pointer drfpf;
  drfpf = ICCDeformableFunctionType::New();

  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                 drfpf.GetPointer() ) );

  typename ICCDeformableFunctionType::Pointer drfpb;
  drfpb = ICCDeformableFunctionType::New();
  this->SetBackwardDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                         drfpb.GetPointer() ) );

  this->SetNumberOfRequiredOutputs( 2 );
  this->SetNthOutput( 0, this->MakeOutput( 0 ) );
  this->SetNthOutput( 1, this->MakeOutput( 1 ) );

  m_UpdateBuffers.reserve(this->GetNumberOfOutputs() );
  m_FFTc2rs.reserve(this->GetNumberOfOutputs() );
  m_InverseUpdateBuffers.reserve(this->GetNumberOfOutputs() );
  for( unsigned int i = 0; i < this->GetNumberOfOutputs(); i++ )
    {
    FFTWComplexToRealImagePointer fftc2r = FFTWComplexToRealImageType::New();
    m_FFTc2rs.push_back(fftc2r);

    DisplacementFieldPointer buffer = DisplacementFieldType::New();
    m_UpdateBuffers.push_back(buffer);

    DisplacementFieldPointer df = DisplacementFieldType::New();
    m_InverseUpdateBuffers.push_back(df);

    DisplacementFieldFFTPointer co = DisplacementFieldFFTType::New();
    m_Coefficients.push_back(co);
    }

  m_sqr11 = ComplexImageType::New();
  m_sqr12 = ComplexImageType::New();
  m_sqr13 = ComplexImageType::New();
  m_sqr22 = ComplexImageType::New();
  m_sqr23 = ComplexImageType::New();
  m_sqr33 = ComplexImageType::New();
  m_SmoothFilter = ComplexImageType::New();

#if 0
  m_WarpedFixedMask = MaskImageType::New();
  m_WarpedMovingMask = MaskImageType::New();

  m_FixedMask = MaskType::New();
  m_MovingMask = MaskType::New();
#endif

//  m_FixedLandmark = LandmarkType::New();
//  m_MovingLandmark = LandmarkType::New();
  m_Alpha = 0.1;
  m_Gamma = 0.0;
  m_Beta = 0.1;
  m_RegularizationWeight =  1.0;
  m_InverseWeight =     1.0;
  m_SimilarityWeight =    1.0;
  m_LandmarkWeight =    0.0;
  m_HarmonicPercent =     0.1;
  m_BackgroundFilledValue =   0.0;
  m_UseConsistentLandmark =   false;
  m_UseConsistentIntensity =  false;
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GenerateInputRequestedRegion()
{
  // call the superclass's implementation
//  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the moving image
  MovingImagePointer movingPtr =
    const_cast<MovingImageType *>( this->GetMovingImage() );

  if( movingPtr )
    {
    movingPtr->SetRequestedRegionToLargestPossibleRegion();
    }

  // just propagate up the output requested region for
  // the fixed image and initial deformation field.
  DisplacementFieldPointer inputPtr0 =
    const_cast<DisplacementFieldType *>( this->GetInput(3) );
  DisplacementFieldPointer outputPtr0 = this->GetOutput(0);

  DisplacementFieldPointer inputPtr1 =
    const_cast<DisplacementFieldType *>( this->GetInput(4) );
  DisplacementFieldPointer outputPtr1 = this->GetOutput(1);

  FixedImagePointer fixedPtr =
    const_cast<FixedImageType *>( this->GetFixedImage() );

  if( inputPtr0 )
    {
    inputPtr0->SetRequestedRegion( outputPtr0->GetRequestedRegion() );
    }

  if( inputPtr1 )
    {
    inputPtr1->SetRequestedRegion( outputPtr1->GetRequestedRegion() );
    }

  if( fixedPtr )
    {
    fixedPtr->SetRequestedRegion( outputPtr0->GetRequestedRegion() );
    }
}

/*
 * Override the default implemenation for the case when the
 * initial deformation is not set.
 * If the initial deformation is not set, the output is
 * fill with zero vectors.
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::CopyInputToOutput()
{
  typename Superclass::InputImageType::ConstPointer  inputPtr0  = this->GetInput(3);
  typename Superclass::InputImageType::ConstPointer  inputPtr1  = this->GetInput(4);

  if( inputPtr0 && inputPtr1 )
    {
//    this->Superclass::CopyInputToOutput();
    for( unsigned int i = 0; i < 2; i++ )
      {
      typename OutputImageType::ConstPointer  input  = this->GetInput(i + 3);
      typename OutputImageType::Pointer      output = this->GetOutput(i);

      if( !input || !output )
        {
        itkExceptionMacro(<< "Either input and/or output is NULL.");
        }

      // Check if we are doing in-place filtering
      if( this->GetInPlace() )
        {
        typename OutputImageType::Pointer tempPtr = output.GetPointer();
        if( tempPtr && tempPtr->GetPixelContainer() == input->GetPixelContainer() )
          {
          // the input and output container are the same - no need to copy
          return;
          }
        }

      ImageRegionConstIterator<OutputImageType> in(input, input->GetRequestedRegion() );
      ImageRegionIterator<OutputImageType>      out(output, output->GetRequestedRegion() );

      while( !out.IsAtEnd() )
        {
        out.Value() =  static_cast<typename TDisplacementField::PixelType>(in.Get() ); // Supports input image adaptors
                                                                                       // only
//			std::cout<<"out.Value:"<<out.Value()<<std::endl;
        ++in;
        ++out;
        }
      }
    }
  else
    {
    typename Superclass::PixelType zeros;
    for( unsigned int j = 0; j < 3; j++ )
      {
      zeros[j] = 0.0;
      }
    for( unsigned int i = 0; i < this->GetNumberOfOutputs(); i++ )
      {
      typename OutputImageType::Pointer output = this->GetOutput(i);

      ImageRegionIterator<OutputImageType> out(output, output->GetRequestedRegion() );

      while( !out.IsAtEnd() )
        {
        out.Value() =  zeros;
        ++out;
        }
      }
    }
}

// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField>
typename ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::ICCDeformableFunctionType
* ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::GetForwardRegistrationFunctionType()
  {
  ICCDeformableFunctionType *drfp =
    dynamic_cast<ICCDeformableFunctionType *>(this->GetDifferenceFunction().GetPointer() );

  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
  }

// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField>
const typename ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::ICCDeformableFunctionType
* ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::GetForwardRegistrationFunctionType() const
  {
  const ICCDeformableFunctionType *drfp =
    dynamic_cast<const ICCDeformableFunctionType *>(this->GetDifferenceFunction().GetPointer() );

  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
  }

// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField>
typename ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::ICCDeformableFunctionType
* ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::GetBackwardRegistrationFunctionType()
  {
  ICCDeformableFunctionType *drfp =
    dynamic_cast<ICCDeformableFunctionType *>(this->GetBackwardDifferenceFunction().GetPointer() );

  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
  }

// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField>
const typename ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::ICCDeformableFunctionType
* ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TField>
::GetBackwardRegistrationFunctionType() const
  {
  const ICCDeformableFunctionType *drfp =
    dynamic_cast<const ICCDeformableFunctionType *>(this->GetBackwardDifferenceFunction().GetPointer() );

  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
  }

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::Initialize()
{
  ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();
  ICCDeformableFunctionType *b = this->GetBackwardRegistrationFunctionType();

  if( this->GetUseConsistentIntensity() )
    {
    f->SetUseConsistentIntensity(this->GetUseConsistentIntensity() );
    b->SetUseConsistentIntensity(this->GetUseConsistentIntensity() );
    // Set different similarity weight at different resolusion
    if( this->GetOutput(0)->GetLargestPossibleRegion().GetSize()[0] > 32 )
      {
      f->SetSimilarityWeight(this->GetSimilarityWeight() / 2.0);
      b->SetSimilarityWeight(this->GetSimilarityWeight() / 2.0);
      }
    else if( this->GetOutput(0)->GetLargestPossibleRegion().GetSize()[0] > 64 )
      {
      f->SetSimilarityWeight(this->GetSimilarityWeight() / 4.0);
      b->SetSimilarityWeight(this->GetSimilarityWeight() / 4.0);
      }
    else
      {
      f->SetSimilarityWeight(this->GetSimilarityWeight() );
      b->SetSimilarityWeight(this->GetSimilarityWeight() );
      }
    }

  if( this->GetUseConsistentLandmark() )
    {
    f->SetUseConsistentLandmark(this->GetUseConsistentLandmark() );
    b->SetUseConsistentLandmark(this->GetUseConsistentLandmark() );
    // Set different landmark matching weight at different resolusion
    f->SetLandmarkWeight(this->GetLandmarkWeight() );
    b->SetLandmarkWeight(this->GetLandmarkWeight() );
#if 0
    if( this->GetOutput(0)->GetLargestPossibleRegion().GetSize()[0] > 32 )
      {
      f->SetSigma(this->GetSimilarityWeight() / 2.0);
      b->SetSigma(this->GetSimilarityWeight() / 2.0);
      }
    else if( this->GetOutput(0)->GetLargestPossibleRegion().GetSize()[0] > 64 )
      {
      f->SetSigma(this->GetSimilarityWeight() / 4.0);
      b->SetSigma(this->GetSimilarityWeight() / 4.0);
      }
    else
      {
      f->SetSigma(this->GetSimilarityWeight() );
      b->SetSigma(this->GetSimilarityWeight() );
      }
#endif
    }

  // Compute D[k], Initialize harmonic
  m_sqr11->CopyInformation(this->GetOutput(0) );
  m_sqr11->SetRegions(this->GetOutput(0)->GetLargestPossibleRegion() );
  m_sqr11->Allocate();

  m_sqr12->CopyInformation(this->GetOutput(0) );
  m_sqr12->SetRegions(this->GetOutput(0)->GetLargestPossibleRegion() );
  m_sqr12->Allocate();

  m_sqr13->CopyInformation(this->GetOutput(0) );
  m_sqr13->SetRegions(this->GetOutput(0)->GetLargestPossibleRegion() );
  m_sqr13->Allocate();

  m_sqr22->CopyInformation(this->GetOutput(0) );
  m_sqr22->SetRegions(this->GetOutput(0)->GetLargestPossibleRegion() );
  m_sqr22->Allocate();

  m_sqr23->CopyInformation(this->GetOutput(0) );
  m_sqr23->SetRegions(this->GetOutput(0)->GetLargestPossibleRegion() );
  m_sqr23->Allocate();

  m_sqr33->CopyInformation(this->GetOutput(0) );
  m_sqr33->SetRegions(this->GetOutput(0)->GetLargestPossibleRegion() );
  m_sqr33->Allocate();

  m_SmoothFilter->CopyInformation(this->GetOutput(0) );
  m_SmoothFilter->SetRegions(this->GetOutput(0)->GetLargestPossibleRegion() );
  m_SmoothFilter->Allocate();

  this->m_sqr11->FillBuffer(0.0);
  this->m_sqr12->FillBuffer(0.0);
  this->m_sqr13->FillBuffer(0.0);
  this->m_sqr22->FillBuffer(0.0);
  this->m_sqr23->FillBuffer(0.0);
  this->m_sqr33->FillBuffer(0.0);

  const float fnx = static_cast<float>(this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[0]);
  const float fny = static_cast<float>(this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[1]);
  const float fnz = static_cast<float>(this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[2]);

  const float spx = static_cast<float>(this->GetFixedImage()->GetSpacing()[0]);
  const float spy = static_cast<float>(this->GetFixedImage()->GetSpacing()[1]);
  const float spz = static_cast<float>(this->GetFixedImage()->GetSpacing()[2]);

  const float delta1 = 1.0F / fnx;
  const float delta2 = 1.0F / fny;
  const float delta3 = 1.0F / fnz;

  const float fnx2 = fnx * fnx;
  const float fny2 = fny * fny;
  const float fnz2 = fnz * fnz;

  const float sumsqdims2alpha = 2.0F * (fnx2 + fny2 + fnz2) * m_Alpha;

  ImageRegionConstIterator<ComplexImageType>
  MyIterator(this->m_sqr11, this->m_sqr11->GetRequestedRegion() );
  for( MyIterator.GoToBegin(); !MyIterator.IsAtEnd(); ++MyIterator )
    {
    const ComplexImageType::IndexType HI = MyIterator.GetIndex();
    const float                       omega1 =
      2.0F * static_cast<float>(vnl_math::pi * static_cast<float>(HI[0]) ) * delta1;
    const float omega2 =
      2.0F * static_cast<float>(vnl_math::pi * static_cast<float>(HI[1]) ) * delta2;
    const float omega3 =
      2.0F * static_cast<float>(vnl_math::pi * static_cast<float>(HI[2]) ) * delta3;

    const float alphaCosOmega2 =
      2.0F * m_Alpha * (fnx2 * vcl_cos(omega1) + fny2 * vcl_cos(omega2)
                        + fnz2 * vcl_cos(omega3) );
    const float b11 =
      ( sumsqdims2alpha + 2.0F * fnx2 * m_Beta + m_Gamma - alphaCosOmega2
        - 2.0F * fnx2 * m_Beta * vcl_cos(omega1) ) / (spx * spx);
    const float b22 =
      (sumsqdims2alpha + 2.0F * fny2 * m_Beta + m_Gamma - alphaCosOmega2
       - 2.0F * fny2 * m_Beta * vcl_cos(omega2) ) / (spy * spy);
    const float b33 =
      (sumsqdims2alpha + 2.0F * fnz2 * m_Beta + m_Gamma - alphaCosOmega2
       - 2.0F * fnz2 * m_Beta * vcl_cos(omega3) ) / (spz * spz);
    const float b12 = (fnx * fny * m_Beta * vcl_sin(omega1) * vcl_sin(omega2) ) / (spx * spy);
    const float b13 = (fnx * fnz * m_Beta * vcl_sin(omega1) * vcl_sin(omega3) ) / (spx * spz);
    const float b23 = (fny * fnz * m_Beta * vcl_sin(omega2) * vcl_cos(omega3) ) / (spy * spz);

    // Square the matrix A=BB (i.e., B'B=BB because B'=B)
    this->m_sqr11->SetPixel(HI,
                            std::complex<float>(b11 * b11 + b12 * b12 + b13 * b13, 0) );
    this->m_sqr12->SetPixel(HI,
                            std::complex<float>(b11 * b12 + b12 * b22 + b13 * b23, 0) );
    this->m_sqr13->SetPixel(HI,
                            std::complex<float>(b11 * b13 + b12 * b23 + b13 * b33, 0) );
    this->m_sqr22->SetPixel(HI,
                            std::complex<float>(b12 * b12 + b22 * b22 + b23 * b23, 0) );
    this->m_sqr23->SetPixel(HI,
                            std::complex<float>(b12 * b13 + b22 * b23 + b23 * b33, 0) );
    this->m_sqr33->SetPixel(HI,
                            std::complex<float>(b13 * b13 + b23 * b23 + b33 * b33, 0) );
    }

  // Compute smooth filter
//  float cutoff=0.0;
  const ComplexImageType::SizeType ImageDims = m_SmoothFilter->GetLargestPossibleRegion().GetSize();
  itk::ImageRegionIteratorWithIndex<ComplexImageType>
  myIterator(m_SmoothFilter, m_SmoothFilter->GetRequestedRegion() );
  for( myIterator.GoToBegin(); !myIterator.IsAtEnd(); ++myIterator )
    {
    const ComplexImageType::IndexType HI = myIterator.GetIndex();
    const unsigned int                u = HI[0];
    const unsigned int                v =
      (static_cast<unsigned int>(HI[1]) < ImageDims[1]
       / 2) ? static_cast<unsigned int>(HI[1]) : ImageDims[1] - 1 - static_cast<unsigned int>(HI[1]);
    const unsigned int w =
      (static_cast<unsigned int>(HI[2]) < ImageDims[2]
       / 2) ? static_cast<unsigned int>(HI[2]) : ImageDims[2] - 1 - static_cast<unsigned int>(HI[2]);

    float radius = 0.0;
    if( m_SmoothFilter->GetLargestPossibleRegion().GetSize()[0] > 64 )
      {
      radius = 0.01 * this->GetHarmonicPercent() * m_SmoothFilter->GetLargestPossibleRegion().GetSize()[0];
      }
    else if( m_SmoothFilter->GetLargestPossibleRegion().GetSize()[0] > 128 )
      {
      radius = 0.005 * this->GetHarmonicPercent() * m_SmoothFilter->GetLargestPossibleRegion().GetSize()[0];
      }
    else
      {
      radius = 0.01 * 1e9 * m_SmoothFilter->GetLargestPossibleRegion().GetSize()[0];
      }

    const float inv_radius = 1.0 / radius;
    const float distFromOrigin = vcl_sqrt(static_cast<float>(u * u + v * v + w * w) );
    float       ftemp = 1.0 / ( 1.0 + vcl_pow( (distFromOrigin * inv_radius), 4 ) );;
    ftemp = (ftemp < 0.01 ) ? 0.0 : ftemp;
    const std::complex<float> cmpxtemp(ftemp, 0.0F);
    myIterator.Set(cmpxtemp);
    }

  typename FFTWRealToComplexImageType::Pointer invFFT12 = FFTWRealToComplexImageType::New();
  invFFT12->SetInput(this->GetOutput(0) );
  invFFT12->Update();
  m_Coefficients[0] = invFFT12->GetOutput();
  m_Coefficients[0]->DisconnectPipeline();

  typename FFTWRealToComplexImageType::Pointer invFFT21 = FFTWRealToComplexImageType::New();
  invFFT21->SetInput(this->GetOutput(1) );
  invFFT21->Update();
  m_Coefficients[1] = invFFT21->GetOutput();
  m_Coefficients[1]->DisconnectPipeline();

  f->SetSmoothFilter(m_SmoothFilter);
  b->SetSmoothFilter(m_SmoothFilter);

  this->Superclass::Initialize();
}

/**
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::InitializeIteration()
{
  // update variables in the equation object
  ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();

  f->SetFixedImage(this->GetFixedImage() );
  f->SetMovingImage(this->GetMovingImage() );
  f->SetDisplacementField(this->GetOutput(0) );
  f->SetCoefficient(this->m_Coefficients[0]);
  f->InitializeIteration();

  ICCDeformableFunctionType *b = this->GetBackwardRegistrationFunctionType();
  b->SetFixedImage( this->GetMovingImage() );
  b->SetMovingImage( this->GetFixedImage() );
  b->SetDisplacementField(this->GetOutput(1) );
  b->SetCoefficient(this->m_Coefficients[1]);
  b->InitializeIteration();

//  Superclass::InitializeIteration();
}

/**
 *
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                TDisplacementField>
::SetMovingImageMask(MaskType *mask)
{
//  m_MovingMask = mask;
  ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();

  f->SetMovingImageMask(mask);
  f->SetBackgroundFilledValue(this->GetBackgroundFilledValue() );
  ICCDeformableFunctionType *b = this->GetBackwardRegistrationFunctionType();
  b->SetFixedImageMask(mask);
  b->SetBackgroundFilledValue(this->GetBackgroundFilledValue() );
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
const typename ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                               TDisplacementField>::MaskType
* ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                  TDisplacementField>
::GetMovingImageMask() const
  {
  const ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();
  return f->GetMovingImageMask();
  }

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                TDisplacementField>
::SetMovingLandmark(PointSetType *lk)
{
  ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();

  f->SetMovingLandmark(lk);
  ICCDeformableFunctionType *b = this->GetBackwardRegistrationFunctionType();
  b->SetFixedLandmark(lk);
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                TDisplacementField>
::SetFixedLandmark(PointSetType *lk)
{
  ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();

  f->SetFixedLandmark(lk);
  ICCDeformableFunctionType *b = this->GetBackwardRegistrationFunctionType();
  b->SetMovingLandmark(lk);
}

/**
 *
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                TDisplacementField>
::SetFixedImageMask(MaskType *mask)
{
//  m_FixedMask = mask;
  ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();

  f->SetFixedImageMask(mask);
  f->SetBackgroundFilledValue(this->GetBackgroundFilledValue() );
  ICCDeformableFunctionType *b = this->GetBackwardRegistrationFunctionType();
  b->SetMovingImageMask(mask);
  f->SetBackgroundFilledValue(this->GetBackgroundFilledValue() );
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
const typename ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                               TDisplacementField>::MaskType
* ICCDeformableRegistrationFilter<TFixedImage, TMovingImage,
                                  TDisplacementField>
::GetFixedImageMask() const
  {
  const ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();
  return f->GetFixedImageMask();
  }

/*
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
double
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GetForwardMetric() const
{
  const ICCDeformableFunctionType *drfp = this->GetForwardRegistrationFunctionType();

  return drfp->GetMetric();
}

/*
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
double
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GetBackwardMetric() const
{
  const ICCDeformableFunctionType *drfp = this->GetBackwardRegistrationFunctionType();

  return drfp->GetMetric();
}

/**
 *  Get Intensity Difference Threshold
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
double
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GetIntensityDifferenceThreshold() const
{
  const ICCDeformableFunctionType *drfp = this->GetForwardRegistrationFunctionType();

  return drfp->GetIntensityDifferenceThreshold();
}

/**
 *  Set Intensity Difference Threshold
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::SetIntensityDifferenceThreshold(double threshold)
{
  ICCDeformableFunctionType *drfpf = this->GetForwardRegistrationFunctionType();

  drfpf->SetIntensityDifferenceThreshold(threshold);
  ICCDeformableFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  drfpb->SetIntensityDifferenceThreshold(threshold);
}

/**
 *  Get Maximum Update Step Length
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
double
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GetMaximumUpdateStepLength() const
{
  const ICCDeformableFunctionType *drfp = this->GetForwardRegistrationFunctionType();

  return drfp->GetMaximumUpdateStepLength();
}

/**
 *  Set Maximum Update Step Length
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::SetMaximumUpdateStepLength(double threshold)
{
  ICCDeformableFunctionType *drfpf = this->GetForwardRegistrationFunctionType();

  drfpf->SetMaximumUpdateStepLength(threshold);

  ICCDeformableFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  drfpb->SetMaximumUpdateStepLength(threshold);
  m_MaximumUpdateStepLength = threshold;
}

/**
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
const double &
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GetForwardRMSChange() const
{
  const ICCDeformableFunctionType *drfp = this->GetForwardRegistrationFunctionType();

  return drfp->GetRMSChange();
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
const double &
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GetBackwardRMSChange() const
{
  const ICCDeformableFunctionType *drfp = this->GetBackwardRegistrationFunctionType();

  return drfp->GetRMSChange();
}

#if 0
template <class TFixedImage, class TMovingImage, class TDisplacementField>
typename ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GradientType
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::GetUseGradientType() const
{
  const ICCDeformableFunctionType *drfp = this->GetForwardRegistrationFunctionType();

  return drfp->GetUseGradientType();
}

/**
 *
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::SetUseGradientType(GradientType gtype)
{
  ICCDeformableFunctionType *drfpf = this->GetForwardRegistrationFunctionType();

  drfpf->SetUseGradientType(gtype);
  ICCDeformableFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  drfpb->SetUseGradientType(gtype);
}

#endif

/**
 *
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::AllocateUpdateBuffer()
{
  // The update buffer looks just like the output.
  for( unsigned int i = 0; i < this->GetNumberOfOutputs(); i++ )
    {
    DisplacementFieldPointer output = this->GetOutput(i);

    m_UpdateBuffers[i]->SetLargestPossibleRegion(output->GetLargestPossibleRegion() );
    m_UpdateBuffers[i]->SetRequestedRegion(output->GetRequestedRegion() );
    m_UpdateBuffers[i]->SetBufferedRegion(output->GetBufferedRegion() );
    m_UpdateBuffers[i]->SetOrigin(output->GetOrigin() );
    m_UpdateBuffers[i]->SetSpacing(output->GetSpacing() );
    m_UpdateBuffers[i]->SetDirection(output->GetDirection() );
    m_UpdateBuffers[i]->Allocate();
    }
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
typename
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>::TimeStepType
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::CalculateChange()
{
  // int threadCount;
  ICCDeformableFunctionType *f = this->GetForwardRegistrationFunctionType();
  void *                     globalData;

  globalData = f->GetGlobalDataPointer();
  m_UpdateBuffers[0] = f->GetUpdateBuffer();
  m_InverseUpdateBuffers[0] = f->GetInverseUpdateBuffer();
  m_Coefficients[0] = f->GetCoefficient();

  ICCDeformableFunctionType *b = this->GetBackwardRegistrationFunctionType();
  void *                     globalData1;
  globalData1 = b->GetGlobalDataPointer();
  m_UpdateBuffers[1] = b->GetUpdateBuffer();
  m_InverseUpdateBuffers[1] = b->GetInverseUpdateBuffer();
  m_Coefficients[1] = b->GetCoefficient();

  f->ComputeMetric(globalData);
  f->ReleaseGlobalDataPointer(globalData);

  b->ComputeMetric(globalData1);
  TimeStepType dt = b->ComputeGlobalTimeStep(globalData1);
  b->ReleaseGlobalDataPointer(globalData1);

  return dt;
}

/**
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ApplyUpdate(const TimeStepType & /* dt */)
{
  // If we smooth the update buffer before applying it, then the are
  // approximating a viscuous problem as opposed to an elastic problem
#if 0
  if( this->GetSmoothUpdateField() )
    {
    this->SmoothUpdateField();
    }
#endif

  // Compute inverse consistency
  typedef SubtractImageFilter<DisplacementFieldType,
                              DisplacementFieldType, DisplacementFieldType> SubtractDisplacementFieldType;
  typename SubtractDisplacementFieldType::Pointer sub12 = SubtractDisplacementFieldType::New();
  sub12->SetInput1(m_UpdateBuffers[0]);
  sub12->SetInput2(m_InverseUpdateBuffers[1]);
  sub12->Update();

  typename SubtractDisplacementFieldType::Pointer sub21 = SubtractDisplacementFieldType::New();
  sub21->SetInput1(m_UpdateBuffers[1]);
  sub21->SetInput2(m_InverseUpdateBuffers[0]);
  sub21->Update();

  float       normalizer_InverseConsistency;
  float       normalizer_Regularization;
  const float fnx = this->GetUpdateBuffers()[0]->GetRequestedRegion().GetSize()[0];
  const float fny = this->GetUpdateBuffers()[0]->GetRequestedRegion().GetSize()[1];
  const float fnz = this->GetUpdateBuffers()[0]->GetRequestedRegion().GetSize()[2];

  if( this->GetUpdateBuffers()[0]->GetRequestedRegion().GetSize()[0] > 64 )
    {
    m_MinJac = 0.1;
    }

  if( this->GetUpdateBuffers()[0]->GetRequestedRegion().GetSize()[0] > 128 )
    {
    m_MinJac = 0.05;
    }

  normalizer_Regularization = 4.0F * m_RegularizationWeight * m_MaximumUpdateStepLength  / (fnx * fny * fnz);
  normalizer_InverseConsistency = 4.0 * m_InverseWeight * m_MaximumUpdateStepLength;

  typename FFTWRealToComplexImageType::Pointer invfftIC12 = FFTWRealToComplexImageType::New();
  invfftIC12->SetInput(sub12->GetOutput() );
  invfftIC12->Update();

  typename FFTWRealToComplexImageType::Pointer invfftIC21 = FFTWRealToComplexImageType::New();
  invfftIC21->SetInput(sub21->GetOutput() );
  invfftIC21->Update();

  DisplacementFieldFFTPointer invfft0 = invfftIC12->GetOutput();
  invfft0->DisconnectPipeline();

  DisplacementFieldFFTPointer invfft1 = invfftIC21->GetOutput();
  invfft1->DisconnectPipeline();

  if( m_InverseWeight > 0.0 )
    {
    ComputeInverseConsistency(invfft0, invfft1, normalizer_InverseConsistency);
    }

  // Todo: Landmark matching
  // Conpute distance with itkKernelTransform and landmark spatial objects

  if( m_RegularizationWeight > 0.0 )
    {
    float for_MinJac, back_MinJac;
    for( unsigned int i = 0; i < this->GetNumberOfOutputs(); i++ )
      {
      typename FFTWComplexToRealImageType::Pointer invFFTtemp = FFTWComplexToRealImageType::New();
      invFFTtemp->SetInput(m_Coefficients[i]);
      invFFTtemp->Update();
      m_UpdateBuffers[i] = invFFTtemp->GetOutput();
      m_UpdateBuffers[i]->DisconnectPipeline();
      }

    for_MinJac = ComputeMinJac(m_UpdateBuffers[0]);
    back_MinJac = ComputeMinJac(m_UpdateBuffers[1]);

    if( for_MinJac < this->GetMinJac() )
      {
      do
        {
        this->ComputeLinearElastic(m_Coefficients[0], normalizer_Regularization);
        typename FFTWComplexToRealImageType::Pointer temp = FFTWComplexToRealImageType::New();
        temp->SetInput(m_Coefficients[0]);
        temp->Update();
        m_UpdateBuffers[0] = temp->GetOutput();
        m_UpdateBuffers[0]->DisconnectPipeline();
        for_MinJac = ComputeMinJac(m_UpdateBuffers[0]);
        std::cout << "for_MinJac:" << for_MinJac << std::endl;
        }
      while( for_MinJac < this->GetMinJac() + 0.015 );
      }
    else
      {
      this->ComputeLinearElastic(m_Coefficients[0], normalizer_Regularization);
      }

    if( back_MinJac < this->GetMinJac() )
      {
      do
        {
        this->ComputeLinearElastic(m_Coefficients[1], normalizer_Regularization);
        typename FFTWComplexToRealImageType::Pointer temp = FFTWComplexToRealImageType::New();
        temp->SetInput(m_Coefficients[1]);
        temp->Update();
        m_UpdateBuffers[1] = temp->GetOutput();
        m_UpdateBuffers[1]->DisconnectPipeline();
        back_MinJac = ComputeMinJac(m_UpdateBuffers[1]);
        std::cout << "back_MinJac:" << back_MinJac << std::endl;
        }
      while( back_MinJac < this->GetMinJac() + 0.015 );
      }
    else
      {
      this->ComputeLinearElastic(m_Coefficients[1], normalizer_Regularization);
      }
    }
//	this->ComputeLinearElastic(m_Coefficients[0]);
//	this->ComputeLinearElastic(m_Coefficients[1]);
  for( unsigned int i = 0; i < this->GetNumberOfOutputs(); i++ )
    {
    m_FFTc2rs[i]->SetInput(m_Coefficients[i]);
    m_FFTc2rs[i]->Update();
    this->GraftNthOutput(i, m_FFTc2rs[i]->GetOutput() );
    this->GetOutput(i)->Modified();
    }

  ICCDeformableFunctionType *drfp = this->GetForwardRegistrationFunctionType();

  this->SetRMSChange( drfp->GetRMSChange() );

  /**
   * Smooth the deformation field
   */
#if 0
  if( this->GetSmoothDisplacementField() )
    {
    this->SmoothDisplacementField();
    }
#endif
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Intensity difference threshold: "
     << this->GetIntensityDifferenceThreshold() << std::endl;
  os << indent << "Use First Order exponential: "
     << this->m_UseFirstOrderExp << std::endl;
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
ProcessObject::DataObjectPointer
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx)
{
  switch( idx )
    {
    case 0:
      {
      return static_cast<DataObject *>(TDisplacementField::New().GetPointer() );
      }
      break;
    case 1:
      {
      return static_cast<DataObject *>(TDisplacementField::New().GetPointer() );
      }
      break;
    default:
      return NULL;
    }
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ComputeLinearElastic(DisplacementFieldFFTPointer& coeff, float normalizer)
{
  ThreadStruct str;

  str.Filter = this;
  str.coeff = coeff;
  str.normalizer_Regularization = normalizer;
//	 str.fft = m_InvFFT12;
//	 std::cout<<"Number:"<<this->GetNumberOfThreads()<<std::endl;
  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->ComputeLinearElasticThreaderCallback, &str);
  this->GetMultiThreader()->SingleMethodExecute();
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
ITK_THREAD_RETURN_TYPE
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ComputeLinearElasticThreaderCallback(void * arg)
{
  ThreadStruct * str;
  int            total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)(arg) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)(arg) )->NumberOfThreads;

  str = (ThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)(arg) )->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if( threadId < total )
    {
    str->Filter->ThreadedComputeLinearElastic(str->coeff, str->normalizer_Regularization, splitRegion, threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ThreadedComputeLinearElastic(DisplacementFieldFFTPointer& coeff, float rho4delta_normalizer,
                               const ThreadRegionType & regionToProcess,
                               int)
{
  ImageRegionIterator<DisplacementFieldFFTType> CoeffsIterator(coeff, regionToProcess);
  for( CoeffsIterator.GoToBegin();
       !CoeffsIterator.IsAtEnd(); // && !Coeffs2Iterator.IsAtEnd() && !Coeffs3Iterator.IsAtEnd();
       ++CoeffsIterator )
    {
    const DisplacementFieldFFTType::IndexType HI = CoeffsIterator.GetIndex();

    // a Gauss-Seidel modification of gradient decent.
    const float dsqr11 = m_sqr11->GetPixel(HI).real();
    const float dsqr12 = m_sqr12->GetPixel(HI).real();
    const float dsqr13 = m_sqr13->GetPixel(HI).real();
    const float dsqr22 = m_sqr22->GetPixel(HI).real();
    const float dsqr23 = m_sqr23->GetPixel(HI).real();
    const float dsqr33 = m_sqr33->GetPixel(HI).real();

    DisplacementFieldFFTType::PixelType pixel = CoeffsIterator.Get();

    const float c1real = pixel[0].real(); // Coeffs1Iterator.Get().real();
    const float c1imag = pixel[0].imag(); // Coeffs1Iterator.Get().imag();
    const float c2real = pixel[1].real(); // Coeffs2Iterator.Get().real();
    const float c2imag = pixel[1].imag(); // Coeffs2Iterator.Get().imag();
    const float c3real = pixel[2].real(); // Coeffs3Iterator.Get().real();
    const float c3imag = pixel[2].imag(); // Coeffs3Iterator.Get().imag();
    pixel[0]
      -= std::complex<float>(
          +rho4delta_normalizer
          * (
            +dsqr11 * c1real
            + dsqr12 * c2real
            + dsqr13 * c3real
            ),
          +rho4delta_normalizer
          * (
            +dsqr11 * c1imag
            + dsqr12 * c2imag
            + dsqr13 * c3imag
            )
          );
    /*coeffs2( i,j,k )*/
    pixel[1]
      -= std::complex<float>(
          +rho4delta_normalizer
          * (
            +dsqr12 * c1real
            + dsqr22 * c2real
            + dsqr23 * c3real
            ),
          +rho4delta_normalizer
          * (
            +dsqr12 * c1imag
            + dsqr22 * c2imag
            + dsqr23 * c3imag
            )
          );
    /*coeffs3( i,j,k )*/
    pixel[2]
      -= std::complex<float>(
          +rho4delta_normalizer
          * (
            +dsqr13 * c1real
            + dsqr23 * c2real
            + dsqr33 * c3real
            ),
          +rho4delta_normalizer
          * (
            +dsqr13 * c1imag
            + dsqr23 * c2imag
            + dsqr33 * c3imag
            )
          );
    CoeffsIterator.Set(pixel);
    }
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ComputeInverseConsistency(DisplacementFieldFFTPointer& inv0, DisplacementFieldFFTPointer& inv1, float normalizer)
{
  ThreadStruct str;

  str.Filter = this;
  str.inverse0 = inv0;
  str.inverse1 = inv1;
  str.normalizer_InverseConsistency = normalizer;
//	 str.coeff = coeff;
//	 std::cout<<"Number:"<<this->GetNumberOfThreads()<<std::endl;
  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->ComputeInverseConsistencyThreaderCallback, &str);
  this->GetMultiThreader()->SingleMethodExecute();
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
ITK_THREAD_RETURN_TYPE
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ComputeInverseConsistencyThreaderCallback(void * arg)
{
  ThreadStruct * str;
  int            total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)(arg) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)(arg) )->NumberOfThreads;

  str = (ThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)(arg) )->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if( threadId < total )
    {
    str->Filter->ThreadedComputeInverseConsistency(str->inverse0, str->inverse1, str->normalizer_InverseConsistency,
                                                   splitRegion,
                                                   threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
void
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ThreadedComputeInverseConsistency(DisplacementFieldFFTPointer& inv0, DisplacementFieldFFTPointer& inv1,
                                    float normalizer, const ThreadRegionType & regionToProcess,
                                    int)
{
  ImageRegionIterator<DisplacementFieldFFTType>      coeffsIter0(m_Coefficients[0], regionToProcess);
  ImageRegionConstIterator<DisplacementFieldFFTType> iter0(inv0, regionToProcess);
  ImageRegionIterator<DisplacementFieldFFTType>      coeffsIter1(m_Coefficients[1], regionToProcess);
  ImageRegionConstIterator<DisplacementFieldFFTType> iter1(inv1, regionToProcess);
  for( iter0.GoToBegin(), coeffsIter0.GoToBegin(),
       iter1.GoToBegin(), coeffsIter1.GoToBegin();
       !iter1.IsAtEnd();
       ++iter0,  ++coeffsIter0,
       ++iter1,  ++coeffsIter1 )
    {
    typename DisplacementFieldFFTType::PixelType pixel0 = coeffsIter0.Get();
    typename DisplacementFieldFFTType::PixelType pixel1 = coeffsIter1.Get();
/*
    typename DisplacementFieldType::IndexType index = dfIter.GetIndex();
        typename MovingImageType::PointType fixedPoint;
        this->GetMovingImage()->TransformIndexToPhysicalPoint(index, fixedPoint);
        typename DisplacementFieldType::PixelType pixel = dfIter.Get();
        if( !this->GetMovingImageMask()->IsInside(fixedPoint) &&  (m_WarpedFixedMask->GetPixel(index) < 1) )
        {
          pixel.Fill(0.0);
          dfIter1.Set(pixel);
        }
*/
    for( unsigned i = 0; i < 3; i++ )
      {
      pixel0[i] = pixel0[i] - iter0.Get()[i] * normalizer;
      pixel1[i] = pixel1[i] - iter1.Get()[i] * normalizer;
      }
    coeffsIter0.Set(pixel0);
    coeffsIter1.Set(pixel1);
    }
  m_Coefficients[0]->Modified();
  m_Coefficients[1]->Modified();
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
unsigned int
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::SplitRequestedRegion(unsigned int i, unsigned int num, OutputImageRegionType& splitRegion)
{
  const typename TFixedImage::SizeType& requestedRegionSize
    = m_Coefficients[0]->GetRequestedRegion().GetSize();

  int splitAxis;
  typename TFixedImage::IndexType splitIndex;
  typename TFixedImage::SizeType splitSize;

  // Initialize the splitRegion to the output requested region
  splitRegion = m_Coefficients[0]->GetRequestedRegion();
  splitIndex = splitRegion.GetIndex();
  splitSize = splitRegion.GetSize();

  // split on the outermost dimension available
  splitAxis = m_Coefficients[0]->GetImageDimension() - 1;
  while( requestedRegionSize[splitAxis] == 1 )
    {
    --splitAxis;
    if( splitAxis < 0 )
      { // cannot split
      itkDebugMacro("  Cannot Split");
      return 1;
      }
    }

  // determine the actual number of pieces that will be generated
  typename TFixedImage::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
  const unsigned int valuesPerThread = (unsigned int)::vcl_ceil(range / (double)num);
  const unsigned int maxThreadIdUsed = (unsigned int)::vcl_ceil(range / (double)valuesPerThread) - 1;

  // Split the region
  if( i < maxThreadIdUsed )
    {
    splitIndex[splitAxis] += i * valuesPerThread;
    splitSize[splitAxis] = valuesPerThread;
    }
  if( i == maxThreadIdUsed )
    {
    splitIndex[splitAxis] += i * valuesPerThread;
    // last thread needs to process the "rest" dimension being split
    splitSize[splitAxis] = splitSize[splitAxis] - i * valuesPerThread;
    }

  // set the split region ivars
  splitRegion.SetIndex( splitIndex );
  splitRegion.SetSize( splitSize );

  itkDebugMacro("  Split Piece: " << splitRegion );

  return maxThreadIdUsed + 1;
}

template <class TFixedImage, class TMovingImage, class TDisplacementField>
double
ICCDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField>
::ComputeMinJac(DisplacementFieldPointer& deffield)
{
  typedef itk::DisplacementFieldJacobianDeterminantFilter<
      TDisplacementField, typename TFixedImage::PixelType, TFixedImage> JacobianFilterType;
  typename JacobianFilterType::Pointer m_JacobianFilter = JacobianFilterType::New();
  m_JacobianFilter->SetUseImageSpacing( true );
  m_JacobianFilter->ReleaseDataFlagOn();
  m_JacobianFilter->SetInput( deffield );
  m_JacobianFilter->UpdateLargestPossibleRegion();

  const unsigned int numPix = m_JacobianFilter->
    GetOutput()->GetLargestPossibleRegion().
    GetNumberOfPixels();

  float *pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
  float *pix_end = pix_start + numPix;

  // Get min jac
  double minJac = *( std::min_element(pix_start, pix_end) );
  return minJac;
}
} // end namespace itk

#endif
