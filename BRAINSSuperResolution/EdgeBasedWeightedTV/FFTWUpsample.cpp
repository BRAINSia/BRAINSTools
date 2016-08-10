//
// Created by Hans Johnson on 7/12/16.
//

// \author Hans J. Johnson
// \date 2016-07-10
// Utility function for upsampling by FFT

#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "FFTWUpsample.h"

#include <itkPeriodicBoundaryCondition.h>

typedef itk::ImageRegionIteratorWithIndex<HalfHermetianImageType> cmplxHHIteratorType;

/**
 * @author Hans J. Johnson
 * @brief This function performs the index-wise movement of coefficients for padding or
 *        downsampling in the frequency domain.
 * @param outputIter The current location in the output coeffs to be filled in
 * @param outSize    Size of the output coefficients
 * @param outFirstSpatialDeminsionIsOdd Do these coeffs represent an image with odd first dimension?
 * @param inputFreqCoeffs The current location in the input coeffs to be filled in
 * @param inSize Size of the input coefficients
 * @param inFirstSpatialDimensionIsOdd Do these coeffs represent an image with odd first dimension?
 */
static void ReshapeConvertIndex2(cmplxHHIteratorType & outputIter,
                                 const HalfHermetianImageType::SizeType &outSize,
                                 const bool outFirstSpatialDeminsionIsOdd,
                                 const HalfHermetianImageType::ConstPointer inputFreqCoeffs,
                                 const HalfHermetianImageType::SizeType & inSize,
                                 const bool inFirstSpatialDimensionIsOdd)
{
  PrecisionType fftScale=1.0F;
  //Return if ths new index is in the outSize.  Otherwise you should truncate
  HalfHermetianImageType::IndexType inHighestFreq;
  HalfHermetianImageType::IndexType outHighestFreq;

  //The nyq_OutOffset[dim] is the offest fromt he outHighestFreq to get the index of the NyquistFrequency
  HalfHermetianImageType::IndexValueType nyq_OutOffset[HalfHermetianImageType::ImageDimension];
  HalfHermetianImageType::IndexValueType nyq_InOffset[HalfHermetianImageType::ImageDimension];

#ifdef USE_HALF_FFTW
  nyq_OutOffset[0]=(outFirstSpatialDeminsionIsOdd)? 0 : 1;
  nyq_InOffset[0]= (inFirstSpatialDimensionIsOdd) ? 0 : 1;

  outHighestFreq[0] = outSize[0] - 1 - nyq_OutOffset[0];
  inHighestFreq[0]  = inSize[0]  - 1 - nyq_InOffset[0];
  for(size_t dim = 1; dim < HalfHermetianImageType::ImageDimension; ++dim)
#else
  for(size_t dim = 0; dim < HalfHermetianImageType::ImageDimension; ++dim)
#endif
  {
    inHighestFreq[dim]  =  (inSize[dim]-1)/2; // if sz(0) = 6 -> HF(0) = 2, and NyquistFrequency=3, isInSizeOdd=false
                                              // if sz(0) = 5 -> HF(0) = 2, and NyquistFrequency=NA, isInSizeOdd=true
    outHighestFreq[dim] = (outSize[dim]-1)/2;
    nyq_OutOffset[dim] = (outSize[dim]+1)%2;//One if even
    nyq_InOffset[dim] =  (inSize[dim] +1)%2;//Zero if odd
  }
  const HalfHermetianImageType::IndexType outputIndex= outputIter.GetIndex();
#ifdef USE_HALF_FFTW
  if ( outputIndex[0] > ( inSize[0] - 1 ) )
  {
     return;
  }
#endif
  HalfHermetianImageType::IndexType inputIndex = outputIndex;

  //First find shifted index location
  for(size_t dim = 0 ; dim < HalfHermetianImageType::ImageDimension; ++dim )
  {
#ifdef USE_HALF_FFTW
    if( dim == 0 && outputIndex[dim] == outHighestFreq[dim] + nyq_OutOffset[dim])
    {
      //No manipulation of index needed
    }
    else
#endif
    if (outputIndex[dim] >= ( outHighestFreq[dim] + nyq_OutOffset[dim]) ) // Shift index to negative quadrant
    {
      inputIndex[dim]  = outputIndex[dim] - outSize[dim];
      if ( inputIndex[dim] < -(inHighestFreq[dim] ) )
      {
        return;
      }
      inputIndex[dim] += inSize[dim];
    }
    else if ( outputIndex[dim] >  ( inHighestFreq[dim] + nyq_InOffset[dim]) )
    {
      return;
    }

     //TODO: HACK FOR debugging
    if(inputIndex[dim] >= inSize[dim]  || inputIndex[dim] < 0 )
    {
      std::cout << "XXX: " << dim << " : " << inputIndex[dim] << std::endl;
      return;
    }
#ifdef USE_HALF_FFTW
    //First process the input nyquist frequency if output index matches, the nyquist frequency must by split.
    if(dim == 0)
    {
      if( (nyq_InOffset[0] == 1 ) && ( outputIndex[0] == (inSize[0] - 1) ) && outSize[0] != inSize[0])
      {
        //For even signals during an upsampling,
        fftScale = .5;
      }
      if( (nyq_OutOffset[0] == 1 ) && ( inputIndex[0] == (outSize[0] - 1) ) && outSize[0] != inSize[0])
      {
        //For even signals during an upsampling,
        fftScale = 2;
      }
    }
#endif

  }

  fftScale = 1.0F; //HACK Need to review this later, Perhaps 0.5 is the right answer to half-hermitian

 //HalfHermetianImageType::IndexType origIndex = outputIter.GetIndex();
 // if (origIndex[0] == 0 && origIndex[1] == 0 && origIndex[2] == 0  )
 //   std::cout << outputIter.GetIndex() << origIndex << "\t" << inputIndex <<  "\t\t\t" << inputFreqCoeffs->GetPixel
 //       (inputIndex) << std::endl;
  //TODO: Use more effiecient form once verified
 // outputIter.Set(outputIter.Value()+fftScale*inputFreqCoeffs->GetPixel(inputIndex));
  outputIter.Set(fftScale*inputFreqCoeffs->GetPixel(inputIndex));
}

#ifdef DEBUG
void testIndexConversion()
{
  std::cout << "======================================" << std::endl;
  {//Test upsampling
    const HalfHermetianImageType::SizeType inSize = { 3, 4, 4};
    const HalfHermetianImageType::SizeType outSize = {5,8,8};
    HalfHermetianImageType::IndexType oIndex;
    size_t count = 0;
    for(size_t k = 0 ; k < inSize[2]; ++k)
      for(size_t j = 0 ; j < inSize[1]; ++j)
        for(size_t i = 0 ; i < inSize[0]; ++i)
        {
          oIndex[0] = i;
          oIndex[1] = j;
          oIndex[2] = k;
          const HalfHermetianImageType::IndexType in=oIndex;
          const bool isInside = ConvertLRI2HRI(oIndex,inSize,outSize);
          if(isInside)
          {
            const HalfHermetianImageType::IndexType out=oIndex;
            std::cout << in << "\t\t" << out << "\t\t" << count << std::endl;
            count++;
          }
        }
  }
  std::cout << "======================================" << std::endl;
  {//Test downsampling
    const HalfHermetianImageType::SizeType inSize = {5,8,8};
    const HalfHermetianImageType::SizeType outSize = { 3, 4, 4};
    HalfHermetianImageType::IndexType oIndex;
    size_t count = 0;
    for(size_t k = 0 ; k < inSize[2]; ++k)
      for(size_t j = 0 ; j < inSize[1]; ++j)
        for(size_t i = 0 ; i < inSize[0]; ++i)
        {
          oIndex[0] = i;
          oIndex[1] = j;
          oIndex[2] = k;
          const HalfHermetianImageType::IndexType in=oIndex;
          const bool isInside = ConvertLRI2HRI(oIndex,inSize,outSize);
          if(isInside)
          {
            const HalfHermetianImageType::IndexType out=oIndex;
            std::cout << in << "\t\t" << out << "\t\t" << count << std::endl;
            count++;
          }
        }
  }
  std::cout << "======================================" << std::endl;
}
#endif



#include <itkFFTWGlobalConfiguration.h>
//This is intended to be called one time
void FFTWInit(const std::string path_for_wisdom)
{
  //Environmental variables
  //itksys::SystemTools::GetEnv("ITK_FFTW_PLAN_RIGOR", "STRING");
  //ITK_FFTW_PLAN_RIGOR   - Defines how aggressive the generation of wisdom should be.
  //         FFTW_ESTIMATE - super fast guess to plan, and often mediocre performance for odd sizes
  //         FFTW_MEASURE - quick planning wiht heuristics, and usually good performance
  //         FFTW_PATIENT -    slow planning with heuristics, and almost always best performance (skip testing
  // odd-ball cases)
  //         FFTW_EXHAUSTIVE - very slow planning by checking odd-ball cases, often does not produce better results
  //ITK_FFTW_READ_WISDOM_CACHE  - Defines if a wisdom file cache should
  //                              be read if found.  (it is "On" by default)
  //ITK_FFTW_WRITE_WISDOM_CACHE  - Defines if generated wisdom file cache
  //                               should be written (it is "Off" by default)
  //ITK_FFTW_WISDOM_CACHE_BASE - Defines the base directory where the
  //                             fftw wisdom cache will be placed,
  //                             this is intended to be used with auto-
  //                             generated cache file names
  //ITK_FFTW_WISDOM_CACHE_FILE - Defines the full name of the cache
  //                             file to be generated.  If this is
  //                             set, then ITK_FFTW_WISDOM_CACHE_BASE
  //                             is ignored.
  if( path_for_wisdom.length() > 0) // If empty, just use the default.
  {
    itk::FFTWGlobalConfiguration::SetWisdomCacheBase(path_for_wisdom);
  }
  itk::FFTWGlobalConfiguration::SetReadWisdomCache(true);
  itk::FFTWGlobalConfiguration::SetWriteWisdomCache(true);
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileDouble();
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileFloat();
  itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_EXHAUSTIVE);
  //itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_MEASURE);
  //std::string temp = itk::FFTWGlobalConfiguration::GetWisdomFileDefaultBaseName();
  //FFTW_MEASURE, FFTW_PATIENT, or FFTW_EXHAUSTIVE
  // We need to ensure that we read in the default wisdom files from the default locations
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileFloat();
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileDouble();
  std::cout << itk::FFTWGlobalConfiguration::GetWisdomCacheBase() << std::endl;
  std::cout << itk::FFTWGlobalConfiguration::GetWisdomFileDefaultBaseName() << std::endl;
}

//This is slow. It would be better to re-use fft filter many times if possible to save memory allocation
//FFTScalar by default is set to 1.0F.
//We often want to scale by 1.0F/static_cast<PrecisionType>(inputImage->GetLargestPossibleRegion().GetNumberOfPixels());
HalfHermetianImageType::Pointer GetForwardFFT(FloatImageType::Pointer inputImage, const PrecisionType FFTScaler)
{
  FloatFFTWFullFFTType::Pointer fft = FloatFFTWFullFFTType::New();
  fft->SetInput(inputImage);
  fft->Update();
  HalfHermetianImageType::Pointer outputCoeffs = fft->GetOutput();
  if (FFTScaler != 1.0F)
  {
    typedef itk::MultiplyImageFilter<HalfHermetianImageType, HalfHermetianImageType> MultType;
    MultType::Pointer multFilt = MultType::New();
    multFilt->SetInput(outputCoeffs);
    multFilt->SetInPlace(true);
    multFilt->SetConstant(FFTScaler);
    multFilt->Update();
    outputCoeffs = multFilt->GetOutput();
  }
  //Bring metatda allong for ride
  outputCoeffs->SetSpacing(   inputImage->GetSpacing() );
  outputCoeffs->SetOrigin(    inputImage->GetOrigin() );
  outputCoeffs->SetDirection( inputImage->GetDirection() );
  return outputCoeffs;
}


FloatImageType::Pointer GetInverseFFT(HalfHermetianImageType::Pointer inputFFTCoeffs,
                                      const bool referenceImageBase_ActualXDimensionIsOdd, const PrecisionType FFTScaler)
{
#ifdef USE_HALF_FFTW
  static FloatFFTWFullIFFTType::Pointer ifft = FloatFFTWFullIFFTType::New();
  ifft->SetInput(inputFFTCoeffs);
#ifdef USE_HALF_FFTW
  ifft->SetActualXDimensionIsOdd(referenceImageBase_ActualXDimensionIsOdd);
#endif
  ifft->Update();
  FloatImageType::Pointer referenceImageBaseOut = ifft->GetOutput();
#else
//#define USE_FULL_COMPLEX_INVERSE // This is mostly working, turning off for testing.
#ifdef USE_FULL_COMPLEX_INVERSE
  ComplexFFTWFullIFFT::Pointer C2CIFFT = ComplexFFTWFullIFFT::New();
  C2CIFFT->SetInput(inputFFTCoeffs);
  C2CIFFT->SetTransformDirection(ComplexFFTWFullIFFT::INVERSE);
  C2CIFFT->Update();
  C2FType::Pointer c2f = C2FType::New();
  c2f->SetInput(C2CIFFT->GetOutput());
  c2f->Update();
  FloatImageType::Pointer referenceImageBaseOut = c2f->GetOutput();
#else
  FloatFFTWFullIFFTType::Pointer ifft = FloatFFTWFullIFFTType::New();
  ifft->SetInput(inputFFTCoeffs);
#ifdef USE_HALF_FFTW
  ifft->SetActualXDimensionIsOdd(referenceImageBase_ActualXDimensionIsOdd);
#endif
  ifft->Update();
  FloatImageType::Pointer referenceImageBaseOut = ifft->GetOutput();
#endif
#endif
  //Bring metatda allong for ride
  if (FFTScaler != 1.0F)
  {
    typedef itk::MultiplyImageFilter<FloatImageType, FloatImageType> MultType;
    MultType::Pointer multFilt = MultType::New();
    multFilt->SetInput(referenceImageBaseOut);
    //HACK multFilt->SetInPlace(true);
    multFilt->SetConstant(FFTScaler);
    multFilt->Update();
    referenceImageBaseOut = multFilt->GetOutput();
  }
  referenceImageBaseOut->SetSpacing(   inputFFTCoeffs->GetSpacing() );
  referenceImageBaseOut->SetOrigin(    inputFFTCoeffs->GetOrigin() );
  referenceImageBaseOut->SetDirection( inputFFTCoeffs->GetDirection() );
  return referenceImageBaseOut;
}


//Create an appropriate memory layout for holding the output of an FFT based on
//the size of the referenceImageBase information.
HalfHermetianImageType::Pointer CreateZeroFFTCoefficients(itk::ImageBase<3>::Pointer referenceImageBase)
{
  FloatImageType::RegionType newRegion = referenceImageBase->GetLargestPossibleRegion();
#ifdef USE_HALF_FFTW
  {
    FloatImageType::SizeType temp_size = referenceImageBase->GetLargestPossibleRegion().GetSize();
    temp_size[0] = (temp_size[0]/2 ) + 1; //Only store 1/2 of the first dimension to take advantage of symmetry
    newRegion.SetSize(temp_size);
  }
#endif
  HalfHermetianImageType::Pointer outputFreqCoeffs = HalfHermetianImageType::New();
  outputFreqCoeffs->CopyInformation(referenceImageBase);
  outputFreqCoeffs->SetRegions(newRegion);
  outputFreqCoeffs->Allocate();
  outputFreqCoeffs->FillBuffer(itk::NumericTraits<HalfHermetianImageType::PixelType>::ZeroValue());
  return outputFreqCoeffs;
}

void MoveFFTCoeffs(HalfHermetianImageType::Pointer outputFreqCoeffs,
                   const bool outFirstSpatialDeminsionIsOdd,
                   HalfHermetianImageType::Pointer inputFreqCoeffs,
                   const bool inFirstSpatialDeminsionIsOdd )
{
  outputFreqCoeffs->FillBuffer(std::complex<PrecisionType >(0.0F,0.0F));
  cmplxHHIteratorType outputHHIter(outputFreqCoeffs, outputFreqCoeffs->GetLargestPossibleRegion());

  const HalfHermetianImageType::SizeType outSize = outputFreqCoeffs->GetLargestPossibleRegion().GetSize();
  const HalfHermetianImageType::SizeType inSize =  inputFreqCoeffs->GetLargestPossibleRegion().GetSize();

  for (outputHHIter.GoToBegin(); !outputHHIter.IsAtEnd(); ++outputHHIter)
  {
    ReshapeConvertIndex2(outputHHIter,outSize,outFirstSpatialDeminsionIsOdd,
                         inputFreqCoeffs.GetPointer(),inSize, inFirstSpatialDeminsionIsOdd);
  }
}

// Transfer FFT coefficients to/from different sized frequency domains (i.e. upsample/downsample image)
// downsampling has the effect of low-pass filtering.
// implement out = p(ind_samples)
// FFTScalar is a scale factor applied uniformly during the reshape process
HalfHermetianImageType::Pointer ReshapeFFT(
  const itk::ImageBase<3>::Pointer outputRealImageBase,
  HalfHermetianImageType::Pointer inputFreqCoeffs,
  const bool inFirstSpatialDeminsionIsOdd)
{
  const bool FirstDimensionIsOdd = outputRealImageBase->GetLargestPossibleRegion().GetSize()[0] %2 == 1;
  HalfHermetianImageType::Pointer outputFreqCoeffs = CreateZeroFFTCoefficients(outputRealImageBase);
  cmplxHHIteratorType outputHHIter(outputFreqCoeffs, outputFreqCoeffs->GetLargestPossibleRegion());

  MoveFFTCoeffs(outputFreqCoeffs,FirstDimensionIsOdd,inputFreqCoeffs,inFirstSpatialDeminsionIsOdd );
  return outputFreqCoeffs;
}

//Produce FFT coefficients from inputImage, but return FFT coefficients that are based on
//the referenceImageBase (up/down sampled coefficients) defines the desired output FFT Space
// low_resolution samples
// inHRRealImage is always a high resolution image
// lpfSamplesDef is always a low-resolution example
// referenceImageBase is always a high-resolution example
HalfHermetianImageType::Pointer A_fhp(FloatImageType::Pointer inHRRealImage,
                                      itk::ImageBase<3>::Pointer desiredOutputRef)
{
  const bool inputFirstDimIsOdd = inHRRealImage->GetLargestPossibleRegion().GetSize()[0] % 2 == 1;
  const PrecisionType FFTScaler =
    1.0F / std::sqrt(static_cast<PrecisionType>(inHRRealImage->GetLargestPossibleRegion().GetNumberOfPixels()));
  HalfHermetianImageType::Pointer inputImage_cmplHH = GetForwardFFT(inHRRealImage, FFTScaler);
  //==========================
  // Transfer FFT coeficients to new space
  HalfHermetianImageType::Pointer outputFreqCoeffs = ReshapeFFT(desiredOutputRef, inputImage_cmplHH, inputFirstDimIsOdd);
  return outputFreqCoeffs;
}

//Produce FFT coefficients from inputImage, but return FFT coefficients that are based on
//the referenceImageBase (up/down sampled coefficients)

FloatImageType::Pointer At_fhp(HalfHermetianImageType::Pointer inLRCoeffs,
                               const bool inputFirstDimIsOdd,
                               itk::ImageBase<3>::Pointer desiredOutputRef)
{
  //TODO: Review this FFTScaler value why sqrt?
  HalfHermetianImageType::Pointer outputFreqCoeffs = ReshapeFFT(desiredOutputRef, inLRCoeffs,inputFirstDimIsOdd );
  const PrecisionType FFTScaler = std::sqrt(
    static_cast<PrecisionType>(desiredOutputRef->GetLargestPossibleRegion().GetNumberOfPixels()));
  const bool referenceImageBase_ActualXDimensionIsOdd = desiredOutputRef->GetLargestPossibleRegion().GetSize()[0] %
                                                             2 == 1;

  FloatImageType::Pointer outHRRealImage = GetInverseFFT(outputFreqCoeffs,
                                                         referenceImageBase_ActualXDimensionIsOdd, FFTScaler);
  return outHRRealImage;
}

// Apply an identity transform to the image
// and rescale image by manipulating the FFT
// coeffiecients into the new shape.
FloatImageType::Pointer IdentityResampleByFFT(FloatImageType::Pointer inOriginalImage,
                                              itk::ImageBase<3>::Pointer desiredOutputRef)
{
  HalfHermetianImageType::Pointer inputImage_cmplHH = GetForwardFFT(inOriginalImage, 1.0F);
  const PrecisionType UpsampleFFTScaler = static_cast<PrecisionType>(desiredOutputRef->GetLargestPossibleRegion().GetNumberOfPixels() )
                                  / static_cast<PrecisionType>(inOriginalImage->GetLargestPossibleRegion().GetNumberOfPixels());
  HalfHermetianImageType::Pointer outputFreqCoeffs = ReshapeFFT(desiredOutputRef, inputImage_cmplHH,
                                                                inOriginalImage->GetLargestPossibleRegion().GetSize()[0] %2 == 1);

  const bool referenceImageBase_ActualXDimensionIsOdd = desiredOutputRef->GetLargestPossibleRegion().GetSize()[0] % 2 == 1;

  FloatImageType::Pointer outHRRealImage = GetInverseFFT(outputFreqCoeffs,
                                                         referenceImageBase_ActualXDimensionIsOdd, UpsampleFFTScaler);
  return outHRRealImage;
}


//AT_fhp
HalfHermetianImageType::Pointer GetLowPassFilterFFT(FloatImageType::Pointer inputImage,
                                                    itk::ImageBase<3>::Pointer referenceImageBase)
{
  FloatFFTWFullFFTType::Pointer fft = FloatFFTWFullFFTType::New();
  fft->SetInput(inputImage);
  fft->Update();
  HalfHermetianImageType::Pointer inputImage_cmplHH = fft->GetOutput();
  HalfHermetianImageType::Pointer smallHH = CreateZeroFFTCoefficients(referenceImageBase);
  MoveFFTCoeffs(smallHH,
                referenceImageBase->GetLargestPossibleRegion().GetSize()[0]%2 ==1,
                inputImage_cmplHH,
                inputImage->GetLargestPossibleRegion().GetSize()[0] %2 ==1);
  MoveFFTCoeffs(inputImage_cmplHH,
                inputImage->GetLargestPossibleRegion().GetSize()[0] %2 ==1,
                smallHH,
                referenceImageBase->GetLargestPossibleRegion().GetSize()[0]%2 ==1 );
  //TODO: This could probably be more efficient.
  //============================
  inputImage_cmplHH->SetSpacing(inputImage->GetSpacing());
  inputImage_cmplHH->SetOrigin(inputImage->GetOrigin());
  inputImage_cmplHH->SetDirection(inputImage->GetDirection());
  return inputImage_cmplHH;
}


CVImageType::Pointer GetGradient(FloatImageType::Pointer inputImage)
{
#if 1
  GradientType::Pointer gradient_of_p = GradientType::New();
  typedef itk::PeriodicBoundaryCondition<FloatImageType> FloatBoundaryType;
  gradient_of_p->OverrideBoundaryCondition( new FloatBoundaryType );
  gradient_of_p->SetInput(inputImage);
  gradient_of_p->SetUseImageDirection(false);
  gradient_of_p->SetUseImageSpacingOff();
  gradient_of_p->Update();
  return gradient_of_p->GetOutput();
#else
  CVImageType::Pointer cvImage = CreateEmptyImage<CVImageType>(inputImage.GetPointer());
  itk::ImageRegionIteratorWithIndex<FloatImageType> inIter(inputImage,inputImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<CVImageType>  outIter(cvImage,cvImage->GetLargestPossibleRegion());
  const FloatImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
  while(!inIter.IsAtEnd())
  {
    const FloatImageType::IndexType idx = inIter.GetIndex();
    FloatImageType::IndexType NextX = idx;
    NextX[0] = (NextX[0] +1) % size[0];
    FloatImageType::IndexType NextY = idx;
    NextY[1] = (NextY[1] +1) % size[1];
    FloatImageType::IndexType NextZ = idx;
    NextZ[2] = (NextZ[2] +1) % size[2];
    CVImageType::PixelType & tmp = outIter.Value();
    tmp[0]=(inputImage->GetPixel(NextX) - inIter.Value());
    tmp[1]=(inputImage->GetPixel(NextY) - inIter.Value());
    tmp[2]=(inputImage->GetPixel(NextZ) - inIter.Value());
    ++inIter;
    ++outIter;
  }
  return cvImage;
#endif
}

DivergenceType::OutputImageType::Pointer GetDivergence(CVImageType::Pointer inputImage)
{
  DivergenceType::Pointer divergence_of_gradient_of_p = DivergenceType::New();
  typedef itk::PeriodicBoundaryCondition<CVImageType> CVBoundaryType;
  divergence_of_gradient_of_p->OverrideBoundaryCondition(new CVBoundaryType);
  divergence_of_gradient_of_p->SetInput(inputImage);
  divergence_of_gradient_of_p->SetUseImageSpacingOff();
  divergence_of_gradient_of_p->Update();
#if 1
  //To make compatible with Matlab, need to negate
  DivergenceType::OutputImageType::Pointer negativeMatlabDiv = divergence_of_gradient_of_p->GetOutput();
  DivergenceType::OutputImageType::Pointer returnImage = opIC(negativeMatlabDiv,negativeMatlabDiv,'*',-1.0F);
  return returnImage;
#else
  return divergence_of_gradient_of_p->GetOutput();
#endif
}


void WriteComplexImages(HalfHermetianImageType::Pointer cmplxIn, const std::string prefix)
{
  typedef itk::ComplexToRealImageFilter<HalfHermetianImageType,FloatImageType> C2RType;
  C2RType::Pointer c2r = C2RType::New();
  c2r->SetInput(cmplxIn);
  c2r->Update();
  WriteFile(c2r->GetOutput(),prefix+"_real.nii");
#if 1
  typedef itk::ComplexToImaginaryImageFilter<HalfHermetianImageType,FloatImageType> C2IType;
  C2IType::Pointer c2i = C2IType::New();
  c2i->SetInput(cmplxIn);
  c2i->Update();
  c2i->GetOutput();
  WriteFile(c2i->GetOutput(),prefix+"_imag.nii");
#endif
}

HalfHermetianImageType::Pointer  ReadComplexImages(const std::string prefix)
{
  typedef itk::ImageFileReader<FloatImageType> FileReaderType;
  FileReaderType::Pointer fr = FileReaderType::New();
  fr->SetFileName(prefix+"_real.nii");
  fr->Update();
  FloatImageType::Pointer re = fr->GetOutput();
  FileReaderType::Pointer fi = FileReaderType::New();
  fi->SetFileName(prefix+"_imag.nii");
  fi->Update();
  FloatImageType::Pointer im = fi->GetOutput();

  HalfHermetianImageType::Pointer cmplx = CreateEmptyImage<HalfHermetianImageType>(re.GetPointer());
  itk::ImageRegionIterator<HalfHermetianImageType> cmplxIter(cmplx,cmplx->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> reIter(re,re->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> imIter(im,im->GetLargestPossibleRegion());
  while(!reIter.IsAtEnd())
  {
    cmplxIter.Set(std::complex<PrecisionType>(reIter.Value(),imIter.Value()));
    ++cmplxIter;
    ++reIter;
    ++imIter;
  }
  return cmplx;
}

#include <itkRescaleIntensityImageFilter.h>
FloatImageType::Pointer NormalizeDataComponent(FloatImageType::Pointer arr)
{
  typedef itk::RescaleIntensityImageFilter<FloatImageType,FloatImageType> RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput(arr);
  rescaler->SetOutputMaximum(1.0);
  rescaler->SetOutputMinimum(0.0);
  rescaler->Update();

  return rescaler->GetOutput();
}
