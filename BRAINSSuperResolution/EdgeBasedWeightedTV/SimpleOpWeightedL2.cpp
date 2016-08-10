

#include "OpWeightedL2.h"
#include "FFTWUpsample.h"
#include <itkTimeProbe.h>



static sitk::Image GetDiracDeltaImage(sitk::Image edgemask) {
  sitk::Image p_image = sitk::Image(edgemask.GetSize(), sitk::sitkFloat32);
  p_image.CopyInformation(edgemask);
  std::vector<uint32_t> idx(3, 0);
  p_image.SetPixelAsFloat(idx, 1.0F);
  return p_image;
}

//Create an appropriate memory layout for holding the output of an FFT based on
//the size of the referenceImageBase information.
sitk::Image CreateZeroFFTCoefficients(sitk::Image &referenceImageBase) {
  sitk::Image outputFreqCoeffs = sitk::Image(referenceImageBase.GetSize(), sitk::sitkComplexFloat32);
  outputFreqCoeffs.CopyInformation(referenceImageBase);
#ifdef USE_HALF_FFTW
  {
    XXX HACK: THINK CAREFULLY TO REVIEW HERE
    std::vector<int32_t> temp_size = referenceImageBase.GetSize();
    temp_size[0] = (temp_size[0]/2 ) + 1; //Only store 1/2 of the first dimension to take advantage of symmetry
    outputFreqCoeffs.SetSize(temp_size);
  }
#endif
  return outputFreqCoeffs;
}


template<typename ITKImageType>
ITKImageType *ConvertToITK(sitk::Image &in) {
  typename ITKImageType::Pointer outITK = dynamic_cast<ITKImageType *>( in.GetITKBase());
  if (outITK.IsNull()) {
    std::cerr << "ERROR: Can not convert from " << in.GetPixelIDTypeAsString() << "\n-- to --\n" << outITK << std::endl;
    exit(-1);
  }
  return outITK.GetPointer();
}

// Transfer FFT coefficients to/from different sized frequency domains (i.e. upsample/downsample image)
// downsampling has the effect of low-pass filtering.
// implement out = p(ind_samples)
// FFTScalar is a scale factor applied uniformly during the reshape process
sitk::Image ReshapeFFT(
    sitk::Image outputRealImageBase,
    sitk::Image inputFreqCoeffs,
    const bool inFirstSpatialDeminsionIsOdd) {
  const bool FirstDimensionIsOdd = outputRealImageBase.GetSize()[0] % 2 == 1;
  sitk::Image outputFreqCoeffs = CreateZeroFFTCoefficients(outputRealImageBase);

  /* NEEED ITK ITERATORS HERE */
  HalfHermetianImageType::Pointer itkoutputFreqCoeffs = ConvertToITK<HalfHermetianImageType>(outputFreqCoeffs);
  HalfHermetianImageType::Pointer itkinputFreqCoeffs = ConvertToITK<HalfHermetianImageType>(inputFreqCoeffs);
  typedef itk::ImageRegionIterator<HalfHermetianImageType> cmplxHHIteratorType;
  cmplxHHIteratorType outputHHIter(itkoutputFreqCoeffs, itkoutputFreqCoeffs->GetLargestPossibleRegion());

  MoveFFTCoeffs(itkoutputFreqCoeffs, FirstDimensionIsOdd, itkinputFreqCoeffs, inFirstSpatialDeminsionIsOdd);
  return outputFreqCoeffs;
}

//Produce FFT coefficients from inputImage, but return FFT coefficients that are based on
//the referenceImageBase (up/down sampled coefficients) defines the desired output FFT Space
// low_resolution samples
// inHRRealImage is always a high resolution image
// lpfSamplesDef is always a low-resolution example
// referenceImageBase is always a high-resolution example
sitk::Image A_fhp(sitk::Image inHRRealImage,
                  sitk::Image desiredOutputRef) {
  const bool inputFirstDimIsOdd = inHRRealImage.GetSize()[0] % 2 == 1;
  const PrecisionType FFTScaler =
      1.0F / std::sqrt(static_cast<PrecisionType>(inHRRealImage.GetNumberOfPixels()));

  sitk::Image inputImage_cmplHH = FFTScaler * sitk::ForwardFFT(inHRRealImage);

  //==========================
  // NEED ITK HERE
  // Transfer FFT coeficients to new space
  //typedef itk::ImageRegionIterator<HalfHermetianImageType> cmplxHHIteratorType;
  sitk::Image outputFreqCoeffs = ReshapeFFT(desiredOutputRef, inputImage_cmplHH, inputFirstDimIsOdd);
  return outputFreqCoeffs;
}

// Apply an identity transform to the image
// and rescale image by manipulating the FFT
// coeffiecients into the new shape.
sitk::Image IdentityResampleByFFT(sitk::Image inOriginalImage,
                                  sitk::Image desiredOutputRef) {
  sitk::Image inputImage_cmplHH = ForwardFFT(inOriginalImage);
  const PrecisionType UpsampleFFTScaler = static_cast<PrecisionType>(desiredOutputRef.GetNumberOfPixels())
                                          / static_cast<PrecisionType>(inOriginalImage.GetNumberOfPixels());
  sitk::Image outputFreqCoeffs = ReshapeFFT(desiredOutputRef, inputImage_cmplHH,
                                            inOriginalImage.GetSize()[0] % 2 == 1);

  const bool referenceImageBase_ActualXDimensionIsOdd = desiredOutputRef.GetSize()[0] % 2 == 1;

  sitk::Image outHRRealImage = sitk::InverseFFT(outputFreqCoeffs) * UpsampleFFTScaler;
  return outHRRealImage;
}

static sitk::Image GetAFP_of_b(sitk::Image norm01_lowres, sitk::Image edgemask) {
  sitk::Image upsampledB = IdentityResampleByFFT(norm01_lowres, edgemask);
  sitk::Image b_FC = A_fhp(upsampledB, norm01_lowres);
  return b_FC;
}

//Produce FFT coefficients from inputImage, but return FFT coefficients that are based on
//the referenceImageBase (up/down sampled coefficients)

sitk::Image At_fhp(sitk::Image inLRCoeffs,
                   const bool inputFirstDimIsOdd,
                   sitk::Image desiredOutputRef) {
  //TODO: Review this FFTScaler value why sqrt?
  sitk::Image outputFreqCoeffs = ReshapeFFT(desiredOutputRef, inLRCoeffs, inputFirstDimIsOdd);
  const PrecisionType FFTScaler = std::sqrt(
      static_cast<PrecisionType>(desiredOutputRef.GetNumberOfPixels()));

  //TODO: HalfHermitian
  const bool referenceImageBase_ActualXDimensionIsOdd = desiredOutputRef.GetSize()[0] % 2 == 1;

  sitk::Image outHRRealImage = sitk::InverseFFT(outputFreqCoeffs) * FFTScaler;
  return outHRRealImage;
}

#include <itkPeriodicBoundaryCondition.h>

sitk::Image GetGradient(sitk::Image inputImage) {
#if 1
  GradientVarVecType::Pointer gradient_of_p = GradientVarVecType::New();
  typedef itk::PeriodicBoundaryCondition<FloatImageType> FloatBoundaryType;
  gradient_of_p->OverrideBoundaryCondition(new FloatBoundaryType);
  FloatImageType::Pointer itkinputImage = dynamic_cast<FloatImageType *>(inputImage.GetITKBase());
  gradient_of_p->SetInput(itkinputImage);
  gradient_of_p->SetUseImageSpacing(false);
  gradient_of_p->SetUseImageDirection(false);
  gradient_of_p->Update();
  VarVecImageType::Pointer outImage = gradient_of_p->GetOutput();
  sitk::Image returnValue = sitk::Image(outImage.GetPointer());
  return returnValue;
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


sitk::Image GetDivergence(sitk::Image inputImage) {
  VarVecImageType::Pointer itkinputImage = dynamic_cast<VarVecImageType *>(inputImage.GetITKBase());

  DivergenceVarVecType::Pointer divergence_of_gradient_of_p = DivergenceVarVecType::New();
  typedef itk::PeriodicBoundaryCondition<VarVecImageType> CVBoundaryType;
  divergence_of_gradient_of_p->OverrideBoundaryCondition(new CVBoundaryType);
  divergence_of_gradient_of_p->SetInput(itkinputImage);
  divergence_of_gradient_of_p->SetUseImageSpacingOff();
  divergence_of_gradient_of_p->Update();
  sitk::Image output = sitk::Image(divergence_of_gradient_of_p->GetOutput());
  //Negate to To make compatible with Matlab, need to negate
  return -1.0F * output;
}

static sitk::Image GetLowpassOperator(sitk::Image norm01_lowres, sitk::Image p_image, const PrecisionType scaler) {
/*
 * Make dirac-delta impulse response filter for convolving with the high-resolution image
The operator A should represent the characteristics of the measurement device as it
takes measurements at lower resolutions.  For example, in MRI as resolution is decreased,
larger anatomical regions are included in each measurement, thus giving a
smoother image as if the high resolution image were passed through a low-pass-filter.
p_image below represents the smoothing convolution filter that is applied to a
high-resolution image to best represent the low-resolution image sampling.
In the case of no downsampling p_image is simply a perfect dirac-delta function.
Convolution in the spatial-domain is multiplication in the frequency domain.

In the case were some downsampling is present between the lowresolution and
the highresolution images, we need the operator A to represent a low-pass filter
that optimally models the lower sampling rate.  The low-pass filter is accomplished
in the frequency domain by zeroing-out frequency components that are not representable
in the frequency domain of the low-resolution image.
 */

  sitk::Image testA_fhp = A_fhp(p_image, norm01_lowres);
  sitk::Image testAtA = At_fhp(testA_fhp, norm01_lowres.GetSize()[0] % 2 == 1, p_image);
  sitk::Image AtAhat = ForwardFFT(testAtA);
  return AtAhat * scaler; // A is the linear measurement operator
}

static sitk::Image ComputeInvTwoMuPlusGamma(sitk::Image edgemask, const PrecisionType gam) {
  sitk::Image mu = 1.0F / (2.0F * edgemask + gam);
  return mu;
}


sitk::Image MultiplyVectorByScalarImage(sitk::Image viSITK, sitk::Image siSITK) {
  sitk::Image outSITK(viSITK);
  typedef itk::VectorImage<PrecisionType, 3> VIType;
  VIType::Pointer out = dynamic_cast<VIType *>(outSITK.GetITKBase());
  typedef itk::Image<PrecisionType, 3> SIType;
  SIType::Pointer si = dynamic_cast<SIType *>(siSITK.GetITKBase());

  itk::ImageRegionIterator<VIType> outIt(out, out->GetLargestPossibleRegion());
  itk::ImageRegionIterator<SIType> siIt(si, si->GetLargestPossibleRegion());

  while (!outIt.IsAtEnd()) {
    outIt.Set(outIt.Get() * siIt.Value());
    ++outIt;
    ++siIt;
  }
  return outSITK;
}

void test() {
  std::vector<uint32_t> size(2, 2);
  sitk::Image vi(size, sitk::sitkVectorFloat32, 5);
  std::vector<uint32_t> idx(2, 1);
  const std::vector<float> oneVoxel = vi.GetPixelAsVectorFloat32(idx);
  std::cout << "XXXXX : " << oneVoxel[3] << std::endl;
}

using sitkHalfHermetianImageType =  sitk::Image;

sitk::Image SimpleOpWeightedL2(sitk::Image &norm01_lowres, sitk::Image &edgemask) {
  const PrecisionType lambda = 1e-3F;
  const int Niter = 100;
  const PrecisionType tol = 1e-8F;

  const PrecisionType gam = 1.0F;

  //The optimal filter for modeling the measurement operator is low pass filter in this case
  // NOTE: That the A operator is a projection operator, so A^{T}A = A, That is to say that applying
  //       the A^{T} to A results in A.
  sitk::Image p_image = GetDiracDeltaImage(edgemask);

  //Make high-res coefficients
  const sitkHalfHermetianImageType b_FC = GetAFP_of_b(norm01_lowres, edgemask);

  //TODO: too many copies of Atb here.
  sitk::Image X = At_fhp(b_FC,
                           edgemask.GetSize()[0] % 2 == 1,
                           edgemask);

  sitk::Image TwoAtb = X * 2.0;

  sitk::Image DX = GetGradient(X);

  sitk::Image DtDhat = ForwardFFT(GetDivergence(  GetGradient( p_image )));

  sitk::Image TwoTimesAtAhat = GetLowpassOperator(norm01_lowres, p_image, 2.0F);
  sitk::Image TwoTimesAtAhatPlusLamGamDtDhat = TwoTimesAtAhat + lambda * gam * DtDhat;

  sitk::Image InvTwoMuPlusGamma = ComputeInvTwoMuPlusGamma(edgemask, gam);

  //const bool edgemask_ActualXDimensionIsOdd = edgemask.GetSize()[0] % 2 == 1;

  sitk::Image SqrtMu = sitk::Sqrt(edgemask);

  std::vector<PrecisionType> resvec(Niter, 0);
  std::vector<PrecisionType> cost(Niter, 0);

  itk::TimeProbe tp;
  tp.Start();
  sitk::Image L = sitk::Image(DX.GetSize(), sitk::sitkVectorFloat32, 3); //Initialize L with all zeros
  L.CopyInformation(DX);

  PrecisionType ifftScaleFactor = DX.GetNumberOfPixels();
  //============================
  //NOTES:  For SimpleITK, math operations on images require both operands are same type!
  //        For SimpleITK & ITK, scalar values are first typecast to double, then typecast to image pixel type
  //         i.e.:  O = I * 7 where I is complex, and 7 is an integer
  //        foreach p in I:
  //             o = p * std::complex<float>(double(7))
  //
  //        VectorImages [*/] scalar in ITK are problematic for multiplication & division, and do not work!
  for (size_t i = 0; i < Niter; ++i) {
    std::cout << "Iteration : " << i << std::endl;
    //Z = 1.0*L+DX
    //Z = MultiplyVectorByConstant((DX + L), gam);
    sitk::Image Z = (DX + L) * gam;
    sitk::Image Y = MultiplyVectorByScalarImage(Z, InvTwoMuPlusGamma);

    itk::Image<float,3>::Pointer m = dynamic_cast<itk::Image<float,3> *>( InvTwoMuPlusGamma.GetITKBase());



    // X Subprob
    // Numerator = 2*Atb+lambda*gam*SRdiv(Y-L))
    sitk::Image NumeratorFC = sitk::ForwardFFT(TwoAtb + lambda * gam * GetDivergence(Y - L));
    sitk::Image tempRatioFC = NumeratorFC / TwoTimesAtAhatPlusLamGamDtDhat;

    X = sitk::InverseFFT(tempRatioFC);
    DX = GetGradient(X);

    sitk::Image residue = (DX - Y);
    //resvec[i] = NormOfRatioImages(residue,Y);

    //
    resvec[i] = 0; //TODO: Figure out the math for here
    if (i == 99) {
      tp.Stop();
      std::cout << " Only iterations " << tp.GetTotal() << tp.GetUnit() << std::endl;
      return X;
    }
    if (i > 99) //HACK: Cutting this function short
    {
      return X;
    }
    L += residue;
    //WDX = opII_CVmult(WDX,SqrtMu,'*',DX);
    //diff = opII(diff,A_fhp(X,norm01_lowres.GetPointer()),'-',b_FC);
    //
    //cost[i] = 0; //TODO: Need to figure out math for here
  }
  return X;
}
