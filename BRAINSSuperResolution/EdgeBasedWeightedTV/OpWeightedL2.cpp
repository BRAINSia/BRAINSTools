#include "SRTypes.h"
#include "FFTWUpsample.h"

#include <cblas.h>

#include <itkTimeProbe.h>

#include "MathUtils.h"


#include <itkSqrtImageFilter.h>
#include <itkComplexToModulusImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkBinaryFunctorImageFilter.h>

//Special override for CVImageType
PrecisionType * GetFirstPointer(CVImageType::Pointer in)
{
  CVType * firstCovariantVectorX = in->GetBufferPointer();
  return  firstCovariantVectorX->GetDataPointer();
}

PrecisionType * GetFirstPointer(HalfHermetianImageType::Pointer in)
{
  return reinterpret_cast<PrecisionType *>(in->GetBufferPointer());
}


template<typename ImagePointerType>
PrecisionType * GetFirstPointer(ImagePointerType in)
{
  return in->GetBufferPointer();
}


//Implement out = c*(a*x + y), y is output, and is corrupted
template<typename ImagePointerType>
void AddAllElements(ImagePointerType &OutImg,
                           const PrecisionType aScaler,
                           ImagePointerType &xImg, ImagePointerType &yImg,
                           const PrecisionType cScaler= 1.0F)
{
  const size_t N = xImg->GetLargestPossibleRegion().GetNumberOfPixels()*xImg->GetNumberOfComponentsPerPixel();
  PrecisionType * x =  GetFirstPointer(xImg);
  PrecisionType * y =  GetFirstPointer(yImg);
  cblas_saxpy(N,aScaler,x,1,y,1);

  PrecisionType * Out=  GetFirstPointer(OutImg);
  if(cScaler != 1.0F)
  {
  cblas_sscal(N,cScaler,y,1);
  }
  if( OutImg.GetPointer() != yImg.GetPointer())
  {
    ImagePointerType temp = OutImg;
    OutImg = yImg;
    yImg = temp; //Sanity Check to induce failures if variable is needed in future.
    //yImg has been corrupted by processing here.
  }
}

//Implement out = x*y, y, and y is corrupted
template<typename ImagePointerType>
void MultiplyVectors(ImagePointerType &OutImg,
                           ImagePointerType &xImg, ImagePointerType &yImg)
{
  const size_t N = xImg->GetLargestPossibleRegion().GetNumberOfPixels()*xImg->GetNumberOfComponentsPerPixel();
  PrecisionType * x_Start = GetFirstPointer(xImg);
  const PrecisionType * x_End = x_Start+N;
  const PrecisionType * y = GetFirstPointer(yImg);
  PrecisionType * o = GetFirstPointer(OutImg);
  //HACK ADD open_MP here
  //HACK end = x+N : while x < end
  for(PrecisionType *x = x_Start; x < x_End; ++x)
  {
    (*o) = (*y) * (*x);
    ++o;
    ++y;
  }
}

template <typename ImageTypePointer>
void Duplicate(ImageTypePointer & Y, ImageTypePointer &YminusL)
{
  const size_t N= Y->GetLargestPossibleRegion().GetNumberOfPixels() *Y->GetNumberOfComponentsPerPixel();
  const PrecisionType * firstInput = GetFirstPointer(Y);
  const PrecisionType * lastInput = firstInput + N;
  PrecisionType * firstOutput = GetFirstPointer(YminusL);
  std::copy(firstInput,lastInput,firstOutput);
}


#include <itkVectorMagnitudeImageFilter.h>

FloatImageType::Pointer ComputeSqrtMu(FloatImageType::Pointer mu)
{
  typedef itk::SqrtImageFilter<FloatImageType,FloatImageType> SqrtType;
  SqrtType::Pointer sqrtFilter = SqrtType::New();
  sqrtFilter->SetInput(mu);
  sqrtFilter->Update();
  return sqrtFilter->GetOutput();
}

// (2*mu+gam)^{-1}
static CVImageType::Pointer ComputeInvTwoMuPlusGamma( FloatImageType::Pointer edgemask,  const PrecisionType gam)
{
  CVImageType::Pointer repMu = CreateEmptyImage<CVImageType>(edgemask);
  itk::ImageRegionIterator<CVImageType> cvIt(repMu,repMu->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> muIt(edgemask,edgemask->GetLargestPossibleRegion());
  CVImageType::PixelType temp;
  while(!muIt.IsAtEnd())
  {
    cvIt.Set( 1.0 / ( 2.0*muIt.Value() + gam) );
    ++cvIt;
    ++muIt;
  }
  return repMu;
}

static HalfHermetianImageType::Pointer GetAFP_of_b(FloatImageType::Pointer norm01_lowres, FloatImageType::Pointer edgemask)
{
  FloatImageType::Pointer upsampledB = IdentityResampleByFFT(norm01_lowres, edgemask.GetPointer());
  HalfHermetianImageType::Pointer b_FC = A_fhp(upsampledB, norm01_lowres.GetPointer());
  return b_FC;
}

static FloatImageType::Pointer GetDiracDeltaImage(FloatImageType::Pointer edgemask)
{
  FloatImageType::Pointer p_image = CreateEmptyImage<FloatImageType>(edgemask);
  FloatImageType::IndexType zeroIdx;
  zeroIdx.Fill(0);
  p_image->SetPixel(zeroIdx,1);
  return p_image;
}

static HalfHermetianImageType::Pointer GetLowpassOperator(FloatImageType::Pointer norm01_lowres,
                                                           FloatImageType::Pointer p_image, const PrecisionType scaler)
{
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

  HalfHermetianImageType::Pointer testA_fhp = A_fhp(p_image,norm01_lowres.GetPointer());
  FloatImageType::Pointer testAtA = At_fhp(testA_fhp,norm01_lowres->GetLargestPossibleRegion().GetSize()[0]%2 == 1, p_image.GetPointer());
  HalfHermetianImageType::Pointer AtAhat = GetForwardFFT(testAtA);
  return opIC(AtAhat,AtAhat,'*',scaler); // A is the linear measurement operator
}

static FloatImageType::Pointer  MakeTwoAtb(FloatImageType::Pointer Atb)
{
FloatImageType::Pointer TwoAtb =DeepImageCopy<FloatImageType>(Atb);
return opIC(TwoAtb,TwoAtb,'*',2.0);
}

/*
OPWEIGHTEDL2: Solves weighted L2 regularized inverse problems.
Minimizes the cost function
X* = argmin_X ||A(X)-b_FC||_2^2 + lambda ||W |D(X)| ||_2^2
where     X* = recovered image
          A  = linear measurement operator
          b_FC  = (noisy) measurements
                 %           W  = diagonal weight matrix built from the edge mask
          |D(X)| = gradient magnitude at each pixel

Inputs:  A = function handle representing the forward
              model/measurement operator
         At = function handle representing the backwards model/
              the transpose of the measurment operator.
              (e.g. if A is a downsampling, At is a upsampling)
         b_FC =  a vector of measurements; should match the
              dimensions of A(X)
         lambda = regularization parameter that balances data fidelity
              and smoothness. set lambda high for more smoothing.
         siz = output image size, e.g. siz = [512,512]
         Niter = is the number of iterations; should be ~100-500

Output:  X = high-resolution output image
         cost = array of cost function value vs. iteration
Define AtA fourier mask
PrecisionType lambda, uvec ind_samples, frowvec res, int Niter, double tol, PrecisionType gam, FloatImageType::Pointer& X, frowvec& cost, frowvec& resvec)
 */
FloatImageType::Pointer OpWeightedL2(FloatImageType::Pointer norm01_lowres, FloatImageType::Pointer edgemask)
{
  const PrecisionType lambda = 1e-3F ;
  const int Niter = 100 ;
  const PrecisionType tol = 1e-8F ;

  const PrecisionType gam = 1.0F ;

  //typedef itk::VectorMagnitudeImageFilter<CVImageType, FloatImageType> GMType;

  //The optimal filter for modeling the measurement operator is low pass filter in this case
  // NOTE: That the A operator is a projection operator, so A^{T}A = A, That is to say that applying
  //       the A^{T} to A results in A.
  FloatImageType::Pointer p_image = GetDiracDeltaImage(edgemask);
  // Precompute

  //Make high-res coefficients
  const HalfHermetianImageType::Pointer b_FC = GetAFP_of_b(norm01_lowres, edgemask);
  //TODO: too many copies of Atb here.
  FloatImageType::Pointer Atb = At_fhp(b_FC,
                                       edgemask->GetLargestPossibleRegion().GetSize()[0]%2 == 1,
                                       edgemask.GetPointer());
  FloatImageType::Pointer TwoAtb = MakeTwoAtb(Atb);
  FloatImageType::Pointer X = DeepImageCopy<FloatImageType>(Atb);
  Atb = nullptr; //Save memory here

  CVImageType::Pointer DX = GetGradient(X);
  CVImageType::Pointer L = CreateEmptyImage<CVImageType>(DX);
  CVImageType::Pointer Y = CreateEmptyImage<CVImageType>(DX);
  //CVImageType::Pointer WDX = CreateEmptyImage<CVImageType>(DX);
  CVImageType::Pointer residue = CreateEmptyImage<CVImageType>(DX);
  CVImageType::Pointer YminusL = CreateEmptyImage<CVImageType>(DX);
  FloatImageType::Pointer tempValue=CreateEmptyImage<FloatImageType>(DX);

  std::vector<PrecisionType> resvec(Niter,0);
  std::vector<PrecisionType> cost(Niter,0);

#ifdef USE_WRITE_DEGUBBING
  itk::ComplexToModulusImageFilter<HalfHermetianImageType,FloatImageType>::Pointer cpx2abs =
  itk::ComplexToModulusImageFilter<HalfHermetianImageType,FloatImageType>::New();
#endif

  CVImageType::Pointer gradIm = GetGradient(p_image);
  FloatImageType::Pointer divIm = GetDivergence(gradIm);
  HalfHermetianImageType::Pointer DtDhat = GetForwardFFT(divIm);
  // TODO:  ALL SAME TO HERE!
  typedef HalfHermetianImageType::PixelType FCType;
  HalfHermetianImageType::Pointer TwoTimesAtAhatPlusLamGamDtDhat = CreateEmptyImage<HalfHermetianImageType>(DtDhat);
  {
    HalfHermetianImageType::Pointer TwoTimesAtAhat = GetLowpassOperator(norm01_lowres,p_image, 2.0F);
    TwoTimesAtAhatPlusLamGamDtDhat = opIC(TwoTimesAtAhatPlusLamGamDtDhat,FCType(lambda*gam),'*',DtDhat);
    //TODO:  Make This an inverse!
    TwoTimesAtAhatPlusLamGamDtDhat = opII(TwoTimesAtAhatPlusLamGamDtDhat,TwoTimesAtAhat,'+',TwoTimesAtAhatPlusLamGamDtDhat);
  }
  p_image = nullptr; //Save memory

  const bool edgemask_ActualXDimensionIsOdd = edgemask->GetLargestPossibleRegion().GetSize()[0] % 2 == 1;

  CVImageType::Pointer InvTwoMuPlusGamma = ComputeInvTwoMuPlusGamma(edgemask,gam);
  //FloatImageType::Pointer SqrtMu = ComputeSqrtMu(edgemask);

#define USE_BLAS_WRAPPERS
#ifdef USE_BLAS_WRAPPERS
#else
  typedef itk::AddImageFilter<CVImageType,CVImageType> CVImageAdder;
  CVImageAdder::Pointer dxPlusL = CVImageAdder::New();
#endif

  itk::TimeProbe tp;
  tp.Start();

  HalfHermetianImageType::Pointer tempRatioFC = CreateEmptyImage<HalfHermetianImageType>(DtDhat);

  for (size_t i=0; i < Niter; ++i)
  {
    std::cout << "Iteration : " << i << std::endl;

#ifdef USE_BLAS_WRAPPERS
    //Z = 1.0*L+DX
    AddAllElements(DX,1.0F,L,DX,gam);//DX destroyed
    CVImageType::Pointer & Z = DX;
#else
    //Z = opII(Z,DX,'+',L);
    dxPlusL->SetInput1(DX);
    dxPlusL->SetInput2(L);
    dxPlusL->SetInPlace(true);
    dxPlusL->Update();
    CVImageType::Pointer Z=dxPlusL->GetOutput();
    MultiplyCVByScalar(Z,gam);
#endif
#ifdef USE_BLAS_WRAPPERS
    //Y=InvTwoMuPlusGamm.*Z
    MultiplyVectors(Y,InvTwoMuPlusGamma,Z);
#else
    Y = opII(Y,Z,'*',InvTwoMuPlusGamma);
#endif

    // X Subprob
    // Numerator = 2*Atb+lambda*gam*SRdiv(Y-L))
#ifdef USE_BLAS_WRAPPERS
    //YminusL = 1.0F* SRdiv( -1.0F*L + Y)
    Duplicate(Y,YminusL);
    AddAllElements(YminusL,-1.0F,L,YminusL,1.0F);
#else
    YminusL=opII(YminusL,Y,'-',L);
#endif
    FloatImageType::Pointer tempNumerator=GetDivergence(YminusL);
#ifdef USE_BLAS_WRAPPERS
    //lambd*gam*tempNumerator+TwoAtb
    Duplicate(TwoAtb,tempValue);
    AddAllElements(tempValue,lambda*gam,tempNumerator,tempValue,1.0F);
    HalfHermetianImageType::Pointer tempNumeratorFC = GetForwardFFT(tempValue);
#else
    tempNumerator=opIC(tempNumerator,lambda*gam,'*',tempNumerator);
    tempNumerator=opII(tempNumerator,TwoAtb,'+',tempNumerator);
    HalfHermetianImageType::Pointer tempNumeratorFC = GetForwardFFT(tempNumerator);
#endif

    //KEEP
    tempRatioFC = opII_scalar(tempRatioFC,tempNumeratorFC,'/',TwoTimesAtAhatPlusLamGamDtDhat);

    X = GetInverseFFT(tempRatioFC,edgemask_ActualXDimensionIsOdd,1.0); //TODO: Determine scale factor here. X

    // should be on same dynamic range as b
    DX = GetGradient(X);
    residue = opII(residue, DX, '-', Y); //TODO:  Implement math graph output here
    L = opII(L,L,'+',residue);

    //
    resvec[i] = 0; //TODO: Figure out the math for here
    if ( i > 900000 )
    {
      tp.Stop();
      std::cout << " Only iterations " << tp.GetTotal() << tp.GetUnit() << std::endl;
      return X;
    }
    if( i > 99 ) //HACK: Cutting this function short
    {
      return X;
    }
    //WDX = opII_CVmult(WDX,SqrtMu,'*',DX);
    //diff = opII(diff,A_fhp(X,norm01_lowres.GetPointer()),'-',b_FC);
    //
    //cost[i] = 0; //TODO: Need to figure out math for here
  }
  return X;
}
