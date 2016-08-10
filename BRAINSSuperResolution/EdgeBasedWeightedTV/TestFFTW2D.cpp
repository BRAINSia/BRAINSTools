//
// Created by Johnson, Hans J on 7/26/16.
//

//HACK Test why this needs to be here.
#include <itkComplexToComplexFFTImageFilter.h>

#include <itkFFTWComplexToComplexFFTImageFilter.h>

#include <itkFFTWForwardFFTImageFilter.h>
#include <itkFFTWInverseFFTImageFilter.h>

#include <itkFFTWRealToHalfHermitianForwardFFTImageFilter.h>
#include <itkFFTWHalfHermitianToRealInverseFFTImageFilter.h>

#include <itkImageRegionIteratorWithIndex.h>

#include <iostream>


template<typename IMTYPE>
void DumpImage2D(IMTYPE out, const std::string ImageName)
{
  const typename IMTYPE::ObjectType::SizeType size = out->GetLargestPossibleRegion().GetSize();
  typename IMTYPE::ObjectType::IndexType idx;
  std::cout << "======= " << ImageName << " ===================================" << std::endl;
  std::cout.precision(2);
  std::cout.width(8);
  for(size_t j = 0; j < size[1]; ++j)
  {
    idx[1]=j;
    for(size_t i = 0; i< size[0]; ++i)
    {
      idx[0]=i;
      std::cout.width(16);
      std::cout.fill(' ');
      std::cout << std::right << std::fixed << out->GetPixel(idx);
    }
    std::cout << "  y=" << j << std::endl;
  }
  std::cout << "==========================================" << std::endl;
}

int main(int argc, char * argv[])
{
  typedef double ScalarType;
  typedef std::complex<ScalarType> ComplexType;
  const size_t DIMENSION =2;

  typedef itk::Image<ComplexType,DIMENSION> CImageType;
  typedef itk::Image<ScalarType ,DIMENSION> SImageType;

  typedef itk::FFTWComplexToComplexFFTImageFilter<CImageType> FFTC2C;
  typedef itk::FFTWForwardFFTImageFilter<SImageType>          FFFTR2C;
  typedef itk::FFTWInverseFFTImageFilter<CImageType>          IFFTC2R;

  CImageType::Pointer CR = CImageType::New();
  {
    CImageType::RegionType region;
    CImageType::SizeType   size;
    size.Fill(9);
    region.SetSize(size);
    CR->SetRegions(region);
  }
  CR->Allocate();
  CR->FillBuffer(itk::NumericTraits<typename CImageType::PixelType>::ZeroValue());
  itk::ImageRegionIteratorWithIndex<CImageType> cIter(CR,CR->GetLargestPossibleRegion());
  while(!cIter.IsAtEnd())
  {
    for(int dim = 0 ; dim < CImageType::ImageDimension; ++dim)
    {
      if( cIter.GetIndex()[dim] < (CR->GetLargestPossibleRegion().GetSize()[dim]*1)/3
        || cIter.GetIndex()[dim] > (CR->GetLargestPossibleRegion().GetSize()[dim]*2)/3 )
      {
        cIter.Set(ComplexType(100.0,0.0));
      }
    }
    ++cIter;
  }
  DumpImage2D(CR,"CR");

  FFTC2C::Pointer forw_c2c = FFTC2C::New();
  forw_c2c->SetInput(CR);
  forw_c2c->SetTransformDirection(FFTC2C::FORWARD);
  forw_c2c->Update();
  CImageType::Pointer forw_c2c_image = forw_c2c->GetOutput();
  DumpImage2D(forw_c2c_image,"forw_c2c_image");

  FFTC2C::Pointer inv_c2c = FFTC2C::New();
  inv_c2c->SetInput(forw_c2c_image);
  inv_c2c->SetTransformDirection(FFTC2C::INVERSE);
  inv_c2c->Update();
  CImageType::Pointer inv_c2c_image = inv_c2c->GetOutput();
  DumpImage2D(inv_c2c_image,"inv_c2c_image");

  IFFTC2R::Pointer ifft_c2r = IFFTC2R::New();
  ifft_c2r->SetInput(forw_c2c_image);
  ifft_c2r->Update();
  SImageType::Pointer ifft_c2r_image = ifft_c2r->GetOutput();
  DumpImage2D(ifft_c2r_image,"ifft_c2r_image");


  return 0;
}