//
// Created by Hans Johnson on 7/12/16.
//

// \author Hans J. Johnson
// \date 2016-07-10
// Test program for evaluating matlab to ITK conversions

#include <iostream>
#include <string>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkShrinkImageFilter.h>

#include "FFTWUpsample.h"

#include "itkPeriodicBoundaryCondition.h"

#include "itkDivergenceImageFilter.h"
//#include <itkGradientImageFilter.h>
#include <rtkForwardDifferenceGradientImageFilter.h>
#include <rtkBackwardDifferenceDivergenceImageFilter.h>

const size_t IM2DSIZE=3;

template<typename IMTYPE>
void DumpImage(IMTYPE out)
{
  const typename IMTYPE::ObjectType::SizeType size = out->GetLargestPossibleRegion().GetSize();
  typename IMTYPE::ObjectType::IndexType idx;
  std::cout << "==========================================" << std::endl;
  for(size_t k = 0; k < size[2]; ++k)
  {
    idx[2] = k;
  for(size_t j = 0; j < size[1]; ++j)
  {
    idx[1]=j;
    for(size_t i = 0; i< size[0]; ++i)
    {
      idx[0]=i;
      std::cout.width(16);
      std::cout.fill(' ');
      std::cout << std::right << out->GetPixel(idx);
    }
    std::cout << "  y=" << j << std::endl;
  }
  std::cout << "z= " << k << ". - . - . - . - . - . - . - . - . - . - . " << std::endl;
  }
  std::cout << "==========================================" << std::endl;
}

static void FillWithIndexValue(HalfHermetianImageType::Pointer img)
{
  const HalfHermetianImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
  HalfHermetianImageType::IndexType idx;

  for(idx[2]= 0; idx[2] < size[2]; ++idx[2]) {
    for(idx[1]= 0; idx[1] < size[1]; ++idx[1]) {
      for (idx[0] = 0; idx[0] < size[0]; ++idx[0]) {
        img->SetPixel(idx,std::complex<float>(idx[0],idx[1]));
      }
    }
  }
}

int main( int argc, char * argv[])
{
  const std::complex<float> Zero(0.8F,0.8F);
  HalfHermetianImageType::RegionType::IndexType startIndex;
  startIndex.Fill(0);
  HalfHermetianImageType::RegionType region;
  region.SetIndex(startIndex);
  HalfHermetianImageType::RegionType::SizeType size;
  size[0]=3;
  size[1]=4;
  size[2]=4;
  HalfHermetianImageType::Pointer inEvenSmall = HalfHermetianImageType::New();
  region.SetSize(size);
  inEvenSmall->SetRegions(region);
  inEvenSmall->Allocate();
  inEvenSmall->FillBuffer(Zero);
  FillWithIndexValue(inEvenSmall);
  //DumpImage(inEvenSmall);

  HalfHermetianImageType::Pointer inEvenLarge = HalfHermetianImageType::New();
  size[0]=4;
  size[1]=6;
  size[2]=6;
  region.SetSize(size);
  inEvenLarge->SetRegions(region);
  inEvenLarge->Allocate();
  inEvenLarge->FillBuffer(Zero);
  FillWithIndexValue(inEvenLarge);
  //DumpImage(inEvenLarge);
  HalfHermetianImageType::Pointer inEven2Large = HalfHermetianImageType::New();
  size[0]=4;
  size[1]=6;
  size[2]=6;
  region.SetSize(size);
  inEven2Large->SetRegions(region);
  inEven2Large->Allocate();
  inEven2Large->FillBuffer(Zero);
  FillWithIndexValue(inEven2Large);
  //DumpImage(inEvenLarge);

  const bool evenLargeFirstIsOdd  = true;
  const bool evenLarge2FirstIsOdd = true;
  const bool evenSmallFirstIsOdd  = true;

  size_t test_num = 1;
  switch ( test_num )
  {
    case 1:
    {
      DumpImage(inEvenSmall);
      MoveFFTCoeffs(inEvenLarge,
                    evenLargeFirstIsOdd,
                    inEvenSmall,
                    evenSmallFirstIsOdd
      );
      std::cout << "Modified Large" << std::endl;
      DumpImage(inEvenLarge);

      MoveFFTCoeffs(inEvenSmall,
                    evenSmallFirstIsOdd,
                    inEvenLarge,
      evenLargeFirstIsOdd);
      std::cout << "Modified Small" << std::endl;
      DumpImage(inEvenSmall);
    }
    break;
    case 2:
    {
      DumpImage(inEvenLarge);
      MoveFFTCoeffs(inEvenSmall,evenSmallFirstIsOdd, inEvenLarge,evenLargeFirstIsOdd );
      std::cout << "Modified Small" << std::endl;
      DumpImage(inEvenSmall);

      MoveFFTCoeffs(inEvenLarge,evenLargeFirstIsOdd, inEvenSmall, evenSmallFirstIsOdd );
      std::cout << "Modified Large" << std::endl;
      DumpImage(inEvenLarge);
    }
    break;
    case 3:
    {
      DumpImage(inEvenLarge);
      MoveFFTCoeffs(inEven2Large, evenLarge2FirstIsOdd, inEvenLarge, evenLargeFirstIsOdd );
      std::cout << "Modified Small" << std::endl;
      DumpImage(inEven2Large);

      MoveFFTCoeffs(inEvenLarge, evenLargeFirstIsOdd, inEven2Large, evenLarge2FirstIsOdd );
      std::cout << "Modified Large" << std::endl;
      DumpImage(inEvenLarge);
    }
    break;
    default:
      break;
  }

  //=================================================
#if 0
  typedef itk::Image<float,3> Image2D;
  Image2D::Pointer im2 = Image2D::New();
  Image2D::RegionType region;
  Image2D::RegionType::SizeType size;
  size.Fill(IM2DSIZE);
  region.SetSize(size);
  Image2D::SpacingType spacing;
  spacing.Fill(1.25);
  im2->SetSpacing(spacing);
  im2->SetRegions(region);
  im2->Allocate();
  Image2D::IndexType idx;
  idx.Fill(IM2DSIZE/2);
  idx.Fill(0);
  im2->SetPixel(idx,1.0F);
  //typedef itk::GradientImageFilter<Image2D> GradientType;
  typedef rtk::ForwardDifferenceGradientImageFilter<Image2D> GradientType;
  GradientType::Pointer grad = GradientType::New();
  grad->SetInput(im2);
  grad ->UseImageSpacingOff();
  typedef itk::PeriodicBoundaryCondition<Image2D> FloatBoundaryType;
  grad->OverrideBoundaryCondition( new FloatBoundaryType );
  grad->Update();
  GradientType::OutputImageType::Pointer gradIm = grad->GetOutput();
  //typedef itk::DivergenceImageFilter<GradientType::OutputImageType> DivergenceType;
  typedef rtk::BackwardDifferenceDivergenceImageFilter<GradientType::OutputImageType> DivergenceType;
  DivergenceType::Pointer div = DivergenceType::New();
  div->SetInput(gradIm);
  div->SetUseImageSpacing(false);
  typedef itk::PeriodicBoundaryCondition<GradientType::OutputImageType> CVBoundaryType;
  div->OverrideBoundaryCondition(new CVBoundaryType);
  div->Update();
  DivergenceType::OutputImageType::Pointer divIm = div->GetOutput();
//  DivergenceType::OutputImageType::Pointer divIm = GetDivergence(gradIm);
  std::cout << "Intensity Image" << std::endl;
  DumpImage(im2);
  DumpImage(gradIm);
  DumpImage(divIm);

  std::cout << "\n\n\n" << std::endl;
  DumpImage(GetGradient(im2.GetPointer()));
  DumpImage(GetDivergence(GetGradient(im2.GetPointer())));
#endif
#if 0

  // testIndexConversion();
  // return EXIT_SUCCESS;
  typedef itk::ImageFileReader<FloatImageType> ReaderType;
  typedef itk::ImageFileWriter<FloatImageType> WriterType;

  const std::string hriFileName = "/Users/johnsonhj/src/20160711_SuperResolution_MatlabExample/input_nrrd/dwi_b0_hr.nii";
  ReaderType::Pointer hriReader = ReaderType::New();
  hriReader->SetFileName(hriFileName);
  hriReader->Update();
  FloatImageType::Pointer hri = hriReader->GetOutput();

  const bool hri_ActualXDimensionIsOdd = (hri->GetLargestPossibleRegion().GetSize()[0] % 2) ? true : false;
  FloatImageType::Pointer test = GetInverseFFT(GetForwardFFT(hri,1.0),hri_ActualXDimensionIsOdd,1.0);
  WriteFile(test.GetPointer(),"/tmp/SR/toFromFFT.nii");

  itk::ShrinkImageFilter<FloatImageType,FloatImageType>::Pointer intensityReader = itk::ShrinkImageFilter<FloatImageType,FloatImageType>::New();
  intensityReader->SetInput(hri);
  intensityReader->SetShrinkFactors(2);
  intensityReader->Update();
  FloatImageType::Pointer lri = intensityReader->GetOutput();

  WriterType::Pointer hriWriter = WriterType::New();
  hriWriter->SetInput(hri);
  hriWriter->SetFileName("/tmp/SR/original.nii.gz");
  hriWriter->Update();
  FloatImageType::Pointer downsampledImage = IdentityResampleByFFT(hri.GetPointer(),lri.GetPointer());
  hriWriter->SetInput(downsampledImage);
  hriWriter->SetFileName("/tmp/SR/downsampled.nii.gz");
  hriWriter->Update();
  FloatImageType::Pointer upsampledImage = IdentityResampleByFFT(downsampledImage.GetPointer(),hri.GetPointer());
  hriWriter->SetInput(upsampledImage);
  hriWriter->SetFileName("/tmp/SR/upsampled.nii.gz");
  hriWriter->Update();

#if 0
  const std::string lriFileName = "/Shared/johnsonhj/HDNI/20160709_SuperResolution_MatlabExample/input_nrrd/dwi_b0_lr.nrrd";//argv[1];
  ReaderType::Pointer intensityReader = ReaderType::New();
  intensityReader->SetFileName(lriFileName);
  intensityReader->Update();
  FloatImageType::Pointer lri = intensityReader->GetOutput();

  HalfHermetianImageType::Pointer upsampledFreq = HalfHermetianImageType::New();
  {
    FloatImageType::RegionType newRegion = hri->GetLargestPossibleRegion();
    {
      FloatImageType::SizeType temp_size=hri->GetLargestPossibleRegion().GetSize();
      temp_size[0] = temp_size[0]/2 + 1; //Only store 1/2 of the first dimension to take advantage of symmetry
      newRegion.SetSize(temp_size);
    }
    upsampledFreq->CopyInformation(hri);
    upsampledFreq->SetRegions(newRegion);
    upsampledFreq->Allocate();
    upsampledFreq->FillBuffer(itk::NumericTraits<typename HalfHermetianImageType::PixelType>::ZeroValue());
  }


#if 1
  std::cout << "===============\n" << upsampledFreq << std::endl;

  TForwardFFT::Pointer TEMP = TForwardFFT::New();
  TEMP->SetInput(hri);
  TEMP->Update();
  upsampledFreq= TEMP->GetOutput();
  {
  itk::ImageFileWriter<FloatImageType>::Pointer cmplxWriter = itk::ImageFileWriter<FloatImageType>::New();
  itk::ComplexToModulusImageFilter<HalfHermetianImageType,FloatImageType>::Pointer tfmCmplx2Real = itk::ComplexToModulusImageFilter<HalfHermetianImageType,FloatImageType>::New();
  tfmCmplx2Real->SetInput(upsampledFreq);
  tfmCmplx2Real->Update();
  cmplxWriter->SetInput(tfmCmplx2Real->GetOutput());
  cmplxWriter->SetFileName("/tmp/SR/ffthri.nii.gz");
  cmplxWriter->Update();
  }
  upsampledFreq->FillBuffer(std::complex<float>(0,0));
  std::cout << "===============\n" << upsampledFreq << std::endl;
#endif
  TForwardFFT::Pointer fft = TForwardFFT::New();
  fft->SetInput( lri );
  const bool lri_ActualXDimensionIsOdd = ( lri->GetLargestPossibleRegion().GetSize()[0] % 2 ) ? true: false;
  fft->Update();
  HalfHermetianImageType::Pointer lri_cmplHH = fft->GetOutput();

  //Rescale FFT for upsampled size and downsampled versions.  N^{DIM} for downsampling and N^{DIM} for upsampling is N^{2*DIM}
  const float FFTScaler = std::pow(
    static_cast<float>(hri->GetLargestPossibleRegion().GetNumberOfPixels()) / static_cast<float>(lri->GetLargestPossibleRegion().GetNumberOfPixels()),
    6.0F);
  typedef itk::ImageRegionIteratorWithIndex<HalfHermetianImageType> cmplxHHIteratorType;
  cmplxHHIteratorType lriIter(lri_cmplHH,lri_cmplHH->GetLargestPossibleRegion());

  const HalfHermetianImageType::SizeType lriSize = lri->GetLargestPossibleRegion().GetSize();
  const HalfHermetianImageType::SizeType hriSize = hri->GetLargestPossibleRegion().GetSize();

  for(lriIter.GoToBegin(); ! lriIter.IsAtEnd(); ++lriIter)
  {
    HalfHermetianImageType::IndexType ioIndex = lriIter.GetIndex();
    const bool isInside = ConvertLRI2HRI( ioIndex , lriSize, hriSize);
    if(isInside)
    {
      upsampledFreq->SetPixel(ioIndex,FFTScaler*lriIter.Value());
    }
  }
  {
  itk::ImageFileWriter<FloatImageType>::Pointer cmplxWriter = itk::ImageFileWriter<FloatImageType>::New();
  itk::ComplexToModulusImageFilter<HalfHermetianImageType,FloatImageType>::Pointer tfmCmplx2Real = itk::ComplexToModulusImageFilter<HalfHermetianImageType,FloatImageType>::New();
  tfmCmplx2Real->SetInput(upsampledFreq);
  tfmCmplx2Real->Update();
  cmplxWriter->SetInput(tfmCmplx2Real->GetOutput());
  cmplxWriter->SetFileName("/tmp/SR/fftlripad.nii.gz");
  cmplxWriter->Update();
  }

  TInverseFFT::Pointer ifft = TInverseFFT::New();

  ifft->SetInput( upsampledFreq );
  const bool hri_ActualXDimensionIsOdd = ( hri->GetLargestPossibleRegion().GetSize()[0] % 2 ) ? true: false;
  //const bool hri_ActualXDimensionIsOdd = TEMP->GetActualXDimensionIsOdd();
  ifft->SetActualXDimensionIsOdd(hri_ActualXDimensionIsOdd);
  ifft->Update();
  FloatImageType::Pointer hriOut = ifft->GetOutput();
  hriOut->CopyInformation(hri);


  const std::string outFileName = "/tmp/SR/out.nrrd";
  WriterType::Pointer hriWriter = WriterType::New();
  hriWriter->SetInput(hriOut);
  hriWriter->SetFileName(outFileName);
  hriWriter->Update();

#endif
#endif
  return EXIT_SUCCESS;
}
