#ifndef FFTWUpsample_h
#define FFTWUpsample_h



#include "SRTypes.h"
#ifdef DEBUG
extern void testIndexConversion(); //A utility function to test
#endif

#include "MathUtils.h"


extern HalfHermetianImageType::Pointer GetLowPassFilterFFT( FloatImageType::Pointer inputImage,   itk::ImageBase<3>::Pointer referenceImageBase );

extern HalfHermetianImageType::Pointer A_fhp(FloatImageType::Pointer inHRRealImage, itk::ImageBase<3>::Pointer desiredOutputRes );
extern FloatImageType::Pointer At_fhp(HalfHermetianImageType::Pointer inLRCoeffs,
                               const bool inputFirstDimIsOdd,
                               itk::ImageBase<3>::Pointer desiredOutputRef);
extern FloatImageType::Pointer IdentityResampleByFFT( FloatImageType::Pointer inputImage,   itk::ImageBase<3>::Pointer referenceImageBase );

extern HalfHermetianImageType::Pointer CreateZeroFFTCoefficients(itk::ImageBase<3>::Pointer referenceImageBase);
extern void FFTWInit(const std::string path_for_wisdom);
extern HalfHermetianImageType::Pointer GetForwardFFT(FloatImageType::Pointer inputImage, PrecisionType FFTScaler=1.0);
extern FloatImageType::Pointer GetInverseFFT(HalfHermetianImageType::Pointer inputImage, const bool
referenceImageBase_ActualXDimensionIsOdd, PrecisionType FFTScaler=1.0);

extern HalfHermetianImageType::Pointer ReshapeFFT(
    const itk::ImageBase<3>::Pointer outputRealImageBase,
    HalfHermetianImageType::Pointer inputFreqCoeffs,
    const bool inFirstSpatialDeminsionIsOdd);

extern CVImageType::Pointer GetGradient(FloatImageType::Pointer inputImage);
extern DivergenceType::OutputImageType::Pointer GetDivergence(CVImageType::Pointer inputImage);

#include <itkImageDuplicator.h>
template<typename ImageType>
typename ImageType::Pointer DeepImageCopy(typename ImageType::Pointer in)
{
  typedef itk::ImageDuplicator<ImageType> ImDupType;
  typename ImDupType::Pointer imDup = ImDupType::New();
  imDup->SetInputImage(in);
  imDup->Update();
  return imDup->GetOutput();
}


#include <itkImageFileWriter.h>
template<typename ImageType>
void WriteFile(ImageType * img, const std::string outfilename)
{
  typedef typename itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outfilename);
  writer->SetInput(img);
  writer->Update();
}

#include <itkComplexToRealImageFilter.h>
#include <itkComplexToImaginaryImageFilter.h>
#include <itkComplexToModulusImageFilter.h>

extern void WriteComplexImages(HalfHermetianImageType::Pointer cmplxIn, const std::string prefix);
extern HalfHermetianImageType::Pointer  ReadComplexImages(const std::string prefix);

extern void MoveFFTCoeffs(HalfHermetianImageType::Pointer outputFreqCoeffs,
                   const bool outFirstSpatialDeminsionIsOdd,
                   HalfHermetianImageType::Pointer inputFreqCoeffs,
                   const bool inFirstSpatialDeminsionIsOdd );

extern FloatImageType::Pointer NormalizeDataComponent(FloatImageType::Pointer arr);

#endif //FFTWUpsample_h
