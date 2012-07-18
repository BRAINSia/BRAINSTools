

#include "itkImage.h"
#include "itkVectorFFTWComplexConjugateToRealImageFilter.h"
#include "itkVectorFFTWRealToComplexConjugateImageFilter.h"
#include "itkVector.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"

int itkVectorFFTWTest(int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "You must supply two output filenames" << std::endl;
    return EXIT_FAILURE;
    }

  typedef  float PixelType;
  const unsigned int dims = 3;
  typedef itk::Image<itk::Vector<PixelType, dims>, dims>                            ImageType;
  typedef itk::VectorFFTWComplexConjugateToRealImageFilter<ImageType::PixelType, 3> FFTWComplexToRealImageType;
  typedef itk::VectorFFTWRealToComplexConjugateImageFilter<ImageType::PixelType, 3> FFTWRealToComplexImageType;
  typedef itk::ImageFileWriter<ImageType>                                           WriterType;

  const ImageType::SizeType  imageSize = {{32, 32, 32}};
  const ImageType::IndexType imageIndex = {{0, 0, 0}};
  ImageType::RegionType      region;
  region.SetSize(imageSize);
  region.SetIndex(imageIndex);
  ImageType::Pointer img = ImageType::New();
  img->SetRegions(region);
  img->Allocate();

  ImageType::PixelType zeros;
  for( unsigned int j = 0; j < 3; j++ )
    {
    zeros[j] = 1.0;
    }

  itk::ImageRegionIterator<ImageType> it(img, img->GetRequestedRegion() );

  while( !it.IsAtEnd() )
    {
    zeros[0] += 1.0;
    zeros[1] += 0.5;
    zeros[2] += -0.1;
    it.Value() =  zeros;
    ++it;
    }

  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetInput( img );
  writer1->SetFileName(argv[1]);
  try
    {
    writer1->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while writing image" << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
    }

  FFTWRealToComplexImageType::Pointer fft = FFTWRealToComplexImageType::New();
  fft->SetInput(img);
  fft->Update();

  FFTWComplexToRealImageType::Pointer invFFT = FFTWComplexToRealImageType::New();
  invFFT->SetInput(fft->GetOutput() );
  invFFT->Update();

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( invFFT->GetOutput() );
  writer2->SetFileName(argv[2]);
  try
    {
    writer2->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while writing image" << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
    }

  return 0;
}
