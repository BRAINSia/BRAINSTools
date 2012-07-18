

#include "itkLargestForegroundFilledMaskImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "CreateMaskCLP.h"
#include "itkIO.h"

int CreateMask(int argc, char *argv[])
{
  PARSE_ARGS;

  const bool debug = true;

  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image:     " <<  inputVolume << std::endl;
    std::cout << "Threshold:           " <<  threshold << std::endl;
    std::cout << "Closeing Size:   " <<  closingSize << std::endl;
    std::cout << "Output Volume:       " <<  outputVolume << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  if( inputVolume.size() == 0 )
    {
    std::cout << "Input Volume is misses!" << std::endl;
    exit(-1);
    }

  if( outputVolume.size() == 0 )
    {
    std::cout << "Output Volume is missed!" << std::endl;
    exit(-1);
    }

  const unsigned int Dimension = 3;
  typedef float                            PixelType;
  typedef itk::Image<PixelType, Dimension> ImageType;

  ImageType::Pointer inputImage
    = itkUtil::ReadImage<ImageType>( inputVolume );

  typedef itk::LargestForegroundFilledMaskImageFilter<ImageType> MaskFilterType;
  MaskFilterType::Pointer LFF = MaskFilterType::New();
  LFF->SetInput(inputImage);
  LFF->SetOtsuPercentileThreshold(threshold);
  LFF->SetClosingSize(closingSize);
  LFF->Update();

  typedef itk::ImageFileWriter<ImageType> ImageWriteType;
  ImageWriteType::Pointer writer = ImageWriteType::New();
  writer->SetInput(LFF->GetOutput() );
  writer->SetFileName(outputVolume);
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  // HACK:  BRAINS2 Masks are currently broken
  // The direction cosines are and the direction labels are not consistently being set.
  // itk::Brains2MaskImageIOFactory::RegisterOneFactory();

  return CreateMask(argc, argv);
}
