#include "StandardizeMaskIntensity.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  typedef itk::Image<float, 3>         ImageType;
  typedef itk::Image<unsigned char, 3> MaskImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName(argv[1]);
  imageReader->Update();
  ImageType::Pointer image = imageReader->GetOutput();

  MaskImageType::Pointer mask = NULL; // itkUtil::ReadImage<MaskImageType>(
                                      // argv[2] );
  if( argc == 4 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName(argv[3]);
    maskReader->Update();
    mask = maskReader->GetOutput();
    }

  const double               lFract = 0.00005;
  const double               uFract = 1.0 - lFract;
  const ImageType::PixelType lTarget = 1;
  const ImageType::PixelType uTarget = 0.95 * MAX_IMAGE_OUTPUT_VALUE;
  const ImageType::PixelType clipMin = 0;
  const ImageType::PixelType clipMax = MAX_IMAGE_OUTPUT_VALUE;

  ImageType::Pointer result = StandardizeMaskIntensity<ImageType, MaskImageType>(image,
                                                                                 mask,
                                                                                 lFract,
                                                                                 uFract,
                                                                                 lTarget,
                                                                                 uTarget,
                                                                                 clipMin,
                                                                                 clipMax);

  if( result.IsNull() )
    {
    return 2;
    }

  typedef itk::ImageFileWriter<ImageType> FloatWriterType;
  FloatWriterType::Pointer writer = FloatWriterType::New();

  writer->SetInput(result);
  writer->SetFileName(argv[2]);
  writer->UseCompressionOn();
  writer->Update();
  return 0;
}
