#include "CleanBrainLabelMap.h"
#include "itkIO.h"

int main(int argc, char * *argv)
{
  if( argc < 3 )
    {
    std::cerr << "Usage: BRAINSCleanMask inputLabelMap outputLabelMap" << std::endl;
    return 1;
    }
  typedef itk::Image<unsigned char, 3> ImageType;

  std::string inputName(argv[1]), outputName(argv[2]);

  ImageType::Pointer input = itkUtil::ReadImage<ImageType>(inputName);
  ImageType::Pointer output = CleanBrainLabelMap<ImageType, ImageType>(input);
  itkUtil::WriteImage<ImageType>(output, outputName);
}
