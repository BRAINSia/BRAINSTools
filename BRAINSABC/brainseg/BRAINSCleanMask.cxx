#include "CleanBrainLabelMap.h"
#include "itkIO.h"
#include "BRAINSCleanMaskCLP.h"
#include "BRAINSThreadControl.h"

int main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);
  if( inputVolume == "" )
    {
    std::cerr << "No input volume name given" << std::endl;
    return 1;
    }
  if( outputVolume == "" )
    {
    std::cerr << "No output volume name given" << std::endl;
    return 1;
    }

  typedef itk::Image<unsigned char, 3> ImageType;

  ImageType::Pointer input;
  try
    {
    input = itkUtil::ReadImage<ImageType>(inputVolume);
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "error reading " << inputVolume << std::endl
              <<  e << std::endl;
    return 1;
    }
  catch( ... )
    {
    std::cerr << "Unable to open " << inputVolume << std::endl;
    exit(1);
    }
  ImageType::Pointer output;
  try
    {
    output = CleanBrainLabelMap<ImageType, ImageType>(input);
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr <<  e << std::endl;
    return 1;
    }
  catch( ... )
    {
    std::cerr << "Error during processing of "
              << inputVolume << std::endl;
    return 1;
    }

  try
    {
    itkUtil::WriteImage<ImageType>(output, outputVolume);
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "error writing " << inputVolume << std::endl
              <<  e << std::endl;
    return 1;
    }
  catch( ... )
    {
    std::cerr << "Unable to write " << outputVolume << std::endl;
    return 1;
    }
  return 0;
}
