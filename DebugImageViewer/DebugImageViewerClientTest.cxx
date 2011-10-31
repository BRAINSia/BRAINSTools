#include <iostream>
#include "DebugImageViewerClient.h"
#include "sstream"
#include "itkIO.h"

int main(int argc, char * *argv)
{
  typedef itk::Image<float, 3> ImageType;

  if( argc != 3 )
    {
    exit(1);
    }

  std::stringstream s(argv[1]);
  int               viewIndex;
  s >> viewIndex;
  ImageType::Pointer img = itkUtil::ReadImage<ImageType>(argv[2]);
  if( img.IsNull() )
    {
    std::cerr << "Can't open " << argv[1] << std::endl;
    exit(1);
    }
  DebugImageViewerClient disp;
  disp.SetEnabled(true);
  disp.SendImage<ImageType>(img, viewIndex);
  exit(0);
}
