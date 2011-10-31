#include <iostream>
#include "DebugImageViewerClient.h"
#include "sstream"
#include <itkVector.h>
#include <itkIO.h>

int main(int argc, char * *argv)
{
  typedef itk::Image<itk::Vector<float, 3>, 3> ImageType;

  if( argc != 2 )
    {
    exit(1);
    }

  ImageType::Pointer img = itkUtil::ReadImage<ImageType>(argv[1]);
  if( img.IsNull() )
    {
    std::cerr << "Can't open " << argv[1] << std::endl;
    exit(1);
    }
  DebugImageViewerClient disp;
  disp.SendImage<ImageType>(img, 0, 0);
  disp.SendImage<ImageType>(img, 1, 1);
  disp.SendImage<ImageType>(img, 2, 2);
  exit(0);
}
