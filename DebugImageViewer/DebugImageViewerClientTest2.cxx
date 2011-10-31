#include <iostream>
#include "DebugImageViewerClient.h"
#include "sstream"
#include "itkIO.h"

int main(int, char * *)
{
  typedef itk::Image<unsigned char, 3> ImageType;

  int viewIndex(0);

  ImageType::RegionType region;
  region.SetSize(0, 16);
  region.SetSize(1, 16);
  region.SetSize(2, 1);
  region.SetIndex(0, 0);
  region.SetIndex(1, 0);
  region.SetIndex(2, 0);

  ImageType::SpacingType spacing;
  spacing[0] = 1.0;
  spacing[1] = 1.0;
  spacing[2] = 1.0;

  ImageType::Pointer img =
    itkUtil::AllocateImageFromRegionAndSpacing<ImageType>
      (region, spacing);
  ImageType::IndexType index;
  index[2] = 0;
  for( unsigned i = 0; i < 16; i++ )
    {
    index[1] = i;
    for( unsigned j = 0; j < 16; j++ )
      {
      index[0] = j;
      if( j == 0 || j == 15 )
        {
        img->SetPixel(index, 255);
        }
      else
        {
        img->SetPixel(index, i);
        }
      }
    }
  DebugImageViewerClient disp;
  disp.SetEnabled(true);
  disp.SendImage<ImageType>(img, viewIndex);
  exit(0);
}
