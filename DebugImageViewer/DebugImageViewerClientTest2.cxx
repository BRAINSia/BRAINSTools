/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include <iostream>
#include "DebugImageViewerClient.h"
#include "sstream"
#include "itkIO.h"

int
main(int, char **)
{
  using ImageType = itk::Image<unsigned char, 3>;

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

  ImageType::Pointer   img = itkUtil::AllocateImageFromRegionAndSpacing<ImageType>(region, spacing);
  ImageType::IndexType index;
  index[2] = 0;
  for (unsigned i = 0; i < 16; i++)
  {
    index[1] = i;
    for (unsigned j = 0; j < 16; j++)
    {
      index[0] = j;
      if (j == 0 || j == 15)
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
