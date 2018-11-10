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

int main(int argc, char * *argv)
{
  using ImageType = itk::Image<float, 3>;

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
