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
#include "itkImage.h"
#include "itkBlendImageFilter.h"
#include "itkRandomImageSource.h"
#include "itkImageRegionConstIterator.h"
#include "itkMath.h"

int
main(int, char **)
{
  using ImageType = itk::Image<float, 2>;

  std::cout << "Create input image using RandomImageSource" << std::endl;
  ImageType::Pointer images[2];
  for (unsigned i = 0; i < 2; i++)
  {
    using SourceType = itk::RandomImageSource<ImageType>;
    SourceType::Pointer      source = SourceType::New();
    ImageType::SizeValueType size[2] = { 64, 64 };
    source->SetSize(size);
    source->SetMin(0.0);
    source->SetMax(1.0);
    source->Update();
    images[i] = source->GetOutput();
  }
  using BlendImageFilterType = itk::BlendImageFilter<ImageType, ImageType>;
  BlendImageFilterType::Pointer filter = BlendImageFilterType::New();
  filter->SetInput1(images[0]);
  filter->SetInput2(images[1]);
  filter->SetBlend1(0.2);
  filter->SetBlend2(0.8);
  filter->Update();
  ImageType::Pointer                       blendImage = filter->GetOutput();
  itk::ImageRegionConstIterator<ImageType> it1(images[0], images[0]->GetLargestPossibleRegion()),
    it2(images[1], images[1]->GetLargestPossibleRegion()), itBlend(blendImage, blendImage->GetLargestPossibleRegion());
  for (; !it1.IsAtEnd() && !it2.IsAtEnd() && !itBlend.IsAtEnd(); ++it1, ++it2, ++itBlend)
  {
    float blend = (it1.Get() * 0.2) + (it2.Get() * 0.8);
    if (std::fabs(blend - itBlend.Get()) > 0.0001)
    {
      std::cerr << "Expected " << blend << " found " << itBlend.Get() << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
