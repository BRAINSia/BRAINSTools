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
/// \file BlendImageFilterGTest.cxx
/// \brief GTest port of BlendImageFilterTest.cxx (issue #403 — GTest migration).
///
/// Verifies that itkBlendImageFilter applies blend weights correctly by
/// comparing each output pixel against the analytic formula:
///   out[i] = in1[i] * blend1 + in2[i] * blend2
///
/// This replaces the original test that was commented out behind
/// BRAINSTools_MAX_TEST_LEVEL > 8.  With GTest the test is always enabled
/// and each pixel failure is reported individually.

#include <gtest/gtest.h>

#include "itkImage.h"
#include "itkBlendImageFilter.h"
#include "itkRandomImageSource.h"
#include "itkImageRegionConstIterator.h"

TEST(BlendImageFilter, BlendWeightsProduceCorrectPixelValues)
{
  using ImageType = itk::Image<float, 2>;

  // Generate two independent random images
  ImageType::Pointer images[2];
  for (unsigned int i = 0; i < 2; ++i)
  {
    using SourceType = itk::RandomImageSource<ImageType>;
    const auto               source = SourceType::New();
    ImageType::SizeValueType size[2] = { 64, 64 };
    source->SetSize(size);
    source->SetMin(0.0f);
    source->SetMax(1.0f);
    source->Update();
    images[i] = source->GetOutput();
  }

  // Run the blend filter with blend1=0.2, blend2=0.8
  using BlendImageFilterType = itk::BlendImageFilter<ImageType, ImageType>;
  const auto filter = BlendImageFilterType::New();
  filter->SetInput1(images[0]);
  filter->SetInput2(images[1]);
  filter->SetBlend1(0.2);
  filter->SetBlend2(0.8);
  filter->Update();

  const ImageType::Pointer blendImage = filter->GetOutput();

  itk::ImageRegionConstIterator<ImageType> it1(images[0], images[0]->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<ImageType> it2(images[1], images[1]->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<ImageType> itBlend(blendImage, blendImage->GetLargestPossibleRegion());

  for (; !it1.IsAtEnd() && !it2.IsAtEnd() && !itBlend.IsAtEnd(); ++it1, ++it2, ++itBlend)
  {
    const float expected = (it1.Get() * 0.2f) + (it2.Get() * 0.8f);
    EXPECT_NEAR(itBlend.Get(), expected, 1e-4f);
  }
}
