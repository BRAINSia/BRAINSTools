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
/*  KENT  -- This CLP Wrapped test program needs to exercise
 *  the itkFindCenterOfBrainCLP.h class.
 *
 * As part of this testing, please move as many of the hard-coded
 * debuggging images out of the hxx files, and make it so that this
 * test program will create those images from the command line.
 *
 * You will have to make some more member variables of the class for
 * some of the intermediate images like "AfterGridComputationsForeground.nii.gz"
 * so that the class can expose them to the test program.
 *
 * Please also write an ADD_TEST section to the CMakeLists.txt file that will execute this test program.
 */
#include "itkFindCenterOfBrainFilter.h"
#include "itkIO.h"
#include <FindCenterOfBrainCLP.h>
#include "itkNumberToString.h"

int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  if (InputVolume.empty())
  {
    std::cerr << "FindCenterOfBrain: missing input image name" << std::endl;
    return 1;
  }
  using ImageType = itk::Image<unsigned char, 3>;
  using FindCenterFilterType = itk::FindCenterOfBrainFilter<ImageType>;
  using MaskImageType = FindCenterFilterType::MaskImageType;

  ImageType::Pointer inputImage = itkUtil::ReadImage<ImageType>(InputVolume);
  if (inputImage.IsNull())
  {
    std::cerr << "FindCenterOfBrain: Can't read input image " << InputVolume << std::endl;
    return 2;
  }

  FindCenterFilterType::Pointer filter = FindCenterFilterType::New();
  filter->SetInput(inputImage);

  MaskImageType::Pointer imageMask;
  if (!ImageMask.empty())
  {
    imageMask = itkUtil::ReadImage<MaskImageType>(ImageMask);
    if (imageMask.IsNull())
    {
      std::cerr << "FindCenterOfBrain: Can't read mask " << ImageMask << std::endl;
      return 3;
    }
    filter->SetImageMask(imageMask);
  }
  filter->SetMaximize(Maximize);
  filter->SetAxis(Axis);
  filter->SetOtsuPercentileThreshold(OtsuPercentileThreshold);
  filter->SetClosingSize(ClosingSize);
  filter->SetHeadSizeLimit(HeadSizeLimit);
  filter->SetHeadSizeEstimate(HeadSizeEstimate);
  filter->SetBackgroundValue(BackgroundValue);
  filter->SetGenerateDebugImages(GenerateDebugImages);
  try
  {
    filter->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 4;
  }
  FindCenterFilterType::PointType center = filter->GetCenterOfBrain();
  itk::NumberToString<double>     doubleConvert;
  std::cout << "Center Of Brain:"
            << " " << doubleConvert(center[0]) << " " << doubleConvert(center[1]) << " " << doubleConvert(center[2])
            << std::endl;
  if (!ClippedImageMask.empty())
  {
    MaskImageType::Pointer clippedMask = const_cast<MaskImageType *>(filter->GetClippedImageMask());
    itkUtil::WriteImage<MaskImageType>(clippedMask, ClippedImageMask);
  }
  if (!GenerateDebugImages)
  {
    return 0;
  }
  if (!DebugDistanceImage.empty())
  {
    FindCenterFilterType::DistanceImagePointer distImage = filter->GetDebugDistanceImage();
    itkUtil::WriteImage<FindCenterFilterType::DistanceImageType>(distImage, DebugDistanceImage);
  }
  if (!DebugGridImage.empty())
  {
    FindCenterFilterType::InputImagePointer gridImage = filter->GetDebugGridImage();
    itkUtil::WriteImage<ImageType>(gridImage, DebugGridImage);
  }
  if (!DebugAfterGridComputationsForegroundImage.empty())
  {
    MaskImageType::Pointer afterImage = filter->GetDebugAfterGridComputationsForegroundImage();
    itkUtil::WriteImage<MaskImageType>(afterImage, DebugAfterGridComputationsForegroundImage);
  }
  if (!DebugClippedImageMask.empty())
  {
    MaskImageType::Pointer clippedMask = filter->GetDebugClippedImageMask();
    itkUtil::WriteImage<MaskImageType>(clippedMask, DebugClippedImageMask);
  }
  if (!DebugTrimmedImage.empty())
  {
    ImageType::Pointer trimmedImage = filter->GetDebugTrimmedImage();
    itkUtil::WriteImage<ImageType>(trimmedImage, DebugTrimmedImage);
  }
  return EXIT_SUCCESS;
}
