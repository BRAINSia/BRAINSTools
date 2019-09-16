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
#if defined(_MSC_VER)
#  pragma warning(disable : 4786)
#endif

#ifdef __BORLANDC__
#  define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAddImageFilter.h"
#include "itkXorImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageDuplicator.h"

#include <BRAINSCommonLib.h>
#include "CleanUpOverlapLabelsCLP.h"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // check inputs
  if (inputBinaryVolumes.size() != outputBinaryVolumes.size())
  {
    std::cout << "* Error! input and output has to be matched! " << std::endl;
    return EXIT_FAILURE;
  }
  using BinaryPixelType = unsigned char;

  constexpr unsigned char Dim = 3;

  using InputBinaryImageType = itk::Image<BinaryPixelType, Dim>;

  using InputBinaryVolumeReaderType = itk::ImageFileReader<InputBinaryImageType>;

  std::vector<InputBinaryImageType::Pointer> inputBinaryVolumeVector;
  std::vector<std::string>::iterator         inputBinaryVolumeStringIt;
  for (inputBinaryVolumeStringIt = inputBinaryVolumes.begin(); inputBinaryVolumeStringIt < inputBinaryVolumes.end();
       ++inputBinaryVolumeStringIt)
  {
    std::cout << "* Read image: " << *inputBinaryVolumeStringIt << std::endl;
    InputBinaryVolumeReaderType::Pointer reader = InputBinaryVolumeReaderType::New();
    reader->SetFileName(*inputBinaryVolumeStringIt);
    reader->Update();

    inputBinaryVolumeVector.push_back(reader->GetOutput());
  }

  // Add all labels

  using ImageDuplicatorType = itk::ImageDuplicator<InputBinaryImageType>;

  ImageDuplicatorType::Pointer imageCopier = ImageDuplicatorType::New();
  imageCopier->SetInputImage(inputBinaryVolumeVector.front());
  imageCopier->Update();

  InputBinaryImageType::Pointer sumVolume = imageCopier->GetOutput();

  using AddImageFilterType = itk::AddImageFilter<InputBinaryImageType, InputBinaryImageType, InputBinaryImageType>;

  AddImageFilterType::Pointer adder = AddImageFilterType::New();
  for (unsigned int i = 1; // should start from second image
       i < inputBinaryVolumeVector.size();
       ++i)
  {
    adder->SetInput1(sumVolume);
    adder->SetInput2(inputBinaryVolumeVector[i]);
    adder->Update();

    sumVolume = adder->GetOutput();
  }

  // threshold summed image

  using ThresholdFilterType = itk::BinaryThresholdImageFilter<InputBinaryImageType, InputBinaryImageType>;

  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
  thresholder->SetInput(sumVolume);
  thresholder->SetOutsideValue(0);
  thresholder->SetInsideValue(1);
  thresholder->SetLowerThreshold(1);

  thresholder->Update();
  sumVolume = thresholder->GetOutput();

  // process one label at a time

  using AndImageFilterType = itk::AndImageFilter<InputBinaryImageType, InputBinaryImageType, InputBinaryImageType>;

  using XORImageFilterType = itk::XorImageFilter<InputBinaryImageType, InputBinaryImageType, InputBinaryImageType>;

  using WriterType = itk::ImageFileWriter<InputBinaryImageType>;
  for (unsigned int i = 0; i < inputBinaryVolumeVector.size(); ++i)
  {
    // write out overlap between current mask and summed mask
    ThresholdFilterType::Pointer maskThresholder = ThresholdFilterType::New();
    maskThresholder->SetInput(inputBinaryVolumeVector[i]);
    maskThresholder->SetOutsideValue(0);
    maskThresholder->SetInsideValue(1);
    maskThresholder->SetLowerThreshold(1);

    AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
    andFilter->SetInput(0, sumVolume);
    andFilter->SetInput(1, maskThresholder->GetOutput());
    andFilter->Update();

    std::cout << "Write binary volume : " << outputBinaryVolumes[i] << std::endl;

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(andFilter->GetOutput());
    writer->SetFileName(outputBinaryVolumes[i]);
    writer->Update();

    // subtract current mask from the summed mask
    XORImageFilterType::Pointer xorFilter = XORImageFilterType::New();
    xorFilter->SetInput1(sumVolume);
    xorFilter->SetInput2(maskThresholder->GetOutput());

    xorFilter->Update();

    sumVolume = xorFilter->GetOutput();
  }

  return EXIT_SUCCESS;
}
