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
/*
 * Author : Eun Young (Regina) Kim
 */

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include <map>

#include "GenerateLabelMapFromProbabilityMapCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>

/*
 * main function
 */

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if (inputVolumes.empty())
  {
    std::cerr << argv[0] << ": Missing required probability maps" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int numberOfProbabilityMaps = inputVolumes.size();

  // Define Image Type
  using ProbabilityMapPixelType = float;
  constexpr unsigned int Dimension = 3;
  using ProbabilityMapImageType = itk::Image<ProbabilityMapPixelType, Dimension>;

  // Label Index should start from zero and increasing order.

  // read in images

  std::vector<ProbabilityMapImageType::Pointer> probabilityImages(numberOfProbabilityMaps);
  for (unsigned int indexInputImages = 0; indexInputImages < numberOfProbabilityMaps; indexInputImages++)
  {
    std::cout << "- Read image::" << inputVolumes[indexInputImages] << std::endl;
    using ProbabilityImageReaderType = itk::ImageFileReader<ProbabilityMapImageType>;

    ProbabilityImageReaderType::Pointer probabilityReader = ProbabilityImageReaderType::New();

    probabilityReader->SetFileName(inputVolumes[indexInputImages]);
    probabilityReader->Update();

    probabilityImages[indexInputImages] = probabilityReader->GetOutput();
  }

  // crerate empty label maps
  std::cout << "Create Label Map" << std::endl;
  using LabelMapPixelType = unsigned int;
  using LabelMapImageType = itk::Image<LabelMapPixelType, Dimension>;

  LabelMapImageType::Pointer labelImage = LabelMapImageType::New();

  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  std::cerr << " ProbImage:  " << probabilityImages[0] << std::endl;
  labelImage->CopyInformation(probabilityImages[0]);
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  labelImage->SetRegions(probabilityImages[0]->GetLargestPossibleRegion());
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  labelImage->Allocate();
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  labelImage->FillBuffer(0);
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;

  std::cout << "start iteration " << std::endl;
  // Iterate Images
  using IteratorType = itk::ImageRegionIterator<ProbabilityMapImageType>;

  IteratorType probabilityIterator(probabilityImages[0], probabilityImages[0]->GetLargestPossibleRegion());

  using LabelIteratorType = itk::ImageRegionIterator<LabelMapImageType>;
  LabelIteratorType labelIterator(labelImage, labelImage->GetLargestPossibleRegion());
  for (probabilityIterator.GoToBegin(), labelIterator.GoToBegin(); !probabilityIterator.IsAtEnd();
       ++probabilityIterator, ++labelIterator)
  {
    ProbabilityMapPixelType max = probabilityIterator.Get();
    LabelMapPixelType       label = 0;
    for (unsigned int indexInputImages = 1; indexInputImages < numberOfProbabilityMaps; indexInputImages++)
    {
      ProbabilityMapPixelType pixelValue =
        probabilityImages[indexInputImages]->GetPixel(probabilityIterator.GetIndex());
      if (pixelValue > max)
      {
        max = pixelValue;
        label = indexInputImages;
      }
    }

    labelIterator.Set(label);
  }

  // Image Writer
  using LabelWriterType = itk::ImageFileWriter<LabelMapImageType>;

  LabelWriterType::Pointer labelWriter = LabelWriterType::New();

  labelWriter->SetFileName(outputLabelVolume);
  labelWriter->SetInput(labelImage);

  labelWriter->Update();

  return 0;
}
