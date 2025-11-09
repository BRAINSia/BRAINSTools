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
 * Author: Ali Ghayoor
 * at SINAPSE Lab,
 * The University of Iowa 2015
 */

#include "GeneratePurePlugMask.h"
#include "GeneratePurePlugMaskCLP.h"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  using FloatImagePointer = FloatImageType::Pointer;
  using InputImageList = std::vector<FloatImagePointer>;
  using LocalReaderType = itk::ImageFileReader<FloatImageType>;
  using MaskImageType = itk::Image<unsigned char, 3>;

  std::vector<std::string> inputFileNames;
  if (inputImageModalities.size() > 1)
  {
    inputFileNames = inputImageModalities;
  }
  else
  {
    std::cerr << "ERROR: At least two image modalities are needed to generate pure plug mask." << std::endl;
    return EXIT_FAILURE;
  }
  const unsigned int numberOfImageModalities = inputFileNames.size(); // number of modality images

  // Read the input modalities and set them in a vector of images
  using LocalReaderPointer = LocalReaderType::Pointer;

  InputImageList inputImageModalitiesList;
  for (unsigned int i = 0; i < numberOfImageModalities; i++)
  {
    std::cout << "Reading image: " << inputFileNames[i] << std::endl;

    LocalReaderPointer imgreader = LocalReaderType::New();
    imgreader->SetFileName(inputFileNames[i].c_str());

    try
    {
      imgreader->Update();
    }
    catch (...)
    {
      std::cout << "ERROR:  Could not read image " << inputFileNames[i] << "." << std::endl;
      return EXIT_FAILURE;
    }

    inputImageModalitiesList.emplace_back(imgreader->GetOutput());
  }

  // define step size based on the number of sub-samples at each direction
  MaskImageType::SizeType numberOfContinuousIndexSubSamples;
  numberOfContinuousIndexSubSamples[0] = numberOfSubSamples[0];
  numberOfContinuousIndexSubSamples[1] = numberOfSubSamples[1];
  numberOfContinuousIndexSubSamples[2] = numberOfSubSamples[2];

  MaskImageType::Pointer mask = GeneratePurePlugMask<FloatImageType, MaskImageType>(
    inputImageModalitiesList, threshold, numberOfContinuousIndexSubSamples, false, verbose);

  using MaskWriterType = itk::ImageFileWriter<MaskImageType>;
  MaskWriterType::Pointer writer = MaskWriterType::New();
  writer->SetInput(mask);
  writer->SetFileName(outputMaskFile);
  writer->UseCompressionOn();
  writer->Update();

  return EXIT_SUCCESS;
}
