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


#include "itkLargestForegroundFilledMaskImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "CreateMaskCLP.h"
#include "itkIO.h"
#include "BRAINSCommonLib.h"

int CreateMask(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const bool debug = true;

  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image:     " <<  inputVolume << std::endl;
    std::cout << "Threshold:           " <<  threshold << std::endl;
    std::cout << "Closeing Size:   " <<  closingSize << std::endl;
    std::cout << "Output Volume:       " <<  outputVolume << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  if( inputVolume.size() == 0 )
    {
    std::cout << "Input Volume is misses!" << std::endl;
    exit(-1);
    }

  if( outputVolume.size() == 0 )
    {
    std::cout << "Output Volume is missed!" << std::endl;
    exit(-1);
    }

  constexpr unsigned int Dimension = 3;
  using PixelType = float;
  using ImageType = itk::Image<PixelType, Dimension>;

  ImageType::Pointer inputImage
    = itkUtil::ReadImage<ImageType>( inputVolume );

  using MaskFilterType = itk::LargestForegroundFilledMaskImageFilter<ImageType>;
  MaskFilterType::Pointer LFF = MaskFilterType::New();
  LFF->SetInput(inputImage);
  LFF->SetOtsuPercentileThreshold(threshold);
  LFF->SetClosingSize(closingSize);
  LFF->Update();

  using ImageWriteType = itk::ImageFileWriter<ImageType>;
  ImageWriteType::Pointer writer = ImageWriteType::New();
  writer->SetInput(LFF->GetOutput() );
  writer->SetFileName(outputVolume);
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  // HACK:  BRAINS2 Masks are currently broken
  // The direction cosines are and the direction labels are not consistently being set.
  // itk::Brains2MaskImageIOFactory::RegisterOneFactory();

  return CreateMask(argc, argv);
}
