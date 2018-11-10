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
 * Author : Ali Ghayoor
 */

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"

#include <iomanip>

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage : " << argv[0] << " <4DImage> <3DComponentsPrefix> " << std::endl;
    return EXIT_FAILURE;
    }

  const std::string input4Dimage(  argv[1]  );
  const std::string outputPrefix( argv[2] );

  using PixelValueType = float;
  using Volume4DType = itk::Image<PixelValueType, 4>;
  using Volume3DType = itk::Image<PixelValueType, 3>;

  std::cout << "- Read image: " << input4Dimage << std::endl;
  using Image4DReaderType = itk::ImageFileReader<Volume4DType>;
  Image4DReaderType::Pointer image4DReader = Image4DReaderType::New();
  image4DReader->SetFileName( input4Dimage );
  try
    {
    image4DReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    }
  Volume4DType::Pointer inputVol = image4DReader->GetOutput();

  // "inputVol" is read as a 4D image. Here we extract each 3D component and write that to disk
  //
  Volume4DType::SizeType inputSize =
    inputVol->GetLargestPossibleRegion().GetSize();

  Volume4DType::IndexType inputIndex =
    inputVol->GetLargestPossibleRegion().GetIndex();

  const unsigned int volumeCount = inputSize[3];

  using ExtractFilterType = itk::ExtractImageFilter< Volume4DType, Volume3DType >;

  for( size_t componentNumber = 0; componentNumber < volumeCount; ++componentNumber )
    {
    Volume4DType::SizeType extractSize = inputSize;
    extractSize[3] = 0;
    Volume4DType::IndexType extractIndex = inputIndex;
    extractIndex[3] = componentNumber;
    Volume4DType::RegionType extractRegion(extractIndex, extractSize);

    ExtractFilterType::Pointer extracter = ExtractFilterType::New();
    extracter->SetExtractionRegion( extractRegion );
    extracter->SetInput( inputVol );
    extracter->SetDirectionCollapseToIdentity();
    extracter->Update();

    // Need to zeropad to ensure that files are printed in order.
    std::stringstream fNumber("");
    fNumber << std::setw(3) << std::setfill('0') << componentNumber;
    const std::string fn = outputPrefix + fNumber.str() + ".nii.gz";

    std::cout << "- Write image: " << fn << std::endl;

    using Image3DWriterType = itk::ImageFileWriter<Volume3DType>;
    Image3DWriterType::Pointer image3DWriter = Image3DWriterType::New();
    image3DWriter->SetFileName( fn );
    image3DWriter->SetInput( extracter->GetOutput() );
    try
      {
      image3DWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown while writing the image" << std::endl;
      std::cerr << excp << std::endl;
      }
    }
  return EXIT_SUCCESS;
}
