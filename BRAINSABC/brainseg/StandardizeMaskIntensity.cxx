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
#include "StandardizeMaskIntensity.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  typedef itk::Image<float, 3>         ImageType;
  typedef itk::Image<unsigned char, 3> MaskImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName(argv[1]);
  imageReader->Update();
  ImageType::Pointer image = imageReader->GetOutput();

  MaskImageType::Pointer mask = NULL; // itkUtil::ReadImage<MaskImageType>(
                                      // argv[2] );
  if( argc == 4 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName(argv[3]);
    maskReader->Update();
    mask = maskReader->GetOutput();
    }

  constexpr double lFract = 0.00005;
  const double               uFract = 1.0 - lFract;
  constexpr ImageType::PixelType lTarget  = 1;
  const ImageType::PixelType uTarget = 0.95 * MAX_IMAGE_OUTPUT_VALUE;
  constexpr ImageType::PixelType clipMin  = 0;
  const ImageType::PixelType clipMax = MAX_IMAGE_OUTPUT_VALUE;

  ImageType::Pointer result = StandardizeMaskIntensity<ImageType, MaskImageType>(image,
                                                                                 mask,
                                                                                 lFract,
                                                                                 uFract,
                                                                                 lTarget,
                                                                                 uTarget,
                                                                                 clipMin,
                                                                                 clipMax);

  if( result.IsNull() )
    {
    return 2;
    }

  typedef itk::ImageFileWriter<ImageType> FloatWriterType;
  FloatWriterType::Pointer writer = FloatWriterType::New();

  writer->SetInput(result);
  writer->SetFileName(argv[2]);
  writer->UseCompressionOn();
  writer->Update();
  return 0;
}
