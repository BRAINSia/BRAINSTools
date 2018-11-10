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
#include "itkVector.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
// #include "itkIterativeInverseDisplacementFieldImageFilter.h"
// #include "itkIterativeInverseDisplacementFieldImageFilter1.h"
#include "itkICCIterativeInverseDisplacementFieldImageFilter.h"

int itkIterativeInverseDisplacementFieldFilterTest2(int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "You must supply two output filenames" << std::endl;
    return EXIT_FAILURE;
    }

  using PixelType = float;
  constexpr unsigned int dims = 3;
  using ImageType = itk::Image<itk::Vector<PixelType, dims>, dims>;
  using WriterType = itk::ImageFileWriter<ImageType>;

  const ImageType::SizeType  imageSize = {{16, 16, 16}};
  const ImageType::IndexType imageIndex = {{0, 0, 0}};
  ImageType::RegionType      region;
  region.SetSize(imageSize);
  region.SetIndex(imageIndex);
  ImageType::Pointer img = ImageType::New();
  img->SetRegions(region);
  img->Allocate();

  ImageType::PixelType zeros;
  zeros[0] = 1.0;
  zeros[1] = 2.0;
  zeros[2] = 3.0;

  itk::ImageRegionIterator<ImageType> it(img, img->GetRequestedRegion() );
  it.GoToBegin();
  while( !it.IsAtEnd() )
    {
    zeros[0] += 0.1;
    zeros[1] += 0.2;
    zeros[2] += -0.1;
    it.Value() =  zeros;
    ++it;
    }

  using InverseDisplacementFieldImageType = itk::ICCIterativeInverseDisplacementFieldImageFilter<ImageType, ImageType>;
  InverseDisplacementFieldImageType::Pointer inverse1 = InverseDisplacementFieldImageType::New();
  inverse1->SetInput(img);
  inverse1->SetStopValue(1.0e-7);
  try
    {
    inverse1->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while inverse image" << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
    }

  InverseDisplacementFieldImageType::Pointer inverse2 = InverseDisplacementFieldImageType::New();
  inverse2->SetInput(inverse1->GetOutput() );
  inverse2->SetStopValue(1.0e-7);
  try
    {
    inverse2->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while inverse image" << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
    }

  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetInput( img );
  writer1->SetFileName(argv[1]);
  try
    {
    writer1->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while writing image" << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
    }

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( inverse2->GetOutput() );
  writer2->SetFileName(argv[2]);
  try
    {
    writer2->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while writing image" << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
    }

  itk::ImageRegionIterator<ImageType> it2(img, img->GetRequestedRegion() );
  itk::ImageRegionIterator<ImageType> it1(inverse2->GetOutput(), inverse2->GetOutput()->GetRequestedRegion() );
  it1.GoToBegin();
  it2.GoToBegin();

  while( !it1.IsAtEnd() )
    {
    std::cout << "1:" << it1.Value() << std::endl;
    ++it1;
    std::cout << "2:" << it2.Value() << std::endl;
    ++it2;
    }

  return 0;
}
