
#include "itkImage.h"
#include "itkVector.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkIterativeInverseDeformationFieldImageFilter.h"
#include "itkMultiThreadIterativeInverseDeformationFieldImageFilter.h"

int itkIterativeInverseDeformationFieldFilterTest(int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "You must supply two output filenames" << std::endl;
    return EXIT_FAILURE;
    }

  typedef  float PixelType;
  const unsigned int dims = 3;
  typedef itk::Image<itk::Vector<PixelType, dims>, dims> ImageType;
  typedef itk::ImageFileWriter<ImageType>                WriterType;

  const ImageType::SizeType  imageSize = {{64, 64, 64}};
  const ImageType::IndexType imageIndex = {{0, 0, 0}};
  ImageType::RegionType      region;
  region.SetSize(imageSize);
  region.SetIndex(imageIndex);
  ImageType::Pointer img = ImageType::New();
  img->SetRegions(region);
  img->Allocate();

  ImageType::PixelType zeros;
  zeros[0] = 15.0;
  zeros[1] = -109.0;
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

  typedef itk::IterativeInverseDeformationFieldImageFilter<ImageType, ImageType> InverseDeformationFieldImageType;
  InverseDeformationFieldImageType::Pointer inverse1 = InverseDeformationFieldImageType::New();
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

  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetInput( inverse1->GetOutput() );
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

  typedef itk::MultiThreadIterativeInverseDeformationFieldImageFilter<ImageType,
                                                                      ImageType> InverseDeformationField1ImageType;
  InverseDeformationField1ImageType::Pointer inverse2 = InverseDeformationField1ImageType::New();
  inverse2->SetInput(img);
  inverse2->SetStopValue(1.0e-7);
  inverse2->Update();
/*
  itk::ImageRegionIterator<ImageType> it1(inverse1->GetOutput(), inverse1->GetOutput()->GetRequestedRegion());
  it1.GoToBegin();

  itk::ImageRegionIterator<ImageType> it2(inverse2->GetOutput(), inverse2->GetOutput()->GetRequestedRegion());
  it2.GoToBegin();
  while( ! it2.IsAtEnd() )
  {
    std::cout<<"1:"<<it1.Value()<<std::endl;
    std::cout<<"2:"<<it2.Value()<<std::endl;
    ++it2;
    ++it1;
  }
*/
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
  return 0;
}
