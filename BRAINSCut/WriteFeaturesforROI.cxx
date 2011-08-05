#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
// Iterator
#include <itkImageRegionConstIterator.h>

int
main(int argc, char * *argv)
{
  // PARSE_ARGS;
  if( argc < 2 )
    {
    std::cout << " Usage:: Need classified Image to calculate " << std::endl
              << argv[0] << "   "
              << " [input(feature) image] [mask image]"
              << std::endl;
    return 0;
    }
  std::cout << "Input File name: " << argv[1]
            << "Input Mask anem : " << argv[2] << std::endl;

  typedef double PixelType;
  const unsigned int Dim = 3;
  typedef itk::Image<PixelType, Dim>      ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  typedef unsigned char                  MaskPixelType;
  typedef itk::Image<MaskPixelType, Dim> MaskType;
  typedef itk::ImageFileReader<MaskType> MaskReaderType;
  // Read Input Image
  ImageReaderType::Pointer InputImage = ImageReaderType::New();
  InputImage->SetFileName(argv[1]);
  InputImage->Update();

  // Read mask
  MaskReaderType::Pointer InputMask = MaskReaderType::New();
  InputMask->SetFileName(argv[2]);
  InputMask->Update();

  // Dilate mask by 3 pixel
  typedef itk::BinaryBallStructuringElement<MaskPixelType, Dim>
    SturucturingElementType;
  typedef itk::BinaryDilateImageFilter<MaskType,
                                       MaskType,
                                       SturucturingElementType> DilateFilterType;
  DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
  SturucturingElementType   structuringElement;
  structuringElement.SetRadius(2);
  structuringElement.CreateStructuringElement();

  binaryDilate->SetKernel(structuringElement);
  binaryDilate->SetInput( InputMask->GetOutput() );
  binaryDilate->SetForegroundValue(1);
  binaryDilate->SetBackgroundValue(0);
  unsigned char fgvalue = binaryDilate->GetForegroundValue();
  std::cout << "Foreground Value:: " << binaryDilate->GetForegroundValue();
  std::cout << "Background Value:: " << binaryDilate->GetBackgroundValue() << std::endl;
  binaryDilate->Update();
  // Write Dilated Image
  typedef itk::ImageFileWriter<MaskType> writerType;
  writerType::Pointer writer = writerType::New();
  writer->UseCompressionOn();
  writer->SetInput( binaryDilate->GetOutput() );
  std::string filename = argv[2];
  filename += "_dilated.nii.gz";
  writer->SetFileName(filename);
  writer->Update();

  // Iterator
  typedef itk::ImageRegionConstIterator<MaskType> ConstIteratorType;
  ConstIteratorType iter( InputMask->GetOutput(),
                          InputMask->GetOutput()->GetLargestPossibleRegion() );
  for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    ImageType::IndexType p_index = iter.GetIndex();
    if( binaryDilate->GetOutput()->GetPixel(p_index) == fgvalue
        && InputMask->GetOutput()->GetPixel(p_index) == fgvalue )
      {
      std::cout << "In  " << InputImage->GetOutput()->GetPixel(p_index)
                << std::endl;
      }
    if( binaryDilate->GetOutput()->GetPixel(p_index) == fgvalue
        && InputMask->GetOutput()->GetPixel(p_index) != fgvalue )
      {
      std::cout << "Out " << InputImage->GetOutput()->GetPixel(p_index)
                << std::endl;
      }
    }

  return 0;
}
