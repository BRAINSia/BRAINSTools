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

/*
 * main function
 */

int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);
  if( inputVolumes.size() < 1 )
    {
    std::cerr << argv[0] << ": Missing required probability maps"
              << std::cerr;
    return EXIT_FAILURE;
    }

  const unsigned int numberOfProbabilityMaps = inputVolumes.size();

  // Define Image Type
  typedef float ProbabilityMapPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image<ProbabilityMapPixelType, Dimension> ProbabilityMapImageType;

  // Label Index should start from zero and increasing order.

  // read in images

  std::vector<ProbabilityMapImageType::Pointer> probabilityImages( numberOfProbabilityMaps );
  for( unsigned int indexInputImages = 0;
       indexInputImages < numberOfProbabilityMaps;
       indexInputImages++ )
    {
    std::cout << "- Read image::"
              << inputVolumes[indexInputImages]
              << std::endl;
    typedef itk::ImageFileReader<ProbabilityMapImageType> ProbabilityImageReaderType;

    ProbabilityImageReaderType::Pointer probabilityReader = ProbabilityImageReaderType::New();

    probabilityReader->SetFileName( inputVolumes[indexInputImages] );
    probabilityReader->Update();

    probabilityImages[indexInputImages] = probabilityReader->GetOutput();
    }

  // crerate empty label maps
  std::cout << "Create Label Map" << std::endl;
  typedef unsigned int                             LabelMapPixelType;
  typedef itk::Image<LabelMapPixelType, Dimension> LabelMapImageType;

  LabelMapImageType::Pointer labelImage = LabelMapImageType::New();

  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  std::cerr << " ProbImage:  " << probabilityImages[0] << std::endl;
  labelImage->CopyInformation( probabilityImages[0] );
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  labelImage->SetRegions( probabilityImages[0]->GetLargestPossibleRegion() );
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  labelImage->Allocate();
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;
  labelImage->FillBuffer(0);
  std::cerr << "here" << __FILE__ << " " << __LINE__ << std::endl;

  std::cout << "start iteration " << std::endl;
  // Iterate Images
  typedef itk::ImageRegionIterator<ProbabilityMapImageType> IteratorType;

  IteratorType probabilityIterator( probabilityImages[0],
                                    probabilityImages[0]->GetLargestPossibleRegion() );

  typedef itk::ImageRegionIterator<LabelMapImageType> LabelIteratorType;
  LabelIteratorType labelIterator( labelImage, labelImage->GetLargestPossibleRegion() );
  for( probabilityIterator.GoToBegin(), labelIterator.GoToBegin();
       !probabilityIterator.IsAtEnd();
       ++probabilityIterator, ++labelIterator )
    {
    ProbabilityMapPixelType max = probabilityIterator.Get();
    LabelMapPixelType       label = 0;
    for( unsigned int indexInputImages = 1; indexInputImages < numberOfProbabilityMaps;
         indexInputImages++ )
      {
      ProbabilityMapPixelType pixelValue =
        probabilityImages[indexInputImages]->GetPixel( probabilityIterator.GetIndex() );
      if(  pixelValue > max )
        {
        max = pixelValue;
        label = indexInputImages;
        }
      }

    labelIterator.Set( label );
    }

  // Image Writer
  typedef itk::ImageFileWriter<LabelMapImageType> LabelWriterType;

  LabelWriterType::Pointer labelWriter = LabelWriterType::New();

  labelWriter->SetFileName( outputLabelVolume );
  labelWriter->SetInput( labelImage );

  labelWriter->Update();

  return 0;
}
