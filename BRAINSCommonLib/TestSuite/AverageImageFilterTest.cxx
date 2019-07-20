
#include <itkImage.h>

#include <itkRandomImageSource.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkSubtractImageFilter.h>
#include <itkStatisticsImageFilter.h>

#include "itkAverageImageFilter.h"

int
main( int, char *[] )
{
  const unsigned int numTestImages( 4 );
  const unsigned int testImageDim( 4 );
  using FloatImage2DType = itk::Image< float, 2 >;

  using AverageImageFilterType = itk::AverageImageFilter< FloatImage2DType, FloatImage2DType >;
  using SubtractFilterType = itk::SubtractImageFilter< FloatImage2DType >;
  using StatFilterType = itk::StatisticsImageFilter< FloatImage2DType >;

  using FloatImageVector = std::vector< FloatImage2DType::Pointer >;


  using FloatImageConstIterator = itk::ImageRegionConstIterator< FloatImage2DType >;
  using FloatImageIterator = itk::ImageRegionIterator< FloatImage2DType >;
  using FloatImageConstIteratorVector = std::vector< FloatImageConstIterator >;

  FloatImage2DType::SizeType randomSize;
  randomSize[1] = randomSize[0] = testImageDim;

  FloatImageVector              inputImages;
  FloatImageConstIteratorVector inputIterators;

  for ( unsigned int i = 0; i < numTestImages; ++i )
  {
    itk::RandomImageSource< FloatImage2DType >::Pointer random;
    random = itk::RandomImageSource< FloatImage2DType >::New();
    random->SetMin( 0.0 );
    random->SetMax( 1000.0 );
    random->SetSize( randomSize );
    random->Update();
    inputImages.push_back( random->GetOutput() );
    FloatImageConstIterator curIt =
      FloatImageConstIterator( inputImages[i], inputImages[i]->GetLargestPossibleRegion() );
    inputIterators.push_back( curIt );
  }
  // do a manual average to compare with the filter output
  FloatImage2DType::Pointer testAvg = FloatImage2DType::New();
  testAvg->CopyInformation( inputImages[0] );
  testAvg->SetRegions( randomSize );
  testAvg->Allocate();
  FloatImageIterator avgIt( testAvg, testAvg->GetLargestPossibleRegion() );
  for ( avgIt.GoToBegin(); !avgIt.IsAtEnd(); ++avgIt )
  {
    double sum = 0;
    for ( unsigned i = 0; i < numTestImages; ++i )
    {
      sum += inputIterators[i].Value();
      ++inputIterators[i];
    }
    avgIt.Set( sum / numTestImages );
  }
  // run the filter
  AverageImageFilterType::Pointer avgFilter = AverageImageFilterType::New();
  for ( unsigned i = 0; i < numTestImages; ++i )
  {
    avgFilter->SetInput( i, inputImages[i] );
  }

  SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
  subFilter->SetInput1( testAvg );
  subFilter->SetInput2( avgFilter->GetOutput() );

  StatFilterType::Pointer statFilter = StatFilterType::New();
  statFilter->SetInput( subFilter->GetOutput() );

  statFilter->Update();

  std::cout << "Min(" << statFilter->GetMinimum() << ") Max(" << statFilter->GetMaximum() << ") Mean ("
            << statFilter->GetMean() << ") Sigma (" << statFilter->GetSigma() << ") Variance("
            << statFilter->GetVariance() << ") Sum (" << statFilter->GetSum() << std::endl;

  const double _eps( 0.00000001 );
  if ( std::abs( statFilter->GetMinimum() ) < _eps && std::abs( statFilter->GetMaximum() ) < _eps &&
       std::abs( statFilter->GetMean() ) < _eps )
  {
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
