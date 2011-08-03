#include "itkImage.h"
#include "itkBlendImageFilter.h"
#include "itkRandomImageSource.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

int main(int, char * *)
{
  typedef itk::Image<float, 2> ImageType;

  std::cout << "Create input image using RandomImageSource" << std::endl;
  ImageType::Pointer images[2];
  for( unsigned i = 0; i < 2; i++ )
    {
    typedef itk::RandomImageSource<ImageType> SourceType;
    SourceType::Pointer      source = SourceType::New();
    ImageType::SizeValueType size[2] = {64, 64};
    source->SetSize( size );
    source->SetMin(0.0);
    source->SetMax(1.0);
    source->Update();
    images[i] = source->GetOutput();
    }
  typedef itk::BlendImageFilter<ImageType, ImageType>
    BlendImageFilterType;
  BlendImageFilterType::Pointer filter =
    BlendImageFilterType::New();
  filter->SetInput1(images[0]);
  filter->SetInput2(images[1]);
  filter->SetBlend1(0.2);
  filter->SetBlend2(0.8);
  filter->Update();
  ImageType::Pointer blendImage = filter->GetOutput();
  itk::ImageRegionConstIterator<ImageType>
  it1(images[0],
      images[0]->GetLargestPossibleRegion() ),
  it2(images[1],
      images[1]->GetLargestPossibleRegion() ),
  itBlend(blendImage,
          blendImage->GetLargestPossibleRegion() );
  for( ; !it1.IsAtEnd() && !it2.IsAtEnd() && !itBlend.IsAtEnd();
       ++it1, ++it2, ++itBlend )
    {
    float blend = (it1.Get() * 0.2) + (it2.Get() * 0.8);
    if( vcl_fabs(blend - itBlend.Get() ) > 0.0001 )
      {
      std::cerr << "Expected " << blend << " found " << itBlend.Get()
                << std::endl;
      exit(1);
      }
    }
  exit(0);
}
