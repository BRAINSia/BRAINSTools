#include <iostream>
#include <itkImage.h>
#include <itkImageFileReader.h>

#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkScalarImageToTextureFeaturesFilter.h>

//#    1) itkEllipseSpatialObject
//#    2) itkSpatialObjectToImageFilter
//#    3) itkSimpleContourExtractorImageFilter

#include <iostream>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkVersorRigid3DTransform.h>
#include <itkImageMaskSpatialObject.h>
#include <itkCastImageFilter.h>
#include <string>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkHistogramToTextureFeaturesFilter.h>

int main(int, char * *)
{
  typedef float                           PixelType;
  typedef itk::Image<PixelType, 3>        ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  std::string filename ("/Users/johnsonhj/Dropbox/DATA/ANTSDENOISE_FAILURE/input.nii.gz");
  ReaderType::Pointer myReader = ReaderType::New();
  myReader->SetFileName(filename);
  myReader->Update();
  ImageType::Pointer myImage = myReader->GetOutput();


#if 1
  typedef itk::Statistics::ScalarImageToTextureFeaturesFilter< ImageType> TextureFilterType;
  TextureFilterType::Pointer texFilter = TextureFilterType::New();
  texFilter->SetInput(myImage);
  texFilter->SetMaskImage( myImage );
  texFilter->SetFastCalculations( false );

  //Test Set/Get Requested features
  typedef TextureFilterType::TextureFeaturesFilterType   TextureFeaturesFilterType;

  TextureFilterType::FeatureNameVectorPointer requestedFeatures =
      TextureFilterType::FeatureNameVector::New();

  requestedFeatures->push_back(TextureFeaturesFilterType::Inertia);
  requestedFeatures->push_back(TextureFeaturesFilterType::ClusterShade);
  texFilter->SetRequestedFeatures(requestedFeatures);

  const TextureFilterType::FeatureNameVector* requestedFeatures2 =
      texFilter->GetRequestedFeatures();

  TextureFilterType::FeatureNameVector::ConstIterator fIt;

  fIt = requestedFeatures2->Begin();
  bool passed = true;
  if ( fIt.Value() != TextureFeaturesFilterType::Inertia )
  {
    std::cerr << "Requested feature name not correctly set" << std::endl;
    passed = false;
  }
  fIt++;

  if ( fIt.Value() != TextureFeaturesFilterType::ClusterShade )
  {
    std::cerr << "Requested feature name not correctly set" << std::endl;
    passed = false;
  }

  texFilter->Update();
  texFilter->GetOutputNames();
  texFilter->GetOutputs();


#endif



  std::cout << "FINSIHED" << std::endl;
  return EXIT_SUCCESS;
}