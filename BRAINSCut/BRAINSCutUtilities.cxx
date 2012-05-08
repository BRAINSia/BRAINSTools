#include "BRAINSCutUtilities.h"

//
// read/warp image
//
WorkingImagePointer
SmoothImage( const WorkingImagePointer image, const float GaussianValue)
{
  if( GaussianValue < 0 + FLOAT_TOLERANCE )
    {
    std::cout << "Gaussian value is less than tolerance. "
              << "No smoothing occurs at this time"
              << std::endl;
    return image;
    }
  /*std::cout<<"Smooth Image with Gaussian value of :: "
           << GaussianValue
           <<std::endl;*/
  typedef itk::SmoothingRecursiveGaussianImageFilter<WorkingImageType, WorkingImageType> SmoothingFilterType;
  SmoothingFilterType::Pointer smoothingFilter = SmoothingFilterType::New();

  smoothingFilter->SetInput( image);
  smoothingFilter->SetSigma( GaussianValue );

  smoothingFilter->Update();

  return smoothingFilter->GetOutput();
}

WorkingImagePointer ReadImageByFilename( const std::string  filename )
{
  WorkingImagePointer readInImage;

  ReadInImagePointer inputImage = itkUtil::ReadImage<ReadInImageType>(filename.c_str() );

  readInImage = itkUtil::ScaleAndCast<ReadInImageType,
                                      WorkingImageType>(inputImage,
                                                        ZeroPercentValue,
                                                        HundredPercentValue);
  return readInImage;
}

/* inline functions */

DisplacementFieldType::Pointer GetDeformationField( std::string filename)
{
  const bool useTransform( filename.find(".mat") != std::string::npos );

  if( useTransform )
    {
    return NULL;
    }
  typedef itk::ImageFileReader<DisplacementFieldType> DeformationReaderType;
  DeformationReaderType::Pointer deformationReader = DeformationReaderType::New();
  deformationReader->SetFileName( filename );
  deformationReader->Update();

  return deformationReader->GetOutput();
}

GenericTransformType::Pointer GetGenericTransform( std::string filename)
{
  const bool useDeformation( filename.find(".mat") == std::string::npos );

  if( useDeformation )
    {
    std::cout << "return null deformation" << std::endl;
    return NULL;
    }
  return itk::ReadTransformFromDisk( filename );
}
