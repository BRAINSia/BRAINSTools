#ifndef BRAINSCutPrimary_h
#define BRAINSCutPrimary_h

#include "BRAINSCutConfiguration.h"

#include "TrainingVectorConfigurationType.h"
#include "TrainingPrameters.h"
#include "ApplyModel.h"
#include "itkIO.h"

#include "GenericTransformImage.h"

#include <itkSmoothingRecursiveGaussianImageFilter.h>

/** include opencv library */
#include "ml.h"
#include "cxcore.h"
#include "cv.h"

typedef CvANN_MLP_Revision OpenCVMLPType;
#include <stdint.h>

/** Training data set definition */

typedef float   scalarType;
typedef CvMat * matrixType;

struct pairedTrainingSetType
  {
  matrixType pairedInput;
  matrixType pairedOutput;
  matrixType pairedOutputRF;
  unsigned int size;
  };

/*
 * constant
 */
static const float        HundredPercentValue = 1.0F;
static const float        ZeroPercentValue = 0.0F;
static const unsigned int LineGuardSize = 1;
static const scalarType   LineGuard = 1234567.0;
static const float        FLOAT_TOLERANCE = 0.01;

/*
* Image Definitions
*/
const unsigned char DIMENSION = 3;

typedef double                                 ReadInPixelType;
typedef itk::Image<ReadInPixelType, DIMENSION> ReadInImageType;
typedef ReadInImageType::Pointer               ReadInImagePointer;

typedef float                                   WorkingPixelType;
typedef itk::Image<WorkingPixelType, DIMENSION> WorkingImageType;
typedef WorkingImageType::Pointer               WorkingImagePointer;

typedef std::vector<WorkingImagePointer> WorkingImageVectorType;

/* Deformations */
typedef float                                         DeformationScalarType;
typedef itk::Vector<DeformationScalarType, DIMENSION> DeformationPixelType;
typedef itk::Image<DeformationPixelType, DIMENSION>   DisplacementFieldType;

typedef WorkingImageType::IndexType WorkingIndexType;

typedef std::vector<WorkingPixelType> InputVectorType;
typedef std::vector<WorkingPixelType> OutputVectorType;

typedef std::map<int, InputVectorType>  InputVectorMapType; // < index ,feature vector > pair
typedef std::map<int, OutputVectorType> OutputVectorMapType;
typedef std::map<int, scalarType>       PredictValueMapType;

const WorkingImageType::IndexType ConstantHashIndexSize = {{255, 255, 255}};

std::string GetAtlasToSubjectRegistrationFilename( DataSet& subject);

std::string GetSubjectToAtlasRegistrationFilename( DataSet& subject);

//
// read/warp image
//
static
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

static
WorkingImagePointer
ReadImageByFilename( const std::string  filename )
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

static
inline
DisplacementFieldType::Pointer
GetDeformationField( std::string filename)
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

static
inline
GenericTransformType::Pointer
GetGenericTransform( std::string filename)
{
  const bool useDeformation( filename.find(".mat") == std::string::npos );

  if( useDeformation )
    {
    std::cout << "return null deformation" << std::endl;
    return NULL;
    }
  return itk::ReadTransformFromDisk( filename );
}

#endif
