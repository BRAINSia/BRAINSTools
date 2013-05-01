#ifndef BRAINSCutUtilities_h
#define BRAINSCutUtilities_h

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

#include "opencv2/flann/flann.hpp"

// typedef CvANN_MLP_Revision OpenCVMLPType;
typedef CvANN_MLP OpenCVMLPType;

/** Training data set definition */

typedef float scalarType;

struct pairedTrainingSetType
  {
  CvMat * pairedInput;
  CvMat * pairedOutput;
  CvMat * pairedOutputRF;
  unsigned int size;
  };

/* normalization type */
enum FeatureNormalizationMethodEnum
  {
  Linear,
  Sigmoid,
  DoubleSigmoid,
  zScore
  };
/*
 * constant
 */
static const float        HundredPercentValue = 1.0F;
static const float        ZeroPercentValue = -1.0F;
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
typedef double                                        DeformationScalarType;
typedef itk::Vector<DeformationScalarType, DIMENSION> DeformationPixelType;
typedef itk::Image<DeformationPixelType, DIMENSION>   DisplacementFieldType;

typedef WorkingImageType::IndexType WorkingIndexType;

typedef std::vector<WorkingPixelType> InputVectorType;
typedef std::vector<WorkingPixelType> OutputVectorType;

// HACK TODO:  Regina int below should be unsigned int to avoid negative index numbers
typedef std::map<int, InputVectorType>  InputVectorMapType; // < index ,feature vector > pair
typedef std::map<int, OutputVectorType> OutputVectorMapType;
typedef std::map<int, scalarType>       PredictValueMapType;

std::string GetAtlasToSubjectRegistrationFilename( DataSet& subject);

std::string GetSubjectToAtlasRegistrationFilename( DataSet& subject);

WorkingImagePointer SmoothImage( const WorkingImagePointer& image, const float GaussianValue);

WorkingImagePointer ReadImageByFilename( const std::string  & filename );

DisplacementFieldType::Pointer GetDeformationField( std::string filename);

GenericTransformType::Pointer GetGenericTransform( std::string filename);

#endif
