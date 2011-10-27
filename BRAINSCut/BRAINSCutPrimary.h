#ifndef BRAINSCutPrimary_h
#define BRAINSCutPrimary_h

#include "NetConfiguration.h"
#include "NeuralParams.h"
#include "itkIO.h"

#include "GenericTransformImage.h"
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
  unsigned int size;
  };

/*
 * constant
 */
static const float        HundreadPercentValue = 1.0F;
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
typedef itk::Image<DeformationPixelType, DIMENSION>   DeformationFieldType;

typedef WorkingImageType::IndexType WorkingIndexType;

typedef std::vector<WorkingPixelType> InputVectorType;
typedef std::vector<WorkingPixelType> OutputVectorType;

typedef map<int, InputVectorType>  InputVectorMapType; // < index ,feature vector > pair
typedef map<int, OutputVectorType> OutputVectorMapType;
typedef map<int, scalarType>       PredictValueMapType;

const WorkingImageType::IndexType ConstantHashIndexSize = {{255, 255, 255}};
/*
 * BRAINSCut Primary Class Starts here
 */

class BRAINSCutPrimary
{
public:
  BRAINSCutPrimary(std::string netConfigurationFilename);

  void SetNetConfiguration();

  NetConfiguration * GetNetConfiguration();

  void SetNetConfigurationFilename(const std::string filename);

  std::string GetNetConfigurationFilename();

  /* Atlas(template) related */
  void SetAtlasDataSet();

  void SetAtlasImage();

  void SetRegionsOfInterestFromNetConfiguration();

  void SetRegistrationParametersFromNetConfiguration();

  void SetRhoPhiThetaFromNetConfiguration();

  void SetANNModelConfiguration();

  void SetGradientSizeFromNetConfiguration();

  WorkingImagePointer ReadImageByFilename( const std::string  filename );

  WorkingImagePointer WarpImageByFilenames( const std::string & deformationFilename, const std::string & inputFilename,
                                            const std::string & referenceFilename );

  /** Get Function */
  DataSet::StringVectorType GetROIIDsInOrder();

  void   GetDeformedSpatialLocationImages( std::map<std::string, WorkingImagePointer>& warpedSpatialLocationImages,
                                           DataSet& subject );

  void GetImagesOfSubjectInOrder( WorkingImageVectorType& subjectImageList, DataSet& subject);

  void GetDeformedROIs( std::map<std::string, WorkingImagePointer>& deformedROIs, DataSet& subject );

  bool GetNormalizationFromNetConfiguration();

  /* Deformation Functions */
  inline DeformationFieldType::Pointer GetDeformationField( std::string filename);

  inline GenericTransformType::Pointer GetGenericTransform( std::string filename);

  inline std::string GetAtlasToSubjectRegistrationFilename( DataSet& subject);

protected:

  NetConfiguration BRAINSCutNetConfiguration;
  NeuralParams *   annModelConfiguration;

  /** atlas data set*/
  DataSet *           atlasDataSet;
  WorkingImagePointer atlasImage;

  /**ProbabilityMaps*/
  ProbabilityMapList *      roiDataList;
  DataSet::StringVectorType roiIDsInOrder;;
  unsigned int              roiCount;

  /** registration data set */
  RegistrationConfigurationParser * registrationParser;
  std::string                       registrationImageTypeToUse;
  std::string                       registrationID;

  /** Spatial Coordinate System Images*/
  WorkingImagePointer rho;
  WorkingImagePointer phi;
  WorkingImagePointer theta;

  unsigned int gradientSize;
private:
  WorkingImageType GetDeformedImage( WorkingImageType image);

  std::string NetConfigurationFilename;
};
#endif
