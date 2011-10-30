#ifndef BRAINSCutApplyModel_h
#define BRAINSCutApplyModel_h

#include "BRAINSCutPrimary.h"
#include "FeatureInputVector.h"

typedef itk::Image<unsigned char, DIMENSION> BinaryImageType;
typedef BinaryImageType::Pointer             BinaryTypePointer;

class BRAINSCutApplyModel : public BRAINSCutPrimary
{
public:
  BRAINSCutApplyModel( std::string netConfigurationFilename);

  void SetApplyDataSetFromNetConfiguration();

  void Apply();

  void ApplyOnSubject( DataSet& subject);

  void SetANNModelFilenameFromNetConfiguration();

  void ReadANNModelFile();

  BinaryTypePointer PostProcessingOfANNContinuousImage( std::string continuousFilname );

  void SetANNOutputThresholdFromNetConfiguration();

  BinaryTypePointer ThresholdImage( WorkingImagePointer image );

  BinaryTypePointer GetOneConnectedRegion( BinaryTypePointer image );

private:
  NetConfiguration::ApplyDataSetListType applyDataSetList;

  bool        normalization;
  std::string ANNModelFilename;

  scalarType annOutputThreshold;

  OpenCVMLPType * openCVANN;

  /* inline functions */

  inline void PredictROI( InputVectorMapType& roiInputFeatureVector, PredictValueMapType& resultOutputVector,
                          const unsigned int roiNumber, const unsigned int inputVectorSize);

  inline scalarType * GetArrayFromVector( scalarType array[], InputVectorType & vector, unsigned int inputVectorSize);

  inline void         GetOpenCVMatrixFromArray( matrixType & matrix, scalarType array[], unsigned int inputVectorSize);

  inline void WritePredictROIProbabilityBasedOnReferenceImage( const PredictValueMapType& predictedOutput,
                                                               const WorkingImagePointer& referenceImage,
                                                               const WorkingImagePointer& roi,
                                                               const std::string imageFIlename );

  inline std::string GetSubjectOutputDirectory( DataSet& subject);

  inline std::string GetContinuousPredictionFilename( DataSet& subject, std::string currentROIName);

  inline std::string GetOutputROIFilename( DataSet& subject, std::string currentROIName);
};

#endif
