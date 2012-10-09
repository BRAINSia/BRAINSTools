#ifndef BRAINSCutApplyModel_h
#define BRAINSCutApplyModel_h

#include "BRAINSCutDataHandler.h"
#include "FeatureInputVector.h"

typedef itk::Image<unsigned char, DIMENSION> BinaryImageType;
typedef BinaryImageType::Pointer             BinaryImagePointer;

class BRAINSCutApplyModel
{
public:
  BRAINSCutApplyModel()
  {
  };
  BRAINSCutApplyModel( BRAINSCutDataHandler& dataHandler );
  ~BRAINSCutApplyModel();

  void Apply();

  void ApplyOnSubject( DataSet& subject);

  void SetComputeSSE( const bool sse );

  void SetMethod( std::string inputMethod);

  void SetNumberOfTrees( const int numberOfTree );

  void SetDepthOfTree( const int depth );

  void ReadANNModelFile();

  void ReadRandomForestModelFile();

  BinaryImagePointer PostProcessingANN( std::string continuousFilename, scalarType threshold );

  BinaryImagePointer PostProcessingRF( std::string labelImageFilename );

  BinaryImagePointer ThresholdImageAtLower( WorkingImagePointer& image, scalarType thresholdValue );

  BinaryImagePointer ThresholdImageAtUpper( WorkingImagePointer& image, scalarType thresholdValue );

  BinaryImagePointer ExtractLabel( BinaryImagePointer image, unsigned char thresholdValue );

  BinaryImagePointer GetOneConnectedRegion( BinaryImagePointer image );

  BinaryImagePointer FillHole( BinaryImagePointer mask);

private:
  BRAINSCutDataHandler                         m_myDataHandler;
  BRAINSCutConfiguration::ApplyDataSetListType applyDataSetList;

  std::string method;
  bool        normalization;
  bool        computeSSE;
  int         trainIteration;

  int numberOfTrees;
  int depthOfTree;

  std::fstream ANNTestingSSEFileStream;

  scalarType      annOutputThreshold;
  scalarType      gaussianSmoothingSigma;
  OpenCVMLPType * m_openCVANN;

  CvRTrees openCVRandomForest;

  /* private functions  */
  std::string GetANNModelBaseName();

  float ComputeSSE(const PredictValueMapType& predictedOutputVector, const std::string roiReferenceFilename );

  /* inline functions */

  void PredictROI( InputVectorMapType& roiInputFeatureVector, PredictValueMapType& resultOutputVector,
                   const unsigned int roiNumber, const unsigned int inputVectorSize) const;

  inline void WritePredictROIProbabilityBasedOnReferenceImage( const PredictValueMapType& predictedOutput,
                                                               const WorkingImagePointer& referenceImage,
                                                               const WorkingImagePointer& roi,
                                                               const std::string imageFilename,
                                                               const WorkingPixelType labelValue = HundredPercentValue);

  inline std::string GetSubjectOutputDirectory( DataSet& subject);

  inline std::string GetContinuousPredictionFilename( DataSet& subject, std::string currentROIName);

  inline std::string GetROIVolumeName( DataSet& subject, std::string currentROIName);
};

#endif
