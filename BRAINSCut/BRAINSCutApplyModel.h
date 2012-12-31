#ifndef BRAINSCutApplyModel_h
#define BRAINSCutApplyModel_h

#include "BRAINSCutDataHandler.h"
#include "FeatureInputVector.h"

typedef itk::Image<unsigned char, DIMENSION> LabelImageType;
typedef LabelImageType::Pointer              LabelImagePointerType;

class BRAINSCutApplyModel
{
public:
  BRAINSCutApplyModel( BRAINSCutDataHandler& dataHandler );
  BRAINSCutApplyModel();
  ~BRAINSCutApplyModel();

  void Apply();

  void ApplyOnSubject( DataSet& subject);

  void SetComputeSSE( const bool sse );

  void SetMethod( std::string inputMethod);

  void SetNumberOfTrees( const int numberOfTree );

  void SetDepthOfTree( const int depth );

  void ReadANNModelFile();

  void ReadRandomForestModelFile();

  LabelImagePointerType PostProcessingANN( std::string continuousFilename, scalarType threshold );

  LabelImagePointerType PostProcessingRF( std::string labelImageFIlename, const unsigned char = 1 );

  LabelImagePointerType CombineLabel( LabelImagePointerType& resultLabel, LabelImagePointerType& currentLabel,
                                      const unsigned char binaryToLabelValue = 0 );

  LabelImagePointerType AmbiguousCountLabel( LabelImagePointerType& resultLabel, LabelImagePointerType& combinedLabel,
                                             LabelImagePointerType& currentLabel );

  LabelImagePointerType ThresholdImageAtLower( WorkingImagePointer& image, scalarType thresholdValue );

  LabelImagePointerType ThresholdImageAtUpper( WorkingImagePointer& image, scalarType thresholdValue );

  LabelImagePointerType ExtractLabel( const LabelImagePointerType& image, unsigned char thresholdValue );

  LabelImagePointerType GetOneConnectedRegion( LabelImagePointerType& image, const unsigned char = 1 );

  LabelImagePointerType FillHole( LabelImagePointerType& mask, const unsigned char = 1 );

  LabelImagePointerType Closing( LabelImagePointerType& mask, const unsigned char = 1 );

  LabelImagePointerType Opening( LabelImagePointerType& image, const unsigned char = 1);

  void WriteLabelMapToBinaryImages( const DataSet& subject, const LabelImagePointerType& labelMapImage );

  WorkingImagePointer ClipImageWithBinaryMask( WorkingImagePointer& image, WorkingImagePointer mask);

protected:
private:
  BRAINSCutDataHandler*                        m_myDataHandler;
  BRAINSCutConfiguration::ApplyDataSetListType m_applyDataSetList;

  std::string m_method;
  std::string m_normalization;
  bool        m_computeSSE;
  int         m_trainIteration;

  int m_numberOfTrees;
  int m_depthOfTree;

  std::fstream m_ANNTestingSSEFileStream;

  scalarType      m_annOutputThreshold;
  scalarType      m_gaussianSmoothingSigma;
  OpenCVMLPType * m_openCVANN;

  CvRTrees* m_openCVRandomForest;

  /* private functions  */
  std::string GetANNModelBaseName();

  float ComputeSSE(const PredictValueMapType& predictedOutputVector, const std::string & roiReferenceFilename );

  /* inline functions */

  void PredictROI( InputVectorMapType& roiInputFeatureVector, PredictValueMapType& resultOutputVector,
                   const unsigned int roiNumber, const unsigned int inputVectorSize) const;

  inline void WritePredictROIProbabilityBasedOnReferenceImage( const PredictValueMapType& predictedOutput,
                                                               const WorkingImagePointer& referenceImage,
                                                               const WorkingImagePointer& roi,
                                                               const std::string & imageFilename, const WorkingPixelType & labelValue =
                                                                 HundredPercentValue);

  inline std::string GetSubjectOutputDirectory( const DataSet& subject);

  inline std::string GetLabelMapFilename( const DataSet& subject );

  inline std::string GetContinuousPredictionFilename( const DataSet& subject, const std::string & currentROIName);

  inline std::string GetROIVolumeName( const DataSet& subject, const std::string & currentROIName);
};

#endif
