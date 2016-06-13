/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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

  LabelImagePointerType PostProcessingANN( const std::string & continuousFilename, scalarType threshold );

  LabelImagePointerType PostProcessingRF( const std::string & labelImageFIlename, const unsigned char );

  LabelImagePointerType CombineLabel( LabelImagePointerType& resultLabel, LabelImagePointerType& currentLabel,
                                      const unsigned char binaryToLabelValue = 0 );

  LabelImagePointerType AmbiguousCountLabel( LabelImagePointerType& resultLabel, LabelImagePointerType& combinedLabel,
                                             LabelImagePointerType& currentLabel );

  LabelImagePointerType ThresholdImageAtLower( WorkingImagePointer& image, scalarType thresholdValue );

  LabelImagePointerType ThresholdImageAtUpper( WorkingImagePointer& image, scalarType thresholdValue );

  LabelImagePointerType ExtractLabel( const LabelImagePointerType& image, unsigned char thresholdValue );

  LabelImagePointerType GetOneConnectedRegion( LabelImagePointerType& image );

  LabelImagePointerType FillHole( LabelImagePointerType& mask, const unsigned char );

  LabelImagePointerType Closing( LabelImagePointerType& mask, const unsigned char );

  LabelImagePointerType Opening( LabelImagePointerType& image, const unsigned char );

  void WriteLabelMapToBinaryImages( const DataSet& subject, const LabelImagePointerType& labelMapImage );

  WorkingImagePointer ClipImageWithBinaryMask( WorkingImagePointer& image, WorkingImagePointer mask);

protected:
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(BRAINSCutApplyModel);
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
  cv::Ptr<OpenCVMLPType> m_openCVANN;

  //CvRTrees* m_openCVRandomForest;
  cv::Ptr<cv::ml::RTrees>  m_openCVRandomForest;

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
