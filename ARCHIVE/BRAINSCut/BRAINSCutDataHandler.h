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
#ifndef BRAINSCutDataHandler_h
#define BRAINSCutDataHandler_h

#include "BRAINSCutUtilities.h"

/*
 * BRAINSCut Primary Class Starts here
 */

class BRAINSCutDataHandler
{
public:
  BRAINSCutDataHandler();
  BRAINSCutDataHandler(const std::string & modelConfigurationFilenameFilename);

  void                     SetNetConfiguration();

  BRAINSCutConfiguration * GetNetConfiguration();

  void        SetNetConfigurationFilename(const std::string & filename);

  std::string GetNetConfigurationFilename();

  void SetAtlasDataSet();

  void SetAtlasImage();

  void SetRegionsOfInterest();

  void         SetRegistrationParameters();

  std::string  GetRegistrationID();

  void         SetRegistrationImageTypeToUse( std::string type );

  std::string  GetRegistrationImageTypeToUse();

  void         SetRhoPhiTheta();

  void         SetGradientSize();

  unsigned int GetGradientSize();

  void SetTrainingVectorConfiguration();

  std::string GetModelBaseName();

  void        SetANNModelFilenameAtIteration( const int iteration);

  void        SetANNTestingSSEFilename();

  std::string GetANNTestingSSEFilename();

  std::string GetANNModelFilenameAtIteration( const int iteration);

  std::string GetANNModelFilename();

  void        SetRandomForestModelFilename( int depth, int nTree );

  void        SetRandomForestModelFilename( std::string name );

  std::string GetRandomForestModelFilename();

  std::string GetRFModelFilename( unsigned int depth, unsigned int NTrees);

  DataSet::StringVectorType GetROIIDsInOrder() const;

  void        SetTrainVectorFilename();

  std::string GetTrainVectorFilename();

  void                  GetDeformedSpatialLocationImages(std::map<std::string,
                                                                  WorkingImagePointer>& warpedSpatialLocationImages,
                                                         DataSet& subject );

  void                  ReadImagesOfSubjectInOrder( WorkingImageVectorType& subjectImageList, DataSet& subject);

  void                  GetDeformedROIs( std::map<std::string,
                                                  WorkingImagePointer>& deformedROIs, DataSet& subject );

  std::string           GetNormalizationMethod();

  void                  SetNormalization();

  std::string           GetAtlasFilename();

  std::string           GetAtlasBinaryFilename();

  int                   GetROIAutoDilateSize();

  unsigned int          GetROICount();

  WorkingImagePointer   GetAtlasImage();

  ProbabilityMapList *  GetROIDataList();

  DataSet *      GetAtlasDataSet();

  BRAINSCutConfiguration::ApplyDataSetListType GetApplyDataSet();

  BRAINSCutConfiguration::TrainDataSetListType GetTrainDataSet();

  int                   GetTrainIteration();

  scalarType            GetANNOutputThreshold();

  scalarType            GetGaussianSmoothingSigma();

  std::string           GetSubjectToAtlasRegistrationFilename( DataSet& subject);

  std::string           GetAtlasToSubjectRegistrationFilename( DataSet& subject);

  void         SetTrainConfiguration( std::string trainParamterName );

  unsigned int GetEpochIteration();

  float        GetDesiredError();

  unsigned int GetMaximumDataSize();

  int          GetANNHiddenNodesNumber();

  float        GetActivationFunctionSlope();

  float        GetActivationFunctionMinMax();

  int          GetMaxDepth();

  int          GetMinSampleCount();

  bool         GetUseSurrogates();

  bool         GetCalcVarImportance();

  int          GetMaxTreeCount();

  WorkingImagePointer GetCandidateRegion( DataSet& subject) const;

protected:
  TrainingVectorConfigurationType * trainingVectorConfiguration;

  /* train parameters */
  TrainingParameters *TrainConfiguration;

  /** atlas data set*/
  DataSet *           m_atlasDataSet;
  std::string         m_atlasFilename;
  std::string         m_atlasBinaryFilename;
  WorkingImagePointer m_atlasImage;

  /**ProbabilityMaps*/
  ProbabilityMapList *      m_roiDataList;
  DataSet::StringVectorType m_roiIDsInOrder;
  unsigned int              roiCount;
  bool                      ROIRegistrationToSubject;

  /** registration data set */
  RegistrationConfigurationParser * registrationParser;
  std::string                       registrationImageTypeToUse;
  std::string                       registrationID;
  int                               roiAutoDilateSize;

  /** Spatial Coordinate System Images*/
  WorkingImagePointer m_rho;
  WorkingImagePointer m_phi;
  WorkingImagePointer m_theta;

  unsigned int m_gradientSize;

  /** vector file name */
  std::string m_trainVectorFilename;
  std::string m_normalization;

  /** model name **/
  std::string ANNModelFilename;
  std::string RandomForestModelFilename;
  std::string ANNTestingSSEFilename;
private:
  std::string              myConfigurationFilename;
  BRAINSCutConfiguration * myConfiguration;

  WorkingImageType GetDeformedImage( WorkingImageType image);
};
#endif
