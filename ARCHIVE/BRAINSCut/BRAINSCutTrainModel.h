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
#include "BRAINSCutExceptionStringHandler.h"
#include "BRAINSCutVectorTrainingSet.h"
#include "BRAINSCutDataHandler.h"
#include "TrainingPrameters.h"

#include "opencv2/core/core_c.h"

class BRAINSCutTrainModel
{
public:
  BRAINSCutTrainModel( BRAINSCutDataHandler & dataHandler );
  virtual ~BRAINSCutTrainModel();

  using CvRTrees = cv::ml::RTrees;

  /** train */
  void InitializeNeuralNetwork();

  void InitializeRandomForest();

  void InitializeTrainDataSet( bool doShuffle);

  void TrainANN();

  void TrainRandomForest();

  void TrainRandomForestAt( const int depth, const int numberOfTree );

  /** inline functions */
  inline void TrainWithUpdate(cv::Ptr<OpenCVMLPType> myTrainer, pairedTrainingSetType& currentTrainData);

  inline void SaveANNTrainModelAtIteration( cv::Ptr<OpenCVMLPType> myTrainer, unsigned int No);

  inline void SaveRFTrainModelAtIteration( cv::Ptr<cv::ml::RTrees> myTrainer, int depth, int NTrees);

  inline void recordRFTrainInformation( int depth, int nTree);

  inline void printANNTrainInformation( cv::Ptr<OpenCVMLPType> myTrainer, unsigned int No );

  void FillANNLayerStructureArray3D( int * const layer ) const;

  // TODO: REGINA all "Get" functions should be const
  /** setting function with net configuration */
  std::string GetModelBasename();

  std::string Getm_ANNVectorFilenamePrefix();

  void SetIteration();

  void SetEpochIteration();

  void SetDesiredError();

  void SetMaximumDataSize();

  void Setm_ANNHiddenNodesNumber();

  void SetActivatioinFunction();

  void SetModelBasename();

  /** default functions to set/get member variables */
  void SetIteration(unsigned int iteration);

  unsigned int GetIteration();

  void SetEpochIteration( unsigned int epochIteration);

  unsigned int GetEpochIteration();

  void SetDesiredError( float desiredError );

  float GetDesiredError();

  void SetMaximumDataSize( unsigned int maximumDataSize);

  unsigned int GetMaximumDataSize();

  void Setm_ANNHiddenNodesNumber( int hiddenNodes);

  int Getm_ANNHiddenNodesNumber();

  void SetActivationFunction( float slope, float minMax);

  float GetActivationSlope();

  float GetActivationMinMax();

  /** random trees */
  void  SetMaxDepth();

  void  SetMinSampleCount();

  void  SetUseSurrogates();

  void  SetCalcVarImportance();

  void  SetMaxTreeCount();

  void  SetRFErrorFilename();

  void  SetRFErrorFile();

  inline void appendToFile( std::string filename, std::string line);

private:
  // TODO:  REGINA:  These all need to be called with "m_" prefix
  BRAINSCutDataHandler m_myDataHandler;

  unsigned int m_trainIteration;
  unsigned int m_trainEpochIteration;
  float        m_trainDesiredError;
  unsigned int m_trainMaximumDataSize;

  int m_ANNHiddenNodesNumber;

  float m_activationSlope;
  float m_activationMinMax;

  /** random tree */
  int  m_trainMaxDepth;
  int  m_trainMinSampleCount;
  bool m_trainUseSurrogates;
  bool m_trainCalcVarImportance;
  int  m_trainMaxTreeCount;

  /** common paramters */
  std::string                 m_modelBasename;
  std::string                 m_ANNVectorFilenamePrefix;
  std::string                 m_RFErrorFilename;
  BRAINSCutVectorTrainingSet* m_trainingDataSet;

  CvMat * m_ANNLayerStructure;
};
