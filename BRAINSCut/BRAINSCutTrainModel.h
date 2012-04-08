#include "BRAINSCutExceptionStringHandler.h"
#include "BRAINSCutVectorTrainingSet.h"
#include "TrainingPrameters.h"

class BRAINSCutTrainModel : public BRAINSCutPrimary
{
public:
  BRAINSCutTrainModel( std::string netConfigurationFilename );

  /** train */
  void InitializeNeuralNetwork();

  void InitializeRandomForest();

  void InitializeTrainDataSet( bool doShuffle);

  void TrainANN();

  void TrainRandomForest();

  void TrainRandomForestAt( const int depth, const int numberOfTree );

  /** inline functions */
  inline void TrainWithUpdate(neuralNetType& myTrainer, bool update, pairedTrainingSetType& currentTrainData);

  inline void SaveANNTrainModelAtIteration( neuralNetType& myTrainer, unsigned int No);

  inline void SaveRFTrainModelAtIteration( CvRTrees& myTrainer, int depth, int NTrees);

  inline void writeRFTrainInformation( CvRTrees& myTrainer, int depth, int nTree);

  inline void printANNTrainInformation( neuralNetType& myTrainer, unsigned int No );

  inline int * GetANNLayerStructureArray();

  /** setting function with net configuration */
  std::string GetModelBasename();

  std::string GetANNVectorFilenamePrefix();

  void SetIterationFromNetConfiguration();

  void SetEpochIterationFromNetConfiguration();

  void SetDesiredErrorFromNetConfiguration();

  void SetMaximumDataSizeFromNetConfiguration();

  void SetANNHiddenNodesNumberFromNetConfiguration();

  void SetActivatioinFunctionFromNetConfiguration();

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

  void SetANNHiddenNodesNumber( int hiddenNodes);

  int GetANNHiddenNodesNumber();

  void SetActivationFunction( float slope, float minMax);

  float GetActivationSlope();

  float GetActivationMinMax();

  /** random trees */
  void  SetMaxDepthFromNetConfiguration();

  void  SetMinSampleCountFromNetConfiguration();

  void  SetUseSurrogatesFromNetConfiguration();

  void  SetCalcVarImportanceFromNetConfiguration();

  void  SetMaxTreeCountFromNetConfiguration();

  void  SetRFErrorFilename();

  void  SetRFErrorFile();

  inline void appendToFile( std::string filename, std::string line);

private:
  /* train parameters */
  TrainingParameters *TrainNetConfiguration;

  unsigned int trainIteration;
  unsigned int trainEpochIteration;
  float        trainDesiredError;
  unsigned int trainMaximumDataSize;

  int ANNHiddenNodesNumber;

  float activationSlope;
  float activationMinMax;

  /** random tree */
  int  trainMaxDepth;
  int  trainMinSampleCount;
  bool trainUseSurrogates;
  bool trainCalcVarImportance;
  int  trainMaxTreeCount;

  /** common paramters */
  std::string                 modelBasename;
  std::string                 ANNVectorFilenamePrefix;
  std::string                 RFErrorFilename;
  BRAINSCutVectorTrainingSet* trainingDataSet;

  matrixType ANNLayerStructure;
};
