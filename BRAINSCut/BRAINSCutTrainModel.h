#include "BRAINSCutExceptionStringHandler.h"
#include "BRAINSCutVectorTrainingSet.h"
#include "BRAINSCutDataHandler.h"
#include "TrainingPrameters.h"

class BRAINSCutTrainModel
{
public:
  BRAINSCutTrainModel( BRAINSCutDataHandler & dataHandler );

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

  void SetIteration();

  void SetEpochIteration();

  void SetDesiredError();

  void SetMaximumDataSize();

  void SetANNHiddenNodesNumber();

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

  void SetANNHiddenNodesNumber( int hiddenNodes);

  int GetANNHiddenNodesNumber();

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
  BRAINSCutDataHandler myDataHandler;

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
