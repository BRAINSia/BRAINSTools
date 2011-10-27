#include "BRAINSCutExceptionStringHandler.h"
#include "BRAINSCutVectorTrainingSet.h"
#include "ANNParams.h"

class BRAINSCutTrainModel : public BRAINSCutPrimary
{
public:
  BRAINSCutTrainModel( std::string netConfigurationFilename );

  /** train */
  void InitializeNeuralNetwork();

  void InitializeTrainDataSet();

  void Train();

  /** inline functions */
  inline void TrainWithUpdate(neuralNetType& myTrainer, bool update, pairedTrainingSetType& currentTrainData);

  inline void SaveTrainModelAtIteration( neuralNetType& myTrainer, unsigned int No);

  inline void printTrainInformation( neuralNetType& myTrainer, unsigned int No );

  inline int * GetANNLayerStructureArray();

  /** setting function with net configuration */
  std::string GetANNModelFilenamePrefix();

  std::string GetANNVectorFilenamePrefix();

  void SetIterationFromNetConfiguration();

  void SetEpochIterationFromNetConfiguration();

  void SetDesiredErrorFromNetConfiguration();

  void SetMaximumDataSizeFromNetConfiguration();

  void SetANNHiddenNodesNumberFromNetConfiguration();

  void SetActivatioinFunctionFromNetConfiguration();

  void SetANNModelFilenamePrefix();

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

private:
  /* train parameters */
  ANNParams *ANNParameterNetConfiguration;

  unsigned int trainIteration;
  unsigned int trainEpochIteration;
  float        trainDesiredError;
  unsigned int trainMaximumDataSize;

  int ANNHiddenNodesNumber;

  float activationSlope;
  float activationMinMax;

  std::string                 ANNModelFilenamePrefix;
  std::string                 ANNVectorFilenamePrefix;
  BRAINSCutVectorTrainingSet* trainingDataSet;

  matrixType ANNLayerStructure;
};
