#include "BRAINSCutTrainModel.h"
#include "NeuralParams.h"
#include "fstream.h"

BRAINSCutTrainModel
::BRAINSCutTrainModel( std::string netConfigurationFIlename)
  : BRAINSCutPrimary( netConfigurationFIlename ),
  trainIteration(0),
  trainEpochIteration(0),
  trainDesiredError(0.0),
  trainMaximumDataSize(0),
  ANNHiddenNodesNumber(0),
  activationSlope(0),
  activationMinMax(0)
{
}

/** train */
void
BRAINSCutTrainModel
::InitializeTrainDataSet( bool doShuffle )
{
  const std::string Local_ANNVectorFilenamePrefix = this->GetANNVectorFilenamePrefix();

  trainingDataSet = new BRAINSCutVectorTrainingSet( Local_ANNVectorFilenamePrefix);
  try
    {
    trainingDataSet->ReadHeaderFileInformation();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error();
    exit(EXIT_FAILURE);
    }
  trainingDataSet->SetRecordSize();
  trainingDataSet->SetBufferRecordSize();
  if( doShuffle )
    {
    trainingDataSet->ShuffleVectors();
    }
  if( trainingDataSet->GetTotalVectorSize() > (int)trainMaximumDataSize )
    {
    unsigned int numberOfSubSet =
      (float)trainingDataSet->GetTotalVectorSize() / (float)trainMaximumDataSize;
    numberOfSubSet = ceil(numberOfSubSet) + 1;
    std::cout << "Divide subset into " << numberOfSubSet << std::endl;
    trainingDataSet->SetNumberOfSubSet( numberOfSubSet );
    }
  else
    {
    trainingDataSet->SetNumberOfSubSet( 1 );
    }
}

/** initialization for the model */
void
BRAINSCutTrainModel
::InitializeNeuralNetwork()
{
  SetANNModelConfiguration();
  TrainNetConfiguration = BRAINSCutNetConfiguration.Get<TrainingParameters>("ANNParameters");
  ANNLayerStructure = cvCreateMat( 1, 3, CV_32SC1);
  SetMaximumDataSizeFromNetConfiguration();

  SetIterationFromNetConfiguration();
  SetEpochIterationFromNetConfiguration();
  SetDesiredErrorFromNetConfiguration();
  SetANNHiddenNodesNumberFromNetConfiguration();
  SetActivatioinFunctionFromNetConfiguration();
  SetModelBasename();
}

void
BRAINSCutTrainModel
::InitializeRandomForest()
{
  TrainNetConfiguration = BRAINSCutNetConfiguration.Get<TrainingParameters>("RandomForestParameters");
  SetANNModelConfiguration();
  SetMaxDepthFromNetConfiguration();
  SetMinSampleCountFromNetConfiguration();
  SetUseSurrogatesFromNetConfiguration();
  SetCalcVarImportanceFromNetConfiguration();
  SetMaxTreeCountFromNetConfiguration();

  SetModelBasename();
  SetRFErrorFilename();
  SetRFErrorFile();
}

inline void
BRAINSCutTrainModel
::TrainWithUpdate( neuralNetType& myTrainer, bool update, pairedTrainingSetType& currentTrainData )
{
  int updateOption = 0;

  if( update )
    {
    updateOption = CvANN_MLP::UPDATE_WEIGHTS;
    }
  // TODO change subset number properly

  myTrainer.train( currentTrainData.pairedInput,
                   currentTrainData.pairedOutput,
                   NULL, // Sample weight
                   0,    // Sample Index,
                   CvANN_MLP_TrainParams( cvTermCriteria( CV_TERMCRIT_ITER
                                                          | CV_TERMCRIT_EPS,
                                                          trainEpochIteration, trainDesiredError),
                                          CvANN_MLP_TrainParams::RPROP,  //
                                          0.1,                           //
                                          FLT_EPSILON ),
                   updateOption);
}

inline void
BRAINSCutTrainModel
::SaveRFTrainModelAtIteration( CvRTrees& myTrainer, int depth, int NTrees)
{
  std::string filename = GetRFModelFilename( depth,
                                             NTrees );

  try
    {
    std::cout << "Save Model File :: " << filename << std::endl;
    myTrainer.save( filename.c_str() );
    }
  catch( ... )
    {
    std::cout << "Fail to save the model file ::" << std::endl
              << filename << std::endl;
    exit(EXIT_FAILURE);
    }
}

inline void
BRAINSCutTrainModel
::SaveANNTrainModelAtIteration( neuralNetType& myTrainer, unsigned int No)
{
  char tempid[10];

  sprintf( tempid, "%09u", No + 1 );
  std::string filename = modelBasename + tempid;

  /** check the directory */
  std::string path = itksys::SystemTools::GetFilenamePath( filename );
  if( !itksys::SystemTools::FileExists( path.c_str(), false ) )
    {
    std::cout << " A directory for training file does not exist. Create as following:: "
              << filename.c_str()
              << std::endl;
    itksys::SystemTools::MakeDirectory( path.c_str() );
    }
  myTrainer.save( filename.c_str() );
}

inline void
BRAINSCutTrainModel
::writeRFTrainInformation( CvRTrees& myTrainer,
                           int depth,
                           int nTree)
{
  char* cError = new char[40];

  sprintf( cError, "%.5g", myTrainer.get_train_error() );

  char* cDepth = new char[10];
  sprintf( cDepth, "%d", depth );

  char* cNTree = new char[10];
  sprintf( cNTree, "%d", nTree);

  std::string line = "error,";
  line += cError;
  line += ", depth, ";
  line += cDepth;
  line += ", number of Tree, ";
  line += cNTree;

  std::cout << line << std::endl;
  appendToFile( RFErrorFilename, line );
}

inline void
BRAINSCutTrainModel
::appendToFile( std::string filename, std::string line )
{
  fstream filestr;

  filestr.open( filename.c_str(), std::ios::app | std::ios::out );
  if( !filestr.good() )
    {
    std::cout << "Cannot open the file :: " << filename << std::endl
              << "Fail to append line  :: " << line << std::endl;
    exit( EXIT_FAILURE );
    }
  filestr << line << std::endl;
  filestr.close();
}

inline void
BRAINSCutTrainModel
::printANNTrainInformation( neuralNetType& myTrainer, unsigned int No )
{
  std::cout << " Error, " << myTrainer.get_MSE()
            << " Iteration, " << No + 1
            << std::endl;
}

inline int *
BRAINSCutTrainModel
::GetANNLayerStructureArray()
{
  int * layer = new int[3];

  layer[0] = trainingDataSet->GetInputVectorSize();
  layer[1] = ANNHiddenNodesNumber;
  layer[2] = trainingDataSet->GetOutputVectorSize();

  return layer;
}

void
BRAINSCutTrainModel
::TrainANN()
{
  neuralNetType * trainner = new neuralNetType();
  int*            layer = GetANNLayerStructureArray();

  cvInitMatHeader( ANNLayerStructure, 1, 3, CV_32SC1, layer );

  trainner->create( ANNLayerStructure,
                    CvANN_MLP::SIGMOID_SYM,
                    activationSlope,
                    activationMinMax);
  for( unsigned int currentIteration = 0;
       currentIteration < trainIteration;
       currentIteration++ )
    {
    unsigned int subSetNo =  currentIteration % trainingDataSet->GetNumberOfSubSet();
    TrainWithUpdate( *trainner,
                     (currentIteration > 0),
                     *(trainingDataSet->GetTrainingSubSet(subSetNo) ) );
    SaveANNTrainModelAtIteration( *trainner, currentIteration );
    printANNTrainInformation( *trainner, currentIteration );
    if( trainner->get_MSE()  < trainDesiredError )
      {
      std::cout << "CONVERGED with " << trainner->get_MSE() << std::endl;
      break;
      }
    }
}

/** random forest training */
void
BRAINSCutTrainModel
::TrainRandomForest()
{
  for( int depth = 1; depth < trainMaxDepth; depth++ )
    {
    for( int nTree = 2; nTree < trainMaxTreeCount; nTree++ )
      {
      TrainRandomForestAt( depth, nTree );
      }
    }
}

void
BRAINSCutTrainModel
::TrainRandomForestAt( const int depth, const int numberOfTree )
{
  CvRTrees   forest;
  CvRTParams randomTreeTrainParamters =
    CvRTParams( depth,
                trainMinSampleCount,
                0.0F,                  // float  _regression_accuracy=0,
                trainUseSurrogates,
                10,                     // int    _max_categories=10,
                0,                      // float* _priors,
                trainCalcVarImportance, // bool   _calc_var_importance=false,
                0,                      // int    _nactive_vars=0,
                numberOfTree,
                0,                     // float  forest_accuracy=0,
                0
                );

  forest.train( trainingDataSet->GetTrainingSubSet(0)->pairedInput,
                CV_ROW_SAMPLE,   // or CV_COL_SAMPLE
                trainingDataSet->GetTrainingSubSet(0)->pairedOutputRF,
                0,
                0,    // CvMat* sampleIdx=0,
                0,    // CvMat* varType=0,
                0,    // CvMat* missingMask=0,
                randomTreeTrainParamters
                );
  writeRFTrainInformation( forest, depth, numberOfTree );
  SaveRFTrainModelAtIteration( forest, depth, numberOfTree);
}

/** setting function with net configuration */
void
BRAINSCutTrainModel
::SetModelBasename()
{
  NeuralParams * model = BRAINSCutNetConfiguration.Get<NeuralParams>("NeuralNetParams");

  modelBasename = model->GetAttribute<StringValue>("TrainingModelFilename");
}

void
BRAINSCutTrainModel
::SetRFErrorFilename()
{
  if( modelBasename == "" )
    {
    std::cout << "Model Basename has to be set first."
              << std::endl
              << __LINE__ << "::" << __FILE__
              << std::endl;
    exit( EXIT_FAILURE );
    }
  RFErrorFilename = modelBasename + "ERROR.txt";
}

void
BRAINSCutTrainModel
::SetRFErrorFile()
{
  fstream  filestr;
  ifstream ifile(RFErrorFilename.c_str() );

  if( ifile )
    {
    std::cout << "file already exist. Append to the file"
              << std::endl;
    ifile.close();
    return;
    }

  filestr.open( RFErrorFilename.c_str(), fstream::out);

  if( !filestr.good() )
    {
    std::cout << "Cannot open the file :: " << RFErrorFilename << std::endl
              << __LINE__ << "::" << __FILE__ << std::endl;
    exit( EXIT_FAILURE );
    }

  filestr << "E, error, D, depth, N, NTree\n";
  filestr.close();
}

std::string
BRAINSCutTrainModel
::GetANNVectorFilenamePrefix()
{
  NeuralParams * model = BRAINSCutNetConfiguration.Get<NeuralParams>("NeuralNetParams");

  return model->GetAttribute<StringValue>("TrainingVectorFilename");
}

void
BRAINSCutTrainModel
::SetIterationFromNetConfiguration()
{
  trainIteration = TrainNetConfiguration->GetAttribute<IntValue>("Iterations");
}

void
BRAINSCutTrainModel
::SetEpochIterationFromNetConfiguration()
{
  trainEpochIteration = TrainNetConfiguration->GetAttribute<IntValue>("EpochIterations");
}

void
BRAINSCutTrainModel
::SetDesiredErrorFromNetConfiguration()
{
  trainDesiredError = TrainNetConfiguration->GetAttribute<FloatValue>("DesiredError");
}

void
BRAINSCutTrainModel
::SetMaximumDataSizeFromNetConfiguration()
{
  trainMaximumDataSize = TrainNetConfiguration->GetAttribute<IntValue>("MaximumVectorsPerEpoch");
}

void
BRAINSCutTrainModel
::SetANNHiddenNodesNumberFromNetConfiguration()
{
  ANNHiddenNodesNumber = TrainNetConfiguration->GetAttribute<IntValue>("NumberOfHiddenNodes");
}

void
BRAINSCutTrainModel
::SetActivatioinFunctionFromNetConfiguration()
{
  activationSlope = TrainNetConfiguration->GetAttribute<FloatValue>("ActivationSlope");
  activationMinMax = TrainNetConfiguration->GetAttribute<FloatValue>("ActivationMinMax");
}

/** default functions to set/get member variables */
void
BRAINSCutTrainModel
::SetIteration(unsigned int iteration)
{
  trainIteration = iteration;
}

unsigned int
BRAINSCutTrainModel
::GetIteration()
{
  return trainIteration;
}

void
BRAINSCutTrainModel
::SetEpochIteration( unsigned int epochIteration)
{
  trainEpochIteration = epochIteration;
}

unsigned int
BRAINSCutTrainModel
::GetEpochIteration()
{
  return trainEpochIteration;
}

void
BRAINSCutTrainModel
::SetDesiredError( float desiredError )
{
  trainDesiredError = desiredError;
}

float
BRAINSCutTrainModel
::GetDesiredError()
{
  return trainDesiredError;
}

void
BRAINSCutTrainModel
::SetMaximumDataSize( unsigned int maximumDataSize)
{
  trainMaximumDataSize = maximumDataSize;
}

unsigned int
BRAINSCutTrainModel
::GetMaximumDataSize()
{
  return trainMaximumDataSize;
}

void
BRAINSCutTrainModel
::SetANNHiddenNodesNumber( int hiddenNodes)
{
  ANNHiddenNodesNumber = hiddenNodes;
}

int
BRAINSCutTrainModel
::GetANNHiddenNodesNumber()
{
  return ANNHiddenNodesNumber;
}

void
BRAINSCutTrainModel
::SetActivationFunction( float slope, float minMax)
{
  activationSlope = slope;
  activationMinMax = minMax;
}

float
BRAINSCutTrainModel
::GetActivationSlope()
{
  return activationSlope;
}

float
BRAINSCutTrainModel
::GetActivationMinMax()
{
  return activationMinMax;
}

/** Set Random Tree */
void
BRAINSCutTrainModel
::SetMaxDepthFromNetConfiguration()
{
  trainMaxDepth = TrainNetConfiguration->GetAttribute<IntValue>("MaxDepth");
}

void
BRAINSCutTrainModel
::SetMinSampleCountFromNetConfiguration()
{
  trainMinSampleCount = TrainNetConfiguration->GetAttribute<IntValue>("MinSampleCount");
}

void
BRAINSCutTrainModel
::SetUseSurrogatesFromNetConfiguration()
{
  trainUseSurrogates = TrainNetConfiguration->GetAttribute<BooleanValue>("UseSurrogates");
}

void
BRAINSCutTrainModel
::SetCalcVarImportanceFromNetConfiguration()
{
  trainCalcVarImportance = TrainNetConfiguration->GetAttribute<BooleanValue>("CalcVarImportance");
}

void
BRAINSCutTrainModel
::SetMaxTreeCountFromNetConfiguration()
{
  trainMaxTreeCount = TrainNetConfiguration->GetAttribute<IntValue>("MaxTreeCount");
}
