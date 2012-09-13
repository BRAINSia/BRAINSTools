#include "BRAINSCutTrainModel.h"
#include "TrainingVectorConfigurationType.h"
#include <fstream>

BRAINSCutTrainModel
::BRAINSCutTrainModel( BRAINSCutDataHandler & dataHandler ) :
  trainIteration(0),
  trainEpochIteration(0),
  trainDesiredError(0.0),
  trainMaximumDataSize(0),
  ANNHiddenNodesNumber(0),
  activationSlope(0),
  activationMinMax(0),
  m_ANNLayerStructure(NULL)
{
  myDataHandler = dataHandler;
}

BRAINSCutTrainModel
::~BRAINSCutTrainModel()
{
  // Must release the memory so that it does not leak
  if( this->m_ANNLayerStructure != NULL )
    {
    cvReleaseMat( &this->m_ANNLayerStructure );
    }
}

/** train */
void
BRAINSCutTrainModel
::InitializeTrainDataSet( bool doShuffle )
{
  myDataHandler.SetTrainVectorFilename();
  const std::string trainVectorFilename = myDataHandler.GetTrainVectorFilename();

  trainingDataSet = new BRAINSCutVectorTrainingSet( trainVectorFilename );
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
    trainingDataSet->RandomizeTrainingVector();
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
  myDataHandler.SetTrainConfiguration( "ANNParameters" );
  myDataHandler.SetTrainingVectorConfiguration();
  if( this->m_ANNLayerStructure != NULL )
    {
    cvReleaseMat( &this->m_ANNLayerStructure );
    }
  this->m_ANNLayerStructure = cvCreateMat( 1, 3, CV_32SC1);

  trainMaximumDataSize = myDataHandler.GetMaximumDataSize();
  trainIteration       = myDataHandler.GetTrainIteration();
  trainEpochIteration  = myDataHandler.GetEpochIteration();
  trainDesiredError    = myDataHandler.GetDesiredError();

  ANNHiddenNodesNumber = myDataHandler.GetANNHiddenNodesNumber();

  activationSlope = myDataHandler.GetActivationFunctionSlope();
  activationMinMax = myDataHandler.GetActivationFunctionMinMax();
}

void
BRAINSCutTrainModel
::InitializeRandomForest()
{
  myDataHandler.SetTrainConfiguration( "RandomForestParameters" );
  myDataHandler.SetTrainingVectorConfiguration();
  trainMaxDepth          = myDataHandler.GetMaxDepth();
  trainMinSampleCount    = myDataHandler.GetMinSampleCount();
  trainUseSurrogates     = myDataHandler.GetUseSurrogates();
  trainCalcVarImportance = myDataHandler.GetCalcVarImportance();
  trainMaxTreeCount      = myDataHandler.GetMaxTreeCount();

  SetRFErrorFilename();
  SetRFErrorFile();
}

void
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

void
BRAINSCutTrainModel
::SaveRFTrainModelAtIteration( CvRTrees& myTrainer, int depth, int NTrees)
{
  std::string filename = myDataHandler.GetRFModelFilename( depth, NTrees );

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

void
BRAINSCutTrainModel
::SaveANNTrainModelAtIteration( neuralNetType& myTrainer, unsigned int No)
{
  char tempid[10];

  sprintf( tempid, "%09u", No );
  std::string filename = myDataHandler.GetModelBaseName() + tempid;

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

void
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

  appendToFile( RFErrorFilename, line );
}

void
BRAINSCutTrainModel
::appendToFile( std::string filename, std::string line )
{
  std::fstream filestr;

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

void
BRAINSCutTrainModel
::printANNTrainInformation( neuralNetType & /*NOT USED myTrainer */, unsigned int No )
{
  std::cout << " Error, " << " NOT_COMPUTED " /* myTrainer.get_MSE() This needs to be instumented again */
            << " Iteration, " << No
            << std::endl;
}

int *
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

  cvInitMatHeader( this->m_ANNLayerStructure, 1, 3, CV_32SC1, layer );

  trainner->create( this->m_ANNLayerStructure,
                    CvANN_MLP::SIGMOID_SYM,
                    activationSlope,
                    activationMinMax);
  for( unsigned int currentIteration = 1;
       currentIteration <= trainIteration;
       currentIteration++ )
    {
    unsigned int subSetNo =  (currentIteration - 1) % trainingDataSet->GetNumberOfSubSet();
    TrainWithUpdate( *trainner,
                     (currentIteration > 1), // after first iteration, update
                     *(trainingDataSet->GetTrainingSubSet(subSetNo) ) );
    SaveANNTrainModelAtIteration( *trainner, currentIteration );
    printANNTrainInformation( *trainner, currentIteration );
    /*
    if( trainner->get_MSE()  < trainDesiredError )
      {
      std::cout << "CONVERGED with " << trainner->get_MSE() << std::endl;
      break;
      }
      */
    }
}

/** random forest training */
void
BRAINSCutTrainModel
::TrainRandomForest()
{
  for( int depth = 1; depth <= trainMaxDepth; depth++ )
    {
    for( int nTree = 2; nTree <= trainMaxTreeCount; nTree++ )
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

void
BRAINSCutTrainModel
::SetRFErrorFilename()
{
  const std::string basename = myDataHandler.GetModelBaseName();

  if( basename.empty() )
    {
    std::cout << "Model Basename has to be set first."
              << std::endl
              << __LINE__ << "::" << __FILE__
              << std::endl;
    exit( EXIT_FAILURE );
    }
  RFErrorFilename = basename + "ERROR.txt";
}

void
BRAINSCutTrainModel
::SetRFErrorFile()
{
  std::fstream  filestr;
  std::ifstream ifile(RFErrorFilename.c_str() );

  if( ifile )
    {
    std::cout << "file already exist. Append to the file"
              << std::endl;
    ifile.close();
    return;
    }

  filestr.open( RFErrorFilename.c_str(), std::fstream::out);

  if( !filestr.good() )
    {
    std::cout << "Cannot open the file :: " << RFErrorFilename << std::endl
              << __LINE__ << "::" << __FILE__ << std::endl;
    exit( EXIT_FAILURE );
    }

  filestr << "E, error, D, depth, N, NTree\n";
  filestr.close();
}
