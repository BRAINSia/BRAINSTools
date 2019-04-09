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
#include "BRAINSCutTrainModel.h"
#include "TrainingVectorConfigurationType.h"
#include <fstream>

BRAINSCutTrainModel
::BRAINSCutTrainModel( BRAINSCutDataHandler& dataHandler ) :
  m_myDataHandler(dataHandler),
  m_trainIteration(0),
  m_trainEpochIteration(0),
  m_trainDesiredError(0.0),
  m_trainMaximumDataSize(0),
  m_ANNHiddenNodesNumber(0),
  m_activationSlope(0),
  m_activationMinMax(0),
  m_trainMaxDepth(0),
  m_trainMinSampleCount(0),
  m_trainUseSurrogates(false),
  m_trainCalcVarImportance(false),
  m_trainMaxTreeCount(0),
  m_modelBasename(""),
  m_ANNVectorFilenamePrefix(""),
  m_RFErrorFilename(""),
  m_trainingDataSet(nullptr),
  m_ANNLayerStructure(nullptr)
{
}

BRAINSCutTrainModel
::~BRAINSCutTrainModel()
{
  // Must release the memory so that it does not leak
  if( this->m_ANNLayerStructure != nullptr )
    {
    cvReleaseMat( &this->m_ANNLayerStructure );
    }
  if( this->m_trainingDataSet != nullptr )
    {
    delete this->m_trainingDataSet;
    }
}

/** train */
void
BRAINSCutTrainModel
::InitializeTrainDataSet( bool doShuffle )
{
  m_myDataHandler.SetTrainVectorFilename();
  const std::string trainVectorFilename = m_myDataHandler.GetTrainVectorFilename();

  this->m_trainingDataSet = new BRAINSCutVectorTrainingSet( trainVectorFilename );
  try
    {
    this->m_trainingDataSet->ReadHeaderFileInformation();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error();
    exit(EXIT_FAILURE);
    }
  this->m_trainingDataSet->SetRecordSize();
  this->m_trainingDataSet->SetBufferRecordSize();
  if( doShuffle )
    {
    this->m_trainingDataSet->RandomizeTrainingVector();
    }
  if( this->m_trainingDataSet->GetTotalVectorSize() > (int)m_trainMaximumDataSize )
    {
    unsigned int numberOfSubSet =
      (float)this->m_trainingDataSet->GetTotalVectorSize() / (float)m_trainMaximumDataSize;
    numberOfSubSet = ceil(numberOfSubSet) + 1;
    std::cout << "Divide subset into " << numberOfSubSet << std::endl;
    this->m_trainingDataSet->SetNumberOfSubSet( numberOfSubSet );
    }
  else
    {
    this->m_trainingDataSet->SetNumberOfSubSet( 1 );
    }
}

/** initialization for the model */
void
BRAINSCutTrainModel
::InitializeNeuralNetwork()
{
  m_myDataHandler.SetTrainConfiguration( "ANNParameters" );
  m_myDataHandler.SetTrainingVectorConfiguration();
  if( this->m_ANNLayerStructure != nullptr )
    {
    cvReleaseMat( &this->m_ANNLayerStructure );
    }
  this->m_ANNLayerStructure = cvCreateMat( 1, 3, CV_32SC1);

  m_trainMaximumDataSize = m_myDataHandler.GetMaximumDataSize();
  m_trainIteration       = m_myDataHandler.GetTrainIteration();
  m_trainEpochIteration  = m_myDataHandler.GetEpochIteration();
  m_trainDesiredError    = m_myDataHandler.GetDesiredError();

  m_ANNHiddenNodesNumber = m_myDataHandler.GetANNHiddenNodesNumber();

  m_activationSlope = m_myDataHandler.GetActivationFunctionSlope();
  m_activationMinMax = m_myDataHandler.GetActivationFunctionMinMax();
}

void
BRAINSCutTrainModel
::InitializeRandomForest()
{
  m_myDataHandler.SetTrainConfiguration( "RandomForestParameters" );
  m_myDataHandler.SetTrainingVectorConfiguration();
  m_trainMaxDepth          = m_myDataHandler.GetMaxDepth();
  m_trainMinSampleCount    = m_myDataHandler.GetMinSampleCount();
  m_trainUseSurrogates     = m_myDataHandler.GetUseSurrogates();
  m_trainCalcVarImportance = m_myDataHandler.GetCalcVarImportance();
  m_trainMaxTreeCount      = m_myDataHandler.GetMaxTreeCount();

  SetRFErrorFilename();
  SetRFErrorFile();
}

void
BRAINSCutTrainModel
::TrainWithUpdate( cv::Ptr<OpenCVMLPType> myTrainer, pairedTrainingSetType& currentTrainData )
{
  // TODO change subset number properly
  myTrainer->setBackpropMomentumScale(0.1);
  myTrainer->setTrainMethod(OpenCVMLPType::RPROP);
  myTrainer->setTermCriteria(cvTermCriteria( CV_TERMCRIT_ITER
                           | CV_TERMCRIT_EPS,
                           m_trainEpochIteration, m_trainDesiredError));
  myTrainer->train( cv::ml::TrainData::create(currentTrainData.pairedInput,
                                              cv::ml::ROW_SAMPLE,
                                              currentTrainData.pairedOutput));
}

void
BRAINSCutTrainModel
::SaveRFTrainModelAtIteration( cv::Ptr<cv::ml::RTrees> myTrainer, int depth, int NTrees)
{
  std::string filename = m_myDataHandler.GetRFModelFilename( depth, NTrees );

  try
    {
    std::cout << "Save Model File :: " << filename << std::endl;
    myTrainer->save( filename.c_str() );
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
::SaveANNTrainModelAtIteration( cv::Ptr<OpenCVMLPType>  myTrainer, unsigned int No)
{
  char tempid[10];

  sprintf( tempid, "%09u", No );
  std::string filename = m_myDataHandler.GetModelBaseName() + tempid;

  /** check the directory */
  std::string path = itksys::SystemTools::GetFilenamePath( filename );
  if( !itksys::SystemTools::FileExists( path.c_str(), false ) )
    {
    std::cout << " A directory for training file does not exist. Create as following:: "
              << filename.c_str()
              << std::endl;
    itksys::SystemTools::MakeDirectory( path.c_str() );
    }
  myTrainer->save( filename.c_str() );
}

void
BRAINSCutTrainModel
::recordRFTrainInformation( int depth,
                           int nTree)
{
  char cDepth[10];
  sprintf( cDepth, "%d", depth );

  char cNTree[10];
  sprintf( cNTree, "%d", nTree);

  std::string line = "error,";
  line += ", depth, ";
  line += cDepth;
  line += ", number of Tree, ";
  line += cNTree;

  appendToFile( m_RFErrorFilename, line );
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
::printANNTrainInformation( cv::Ptr<OpenCVMLPType> /*NOT USED myTrainer */, unsigned int No )
{
  std::cout << " Error, " << " NOT_COMPUTED " /* myTrainer.get_MSE() This needs to be instumented again */
            << " Iteration, " << No
            << std::endl;
}

void
BRAINSCutTrainModel
::FillANNLayerStructureArray3D( int * const layer ) const
{
  layer[0] = this->m_trainingDataSet->GetInputVectorSize();
  layer[1] = m_ANNHiddenNodesNumber;
  layer[2] = this->m_trainingDataSet->GetOutputVectorSize();
}

void
BRAINSCutTrainModel
::TrainANN()
{
  cv::Ptr<OpenCVMLPType>  trainner =  OpenCVMLPType::create();
  int             layer[3];

  this->FillANNLayerStructureArray3D(layer);
  constexpr int nLayer = 3;
  //cv::Mat annLayer = cv::Mat(nLayer, layer, CV_32SC1, cv::Scalar::all(0.0F));
  cv::Mat annLayer = cv::Mat(1, nLayer, CV_32SC1, layer);
  trainner->setActivationFunction(OpenCVMLPType::SIGMOID_SYM);
  std::cout<<"ANNLayer::"
           <<annLayer
           <<std::endl;
  trainner->setLayerSizes(annLayer);
  std::cout<<"getLayerSizes::"
           <<trainner->getLayerSizes()
           <<std::endl;
  for( unsigned int currentIteration = 1;
       currentIteration <= m_trainIteration;
       currentIteration++ )
    {
    unsigned int subSetNo =  (currentIteration - 1) % this->m_trainingDataSet->GetNumberOfSubSet();
    TrainWithUpdate( trainner,
                     *(this->m_trainingDataSet->GetTrainingSubSet(subSetNo) ) );
    SaveANNTrainModelAtIteration( trainner, currentIteration );
    printANNTrainInformation( trainner, currentIteration );
    /*
    if( trainner->get_MSE()  < m_trainDesiredError )
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
  for( int depth = 1; depth <= m_trainMaxDepth; depth++ )
    {
    for( int nTree = 2; nTree <= m_trainMaxTreeCount; nTree++ )
      {
      TrainRandomForestAt( depth, nTree );
      }
    }
}

void
BRAINSCutTrainModel
::TrainRandomForestAt( const int depth, const int numberOfTree )
{
  cv::Ptr<cv::ml::RTrees> forest = cv::ml::RTrees::create();
  forest->setMaxDepth(depth);
  forest->setMinSampleCount(m_trainMinSampleCount);
  forest->setUseSurrogates(m_trainUseSurrogates);
  forest->setCalculateVarImportance(m_trainCalcVarImportance);
  forest->setTermCriteria( cv::TermCriteria(
                              cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS,
                              numberOfTree,
                              0.1F));
  forest->train(cv::ml::TrainData::create( this->m_trainingDataSet->GetTrainingSubSet(0)->pairedInput,
                                          cv::ml::ROW_SAMPLE,
                                           this->m_trainingDataSet->GetTrainingSubSet(0)->pairedOutputRF));

  recordRFTrainInformation( depth, numberOfTree );
  SaveRFTrainModelAtIteration( forest, depth, numberOfTree);
}

void
BRAINSCutTrainModel
::SetRFErrorFilename()
{
  const std::string basename = m_myDataHandler.GetModelBaseName();

  if( basename.empty() )
    {
    std::cout << "Model Basename has to be set first."
              << std::endl
              << __LINE__ << "::" << __FILE__
              << std::endl;
    exit( EXIT_FAILURE );
    }
  m_RFErrorFilename = basename + "ERROR.txt";
}

void
BRAINSCutTrainModel
::SetRFErrorFile()
{
  std::fstream  filestr;
  std::ifstream ifile(m_RFErrorFilename.c_str() );

  if( ifile )
    {
    std::cout << "file already exist. Append to the file"
              << std::endl;
    ifile.close();
    return;
    }

  filestr.open( m_RFErrorFilename.c_str(), std::fstream::out);

  if( !filestr.good() )
    {
    std::cout << "Cannot open the file :: " << m_RFErrorFilename << std::endl
              << __LINE__ << "::" << __FILE__ << std::endl;
    exit( EXIT_FAILURE );
    }

  filestr << "E, error, D, depth, N, NTree\n";
  filestr.close();
}
