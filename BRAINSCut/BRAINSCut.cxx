#include <string>
#include <iostream>
#include "GenericTransformImage.h"
#include "itkScaleVersor3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkTransformFactory.h"

#include "BRAINSThreadControl.h"
#include "BRAINSCutGenerateProbability.h"
#include "BRAINSCutCreateVector.h"
#include "BRAINSCutVectorTrainingSet.h"
#include "BRAINSCutTrainModel.h"
#include "BRAINSCutApplyModel.h"
#include "BRAINSCutGenerateRegistrations.h"

#include "BRAINSCutCLP.h"

int main(int argc, char * *argv)
{
  PARSE_ARGS;

  /* Solution from Kent
   * ITK4 resigration initilization
   */
  // Call register default transforms
  // itk::TransformFactoryBase::RegisterDefaultTransforms();

  if( !netConfiguration.empty() && modelConfigurationFilename.empty() )
    {
    modelConfigurationFilename = netConfiguration;
    }

  // Data handler
  if( !itksys::SystemTools::FileExists( modelConfigurationFilename.c_str() ) )
    {
    std::string errorMsg = " File does not exist! :";
    errorMsg += modelConfigurationFilename;
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }
  BRAINSCutDataHandler dataHandler( modelConfigurationFilename );

  BRAINSCutGenerateRegistrations registrationGenerator( dataHandler );
  const bool                     applyDataSetOff = false;
  const bool                     applyDataSetOn = true;
  const bool                     shuffleTrainVector = (NoTrainingVectorShuffling != true );

  std::cout << "shuffleTrainVector::" << shuffleTrainVector << std::endl;

  if( generateProbability )
    {
    registrationGenerator.SetAtlasToSubjectRegistrationOn( false );
    registrationGenerator.SetDataSet( applyDataSetOff );
    registrationGenerator.GenerateRegistrations();

    BRAINSCutGenerateProbability testBRAINSCutClass( dataHandler );
    testBRAINSCutClass.GenerateProbabilityMaps();
    }
  if( createVectors )
    {
    registrationGenerator.SetAtlasToSubjectRegistrationOn( true );
    registrationGenerator.SetDataSet( applyDataSetOff );
    registrationGenerator.GenerateRegistrations();

    BRAINSCutCreateVector testCreateVector( dataHandler);
    testCreateVector.SetTrainingDataSet();
    testCreateVector.CreateVectors();
    }
  if( trainModel )
    {
    if( method == "ANN" )
      {
      try
        {
        BRAINSCutTrainModel ANNTrain( dataHandler );
        ANNTrain.InitializeNeuralNetwork();
        ANNTrain.InitializeTrainDataSet( shuffleTrainVector );
        ANNTrain.TrainANN();
        }
      catch( BRAINSCutExceptionStringHandler& e )
        {
        std::cout << e.Error();
        }
      }
    else if( method == "RandomForest" )
      {
      BRAINSCutTrainModel RandomForestTrain( dataHandler );
      RandomForestTrain.InitializeRandomForest();
      RandomForestTrain.InitializeTrainDataSet( shuffleTrainVector);

      // these set has to be **AFTER** InitializeTrainDataSet
      if( numberOfTrees > 0 && randomTreeDepth > 0 )
        {
        RandomForestTrain.TrainRandomForestAt( randomTreeDepth, numberOfTrees );
        }
      else
        {
        RandomForestTrain.TrainRandomForest();
        }
      }
    else
      {
      std::cout << "No proper method found to train"
                << std::endl;
      exit(EXIT_FAILURE);
      }
    }
  if( applyModel )
    {
    try
      {
      registrationGenerator.SetAtlasToSubjectRegistrationOn( true );
      registrationGenerator.SetDataSet( applyDataSetOn );
      registrationGenerator.GenerateRegistrations();

      dataHandler.SetRandomForestModelFilename( modelFilename );

      BRAINSCutApplyModel ApplyModule( dataHandler );

      ApplyModule.SetMethod( method );
      ApplyModule.SetComputeSSE( computeSSEOn );

      if( method == "RandomForest" )
        {
        ApplyModule.SetDepthOfTree( randomTreeDepth );
        ApplyModule.SetNumberOfTrees( numberOfTrees );
        }
      ApplyModule.Apply();
      }
    catch( BRAINSCutExceptionStringHandler& e )
      {
      std::cout << e.Error();
      }
    }

  return EXIT_SUCCESS;
}
