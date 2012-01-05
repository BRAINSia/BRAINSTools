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

  BRAINSCutGenerateRegistrations registrationGenerator( netConfiguration );
  const bool                     applyDataSetOff = false;
  const bool                     applyDataSetOn = true;
  const bool                     shuffleTrainVector = (NoTrainingVectorShuffling != true );

  std::cout << "shuffleTrainVector::" << shuffleTrainVector << std::endl;

  if( generateProbability )
    {
    registrationGenerator.SetAtlasToSubjectRegistrationOn( false );
    registrationGenerator.SetSubjectDataSet( applyDataSetOff );
    registrationGenerator.GenerateRegistrations();

    BRAINSCutGenerateProbability testBRAINSCutClass( netConfiguration );
    testBRAINSCutClass.GenerateProbabilityMaps();
    }
  if( createVectors )
    {
    registrationGenerator.SetAtlasToSubjectRegistrationOn( true );
    registrationGenerator.SetSubjectDataSet( applyDataSetOff );
    registrationGenerator.GenerateRegistrations();

    BRAINSCutCreateVector testCreateVector( netConfiguration );
    testCreateVector.SetTrainingDataSetFromNetConfiguration();
    testCreateVector.CreateVectors();
    }
  if( trainModel )
    {
    if( method == "ANN" )
      {
      try
        {
        BRAINSCutTrainModel ANNTrain( netConfiguration);
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
      BRAINSCutTrainModel RandomForestTrain( netConfiguration );
      RandomForestTrain.InitializeRandomForest();
      RandomForestTrain.InitializeTrainDataSet( shuffleTrainVector);

      /* these set has to be **AFTER** InitializeTrainDataSet */
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
      registrationGenerator.SetSubjectDataSet( applyDataSetOn );
      registrationGenerator.GenerateRegistrations();

      BRAINSCutApplyModel ApplyModule( netConfiguration );

      ApplyModule.SetMethod( method );
      ApplyModule.SetComputeSSE( computeSSEOn );
      /* these set has to be **AFTER** InitializeTrainDataSet */
      if( numberOfTrees > 0 && randomTreeDepth > 0 )
        {
        ApplyModule.SetRandomForestModelFilename( randomTreeDepth, numberOfTrees );
        }
      else
        {
        ApplyModule.SetRandomForestModelFilenameFromNetConfiguration();
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
