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

#include "BRAINSCutCLP.h"

int main(int argc, char * *argv)
{
  PARSE_ARGS;

  /* Solution from Kent
   * ITK4 resigration initilization
   */
  // Call register default transforms
  itk::TransformFactoryBase::RegisterDefaultTransforms();

  if( generateProbability )
    {
    BRAINSCutGenerateProbability testBRAINSCutClass( netConfiguration );
    testBRAINSCutClass.GenerateProbabilityMaps();
    }
  if( createVectors )
    {
    BRAINSCutCreateVector testCreateVector( netConfiguration );
    testCreateVector.SetTrainingDataSetFromNetConfiguration();
    testCreateVector.CreateVectors();
    }
  if( trainModel )
    {
    try
      {
      BRAINSCutTrainModel testTrain( netConfiguration);
      testTrain.InitializeNeuralNetwork();
      testTrain.InitializeTrainDataSet();
      testTrain.Train();
      }
    catch( BRAINSCutExceptionStringHandler& e )
      {
      std::cout << e.Error();
      }
    }
  if( applyModel )
    {
    try
      {
      BRAINSCutApplyModel applyTest( netConfiguration );

      applyTest.SetRegionsOfInterestFromNetConfiguration();
      applyTest.SetRegistrationParametersFromNetConfiguration();
      applyTest.SetAtlasDataSet();
      applyTest.SetAtlasImage();
      applyTest.SetRhoPhiThetaFromNetConfiguration();

      applyTest.SetANNModelConfiguration();
      applyTest.SetANNTestingSSEFilename();
      applyTest.SetGradientSizeFromNetConfiguration();
      applyTest.SetANNOutputThresholdFromNetConfiguration();
      applyTest.SetTrainIterationFromNetConfiguration();
      applyTest.SetComputeSSE(false);

      applyTest.SetApplyDataSetFromNetConfiguration();
      applyTest.SetANNModelFilenameFromNetConfiguration();
      applyTest.SetComputeSSE( computeSSEOn );

      applyTest.Apply();
      }
    catch( BRAINSCutExceptionStringHandler& e )
      {
      std::cout << e.Error();
      }
    }

  return EXIT_SUCCESS;
}
