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
  itk::TransformFactoryBase::RegisterDefaultTransforms();

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  BRAINSCutGenerateRegistrations registrationGenerator( netConfiguration );
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  const bool applyDataSetOff = false;
  const bool applyDataSetOn = true;

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
      std::cout << __LINE__ << "::" << __FILE__ << std::endl;
      registrationGenerator.SetAtlasToSubjectRegistrationOn( true );
      registrationGenerator.SetSubjectDataSet( applyDataSetOn );
      registrationGenerator.GenerateRegistrations();

      BRAINSCutApplyModel applyTest( netConfiguration );

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
