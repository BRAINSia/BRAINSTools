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

  /** Test of Class
   */
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

      applyTest.SetApplyDataSetFromNetConfiguration();
      applyTest.SetANNModelFilenameFromNetConfiguration();

      applyTest.Apply();
      }
    catch( BRAINSCutExceptionStringHandler& e )
      {
      std::cout << e.Error();
      }
    }

/*
  BRAINSUtils::SetThreadCount(numberOfThreads);

  // Apparently when you register one transform, you need to register all your
  // transforms.
  itk::AddExtraTransformRegister();

  int status = 0;
  if( netConfiguration == "" )
    {
    std::cerr << "No XML Configuration file given" << std::endl;
    std::cerr.flush();
    }
  if( generateProbability )
    {
    extern int GenerateProbability(const std::string & xmlFile,
                                   int verbose,
                                   bool validate);
    status = GenerateProbability(netConfiguration,
                                 verbose,
                                 validate);
    if( verbose > 9 )
      {
      std::cout << " Status Returned from GenerateProbability is "
                << status
                << std::endl;
      }
    }
  if( status == 0 && createVectors )
    {
    extern int CreateVectors(const std::string & xmlFile,
                             bool histogramEqualization,
                             int verbose);
    status = CreateVectors(netConfiguration,
                           histogramEqualization,
                           verbose);
    if( verbose > 9 )
      {
      std::cout << " Status Returned from CreateVectors is "
                << status
                << std::endl;
      }
    }
  if( status == 0 && trainModel )
    {
    extern int TrainModel(const std::string & xmlFile,
                          int verbose,
                          int start_iteration
                          );
    status = TrainModel(netConfiguration,
                        verbose,
                        trainModelStartIndex);
    if( verbose > 9 )
      {
      std::cout << " Status Returned from CreateVectors is "
                << status
                << std::endl;
      }
    }
  if( status == 0 && applyModel )
    {
    extern int ApplyModel(const std::string & xmlFile,
                          int verbose,
                          bool histogramEqualization,
                          bool multiStructureThreshold);
    status = ApplyModel(netConfiguration,
                        verbose,
                        histogramEqualization,
                        multiStructureThreshold
                        );
    }
  exit(status);
  */
}
