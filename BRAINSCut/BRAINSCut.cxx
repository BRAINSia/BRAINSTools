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
  BRAINSRegisterAlternateIO();

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
  BRAINSCutDataHandler m_dataHandler( modelConfigurationFilename );

  BRAINSCutGenerateRegistrations m_registrationGenerator( m_dataHandler );
  const bool                     m_applyDataSetOff = false;
  const bool                     m_shuffleTrainVector = (NoTrainingVectorShuffling != true );

  std::cout << "m_shuffleTrainVector::" << m_shuffleTrainVector << std::endl;

  if( generateProbability )
    {
    m_registrationGenerator.SetAtlasToSubjectRegistrationOn( false );
    m_registrationGenerator.SetDataSet( m_applyDataSetOff );
    m_registrationGenerator.GenerateRegistrations();

    BRAINSCutGenerateProbability m_probabilityMapGenerator( m_dataHandler );
    m_probabilityMapGenerator.GenerateProbabilityMaps();
    }
  if( createVectors )
    {
    m_registrationGenerator.SetAtlasToSubjectRegistrationOn( true );
    m_registrationGenerator.SetDataSet( m_applyDataSetOff );
    m_registrationGenerator.GenerateRegistrations();

    BRAINSCutCreateVector m_trainingVectorCreator( m_dataHandler);
    m_trainingVectorCreator.SetTrainingDataSet();
    m_trainingVectorCreator.CreateVectors();
    }
  if( trainModel )
    {
    if( method == "ANN" )
      {
      try
        {
        BRAINSCutTrainModel m_ANNTrainer( m_dataHandler );
        m_ANNTrainer.InitializeNeuralNetwork();
        m_ANNTrainer.InitializeTrainDataSet( m_shuffleTrainVector );
        m_ANNTrainer.TrainANN();
        }
      catch( BRAINSCutExceptionStringHandler& e )
        {
        std::cout << e.Error();
        }
      }
    else if( method == "RandomForest" )
      {
      BRAINSCutTrainModel m_RandomForestTrainer( m_dataHandler );
      m_RandomForestTrainer.InitializeRandomForest();
      m_RandomForestTrainer.InitializeTrainDataSet( m_shuffleTrainVector);

      // these set has to be **AFTER** InitializeTrainDataSet
      if( numberOfTrees > 0 && randomTreeDepth > 0 )
        {
        m_RandomForestTrainer.TrainRandomForestAt( randomTreeDepth, numberOfTrees );
        }
      else
        {
        m_RandomForestTrainer.TrainRandomForest();
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
      const bool m_applyDataSetOn = true;
      m_registrationGenerator.SetAtlasToSubjectRegistrationOn( true );
      m_registrationGenerator.SetDataSet( m_applyDataSetOn );
      m_registrationGenerator.GenerateRegistrations();

      m_dataHandler.SetRandomForestModelFilename( modelFilename );

      BRAINSCutApplyModel m_ModelApplier( m_dataHandler );

      m_ModelApplier.SetMethod( method );
      m_ModelApplier.SetComputeSSE( computeSSEOn );

      if( method == "RandomForest" )
        {
        m_ModelApplier.SetDepthOfTree( randomTreeDepth );
        m_ModelApplier.SetNumberOfTrees( numberOfTrees );
        }
      m_ModelApplier.Apply();
      }
    catch( BRAINSCutExceptionStringHandler& e )
      {
      std::cout << e.Error();
      }
    }

  return EXIT_SUCCESS;
}
