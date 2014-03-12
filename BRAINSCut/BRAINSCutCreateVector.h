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
#ifndef BRAINSCutCreateVectorModel_h
#define BRAINSCutCreateVectorModel_h

#include "BRAINSCutDataHandler.h"
#include "FeatureInputVector.h"

class BRAINSCutCreateVector
{
public:
  BRAINSCutCreateVector( BRAINSCutDataHandler dataHandler );

  void SetTrainingDataSet();

  void SetTrainingVectorFilename();

  void CreateVectors();

  int  CreateSubjectVectors( DataSet& subject, std::ofstream& outputStream);

  void WriteCurrentVectors( InputVectorMapType& pairedInput, OutputVectorMapType& pairedOutput,
                            std::ofstream& outputStream );

  void WriteHeaderFile( std::string vectorFilename, int m_inputVectorSize, int m_outputVectorSize,
                        int numberOfInputVector);

private:
  int m_inputVectorSize;
  int m_outputVectorSize;
  BRAINSCutDataHandler                         m_myDataHandler;
  BRAINSCutConfiguration::TrainDataSetListType m_trainDataSetList;

  OutputVectorMapType GetPairedOutput( std::map<std::string,
                                                WorkingImagePointer>& deformedROIs, std::string roiName,
                                       std::string subjectROIBinaryFilename, int roiNumber);

  inline std::string GetROIBinaryFilename( DataSet& subject, std::string roiName);

  inline scalarType GetBinaryValue( scalarType value);
};
#endif
