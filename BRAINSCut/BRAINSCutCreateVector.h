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
