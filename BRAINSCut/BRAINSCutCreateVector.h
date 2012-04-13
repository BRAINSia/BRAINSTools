#ifndef BRAINSCutCreateVectorModel_h
#define BRAINSCutCreateVectorModel_h

#include "BRAINSCutPrimary.h"
#include "FeatureInputVector.h"

class BRAINSCutCreateVector : public BRAINSCutPrimary
{
public:
  BRAINSCutCreateVector(std::string netConfigurationFilename);

  void SetTrainingDataSetFromNetConfiguration();

  void SetTrainingVectorFilenameFromNetConfiguration();

  void SetNormalizationFromNetConfiguration();

  void CreateVectors();

  int  CreateSubjectVectors( DataSet& subject, std::ofstream& outputStream);

  void WriteCurrentVectors( InputVectorMapType& pairedInput, OutputVectorMapType& pairedOutput,
                            std::ofstream& outputStream );

  void WriteHeaderFile( int inputVectorSize, int outputVectorSize, int numberOfInputVector);

private:
  BRAINSCutConfiguration::TrainDataSetListType trainDataSetList;
  bool                                         normalization;

  std::string vectorFilename;

  int inputVectorSize;
  int outputVectorSize;
  OutputVectorMapType GetPairedOutput( std::map<std::string,
                                                WorkingImagePointer>& deformedROIs, std::string roiName,
                                       std::string subjectROIBinaryFilename, int roiNumber);

  inline std::string GetROIBinaryFilename( DataSet& subject, std::string roiName);

  inline scalarType GetBinaryValue( scalarType value);
};
#endif
