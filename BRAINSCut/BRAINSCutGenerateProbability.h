#ifndef BRAINSCutGenerateProbability_h
#define BRAINSCutGenerateProbability_h

#include "BRAINSCutPrimary.h"
#include "NetConfiguration.h"
#include <itkIO.h>

class BRAINSCutGenerateProbability : private BRAINSCutPrimary
{
public:
  BRAINSCutGenerateProbability(std::string netConfigurationFilename);

  NetConfiguration * GetNetConfiguration();

  void SetNetConfiguration( NetConfiguration * netConfiguration);

  void SetNetConfigurationFilename(std::string filename);

  void SetNetConfiguration();

  void SetTrainingDataSetsList();

  std::string GetNetConfigurationFilename();

  void GenerateProbabilityMaps();

  WorkingImagePointer SmoothImage( const WorkingImagePointer image, const float GaussianValue);

private:

  /** DataSets */
  std::list<DataSet *> trainingDataSetList;

  void GenerateSymmetricalSphericalCoordinateImage();
};
#endif
