#ifndef NetConfiguration_h
#define NetConfiguration_h
#include <XMLElementParser.h>
#include <DataSet.h>
#include <RegistrationConfigurationParser.h>
#include <ProbabilityMapParser.h>
#include <list>
#include <SpatialLocationType.h>

class NetConfiguration : public XMLElementParser
{
public:
  typedef XMLElementParser SuperClass;
  NetConfiguration();

  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    return indent + 2;
  }

  typedef std::list<DataSet *> TrainDataSetListType;
  typedef std::list<DataSet *> ApplyDataSetListType;

  void AddDataSet(DataSet *newSet);

  DataSet * GetAtlasDataSet() const;

  TrainDataSetListType GetTrainDataSets() const;

  ApplyDataSetListType GetApplyDataSets() const;

  const DataSet::StringVectorType GetImageTypes() const;

private:
  std::string ImageTypeToUse;
  std::string RegistrationID;
};

#endif // NetConfiguration_h
