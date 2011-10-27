#include <XMLElementParser.h>
#include <DataSet.h>
#include <RegistrationConfigurationParser.h>
#include <ProbabilityMapParser.h>
#include <list>
#include <SpatialLocationType.h>

#include <NetConfiguration.h>

NetConfiguration::NetConfiguration() : XMLElementParser("BRAINSCutNetConfiguration")
{
  this->Add(new DataSetList, "DataSetList");
  this->Add(new ProbabilityMapList, "ProbabilityMapList");
  this->Add(new RegistrationConfigurationParser, "RegistrationConfiguration");
}

void
NetConfiguration::AddDataSet(DataSet *newSet)
{
  DataSetList *set = this->Get<DataSetList>("DataSetList");

  set->Add( newSet,
            newSet->GetAttribute<StringValue>("Name") );
}

DataSet * NetConfiguration::GetAtlasDataSet() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  for( XMLElementParser::const_iterator it = set->begin();
       it != set->end();
       ++it )
    {
    DataSet *current = dynamic_cast<DataSet *>( it->second );
    if( current->GetAttribute<StringValue>("Type") == "Atlas" )
      {
      return current;
      }
    }
  return 0;
}

NetConfiguration::TrainDataSetListType
NetConfiguration::GetTrainDataSets() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  std::list<DataSet *> rval;

  for( XMLElementParser::const_iterator it = set->begin();
       it != set->end();
       ++it )
    {
    DataSet *   current = dynamic_cast<DataSet *>( it->second );
    std::string type( current->GetAttribute<StringValue>("Type") );
    if( type != "Atlas" && type != "Apply" )
      {
      rval.push_back(current);
      }
    }
  return rval;
}

NetConfiguration::ApplyDataSetListType
NetConfiguration::GetApplyDataSets() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  std::list<DataSet *> rval;

  for( XMLElementParser::const_iterator it = set->begin();
       it != set->end();
       ++it )
    {
    DataSet *current = dynamic_cast<DataSet *>( it->second );
    if( current->GetAttribute<StringValue>("Type") == "Apply" )
      {
      rval.push_back(current);
      }
    }
  return rval;
}

const DataSet::StringVectorType
NetConfiguration::GetImageTypes() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  if( set->size() == 0 )
    {
    return DataSet::StringVectorType();
    }
  return dynamic_cast<const DataSet *>( set->begin()->second )->GetImageTypes();
}
