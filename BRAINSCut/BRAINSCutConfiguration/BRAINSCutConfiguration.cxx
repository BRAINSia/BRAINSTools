#include <DataSet.h>
#include <RegistrationConfigurationParser.h>
#include <ProbabilityMapParser.h>

#include <SpatialLocationType.h>

#include <list>
#include <BRAINSCutConfiguration.h>

BRAINSCutConfiguration::BRAINSCutConfiguration() : ElementParser("BRAINSCutBRAINSCutConfiguration")
{
  this->Add(new DataSetList, "DataSetList");
  this->Add(new ProbabilityMapList, "ProbabilityMapList");
  this->Add(new RegistrationConfigurationParser, "RegistrationConfiguration");
}

void
BRAINSCutConfiguration::AddDataSet(DataSet *newSet)
{
  DataSetList *set = this->Get<DataSetList>("DataSetList");

  set->Add( newSet,
            newSet->GetAttribute<StringValue>("Name") );
}

DataSet * BRAINSCutConfiguration::GetAtlasDataSet() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  for( ElementParser::const_iterator it = set->begin();
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

BRAINSCutConfiguration::TrainDataSetListType
BRAINSCutConfiguration::GetTrainDataSets() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  std::list<DataSet *> rval;

  for( ElementParser::const_iterator it = set->begin();
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

BRAINSCutConfiguration::ApplyDataSetListType
BRAINSCutConfiguration::GetApplyDataSets() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  std::list<DataSet *> rval;

  for( ElementParser::const_iterator it = set->begin();
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
BRAINSCutConfiguration::GetImageTypes() const
{
  const DataSetList *set = this->Get<DataSetList>("DataSetList");

  if( set->size() == 0 )
    {
    return DataSet::StringVectorType();
    }
  return dynamic_cast<const DataSet *>( set->begin()->second )->GetImageTypes();
}

// Set/Get Functions
//
void
BRAINSCutConfiguration::SetImageTypeToUse( std::string imageTypeToUse )
{
  ImageTypeToUse = imageTypeToUse;
}

std::string
BRAINSCutConfiguration::GetImageTypeToUse()
{
  return ImageTypeToUse;
}

void
BRAINSCutConfiguration::SetRegistrationID( std::string registrationID )
{
  RegistrationID = registrationID;
}

std::string
BRAINSCutConfiguration::GetRegistrationID()
{
  return RegistrationID;
}
