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
#include <DataSet.h>
#include <RegistrationConfigurationParser.h>
#include <ProbabilityMapParser.h>

#include <SpatialLocationType.h>

#include <list>
#include <BRAINSCutConfiguration.h>

BRAINSCutConfiguration::BRAINSCutConfiguration()
  : ElementParser("BRAINSCutBRAINSCutConfiguration")
{
  this->Add(new DataSetList, "DataSetList");
  this->Add(new ProbabilityMapList, "ProbabilityMapList");
  this->Add(new RegistrationConfigurationParser, "RegistrationConfiguration");
}

void
BRAINSCutConfiguration::AddDataSet(DataSet * newSet)
{
  DataSetList * set = this->Get<DataSetList>("DataSetList");

  set->Add(newSet, newSet->GetAttribute<StringValue>("Name"));
}

DataSet *
BRAINSCutConfiguration::GetAtlasDataSet() const
{
  const DataSetList * set = this->Get<DataSetList>("DataSetList");

  for (ElementParser::const_iterator it = set->begin(); it != set->end(); ++it)
  {
    DataSet * current = dynamic_cast<DataSet *>(it->second);
    if (current->GetAttribute<StringValue>("Type") == "Atlas")
    {
      return current;
    }
  }
  return nullptr;
}

BRAINSCutConfiguration::TrainDataSetListType
BRAINSCutConfiguration::GetTrainDataSets() const
{
  const DataSetList * set = this->Get<DataSetList>("DataSetList");

  std::list<DataSet *> rval;

  for (ElementParser::const_iterator it = set->begin(); it != set->end(); ++it)
  {
    DataSet *   current = dynamic_cast<DataSet *>(it->second);
    std::string type(current->GetAttribute<StringValue>("Type"));
    if (type != "Atlas" && type != "Apply")
    {
      rval.push_back(current);
    }
  }
  return rval;
}

BRAINSCutConfiguration::ApplyDataSetListType
BRAINSCutConfiguration::GetApplyDataSets() const
{
  const DataSetList * set = this->Get<DataSetList>("DataSetList");

  std::list<DataSet *> rval;

  for (ElementParser::const_iterator it = set->begin(); it != set->end(); ++it)
  {
    DataSet * current = dynamic_cast<DataSet *>(it->second);
    if (current->GetAttribute<StringValue>("Type") == "Apply")
    {
      rval.push_back(current);
    }
  }
  return rval;
}

const DataSet::StringVectorType
BRAINSCutConfiguration::GetImageTypes() const
{
  const DataSetList * set = this->Get<DataSetList>("DataSetList");

  if (set->size() == 0)
  {
    return DataSet::StringVectorType();
  }
  return dynamic_cast<const DataSet *>(set->begin()->second)->GetImageTypes();
}

// Set/Get Functions
//
std::string
BRAINSCutConfiguration::GetImageTypeToUse()
{
  return ImageTypeToUse;
}

std::string
BRAINSCutConfiguration::GetRegistrationID()
{
  return RegistrationID;
}
