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
#ifndef ProcessDescription_h
#define ProcessDescription_h
#include <CompoundObjectBase.h>
#include <DataSet.h>
#include <RegistrationParams.h>
#include <ProbabilityMap.h>
#include <list>

class ProcessDescription : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    return indent + 2;
  }

  typedef std::list<DataSet *> TrainDataSetListType;
  ProcessDescription() : CompoundObjectBase("BRAINSCutProcessDescription")
  {
    this->Add(new DataSetList, "DataSetList");
    this->Add(new ProbabilityMapList, "ProbabilityMapList");
    this->Add(new RegistrationParams, "RegistrationParams");
  }

  void AddDataSet(DataSet *newSet)
  {
    DataSetList *set = this->Get<DataSetList>("DataSetList");

    set->Add( newSet,
              newSet->GetAttribute<StringValue>("Name") );
  }

  DataSet * GetAtlasDataSet() const
  {
    const DataSetList *set = this->Get<DataSetList>("DataSetList");

    for( CompoundObjectBase::const_iterator it = set->begin();
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

  TrainDataSetListType GetTrainDataSets() const
  {
    const DataSetList *set = this->Get<DataSetList>("DataSetList");

    std::list<DataSet *> rval;

    for( CompoundObjectBase::const_iterator it = set->begin();
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

  TrainDataSetListType GetApplyDataSets() const
  {
    const DataSetList *set = this->Get<DataSetList>("DataSetList");

    std::list<DataSet *> rval;

    for( CompoundObjectBase::const_iterator it = set->begin();
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

  const DataSet::TypeVector ImageTypes() const
  {
    const DataSetList *set = this->Get<DataSetList>("DataSetList");

    if( set->size() == 0 )
      {
      return DataSet::TypeVector();
      }
    return dynamic_cast<const DataSet *>( set->begin()->second )->ImageTypes();
  }
};

#endif // ProcessDescription_h
