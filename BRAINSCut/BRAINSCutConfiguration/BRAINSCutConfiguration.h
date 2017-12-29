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
#ifndef BRAINSCutConfiguration_h
#define BRAINSCutConfiguration_h

#include <DataSet.h>
#include <RegistrationConfigurationParser.h>
#include <ProbabilityMapParser.h>
#include <list>
#include <SpatialLocationType.h>

class BRAINSCutConfiguration : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  BRAINSCutConfiguration();

  virtual int PrintSelf(std::ostream & os, int indent) const override
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

  // Set/Get Functions
  //
  void SetImageTypeToUse( std::string imageTypeToUse );

  std::string GetImageTypeToUse();

  void SetRegistrationID( std::string registrationID );

  std::string GetRegistrationID();

private:
  std::string ImageTypeToUse;
  std::string RegistrationID;
};

#endif // BRAINSCutConfiguration_h
