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
#ifndef DataSet_h
#define DataSet_h
// #include "StringValue.h"
#include "ImageDescription.h"
#include "MaskType.h"
#include "ElementParser.h"
#include "RegistrationType.h"
#include "SpatialLocationType.h"
#include <vector>
//
// subclass ElementParser, overriding the
// constructor to initialize the model structure.
class DataSet : public ElementParser
{
public:
  typedef ElementParser SuperClass;

  int PrintSelf(std::ostream & os, int indent) const override
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== DataSet ===" << std::endl;
    return indent + 2;
  }

  typedef SuperClass::StringVectorType StringVectorType;

  DataSet() : ElementParser("DataSet")
  {
    this->Add(new StringValue("Name", ""), "Name");
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("OutputDir", "na"), "OutputDir");
    this->Add(new ImageList, "ImageList");
    this->Add(new SpatialLocationList, "SpatialLocationList");
    this->Add(new MaskList, "MaskList");
    this->Add(new RegistrationList, "RegistrationList");
  }

  //
  // get the filename out of the list specified by listname.
  // templated over the list type and element typee.
  //
  template <typename ListType, typename ElementType>
  std::string GetFilenameByType(const char *listName,  const char *type) const
  {
    const ListType *   list = this->Get<ListType>(listName);
    const ElementType *element    = list->template GetMatching<ElementType>("Type", type);

    return element == nullptr ? std::string() :
           element->template GetAttribute<StringValue>("Filename");
  }

  //
  // get the image with the type 'type'
  const std::string GetImageFilenameByType(const char *type) const
  {
    return this->GetFilenameByType<ImageList, ImageDescription>("ImageList",
                                                                type);
  }

  const std::string GetImageFilenameByType(const std::string & type) const
  {
    return GetImageFilenameByType( type.c_str() );
  }

  //
  // get the SpatialLocation with the type 'type'
  const std::string GetSpatialLocationFilenameByType(const char *type) const
  {
    return this->GetFilenameByType<SpatialLocationList, SpatialLocationType>("SpatialLocationList",
                                                                             type);
  }

  const std::string GetSpatialLocationFilenameByType(const std::string & type) const
  {
    return GetSpatialLocationFilenameByType( type.c_str() );
  }

  //
  // get the Mask with the type 'type'
  const std::string GetMaskFilenameByType(const char *type) const
  {
    return this->GetFilenameByType<MaskList, MaskType>("MaskList",
                                                       type);
  }

  const std::string GetMaskFilenameByType(const std::string & type) const
  {
    return GetMaskFilenameByType( type.c_str() );
  }

  //
  // get the Registration with the ID 'id'
  const RegistrationType * GetRegistrationWithID(const char *id) const
  {
    const RegistrationList *list =
      this->Get<RegistrationList>("RegistrationList");

    return list->GetMatching<RegistrationType>("ID", id);
  }

  const RegistrationType * GetRegistrationWithID(const std::string & id) const
  {
    return GetRegistrationWithID( id.c_str() );
  }

  //
  // get all potential image types
  const StringVectorType GetImageTypes() const
  {
    const ImageList *imList = this->Get<ImageList>("ImageList");

    return imList->CollectAttValues<ImageDescription>("Type");
  }

  //
  // get all spatialLocation types
  const StringVectorType GetSpatialLocationTypes() const
  {
    const SpatialLocationList *spatialLocationList = this->Get<SpatialLocationList>("SpatialLocationList");

    return spatialLocationList->CollectAttValues<SpatialLocationType>("Type");
  }

  //
  // get all mask types
  const StringVectorType GetMaskTypes() const
  {
    const MaskList *maskList = this->Get<MaskList>("MaskList");

    return maskList->CollectAttValues<MaskType>("Type");
  }
};

class DataSetList : public ElementParser
{
public:
  DataSetList() : ElementParser("DataSetList")
  {
  }
};

#endif
