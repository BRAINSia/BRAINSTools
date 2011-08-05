#ifndef DataSet_h
#define DataSet_h
#include "StringValue.h"
#include "ImageDescription.h"
// #include "LandmarkType.h"
#include "MaskType.h"
#include "CompoundObjectBase.h"
#include "RegistrationType.h"
#include <vector>
//
// subclass CompoundObjectBase, overriding the
// constructor to initialize the model structure.
class DataSet : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== DataSet ===" << std::endl;
    return indent + 2;
  }

  typedef SuperClass::StringVectorType TypeVector;

  DataSet() : CompoundObjectBase("DataSet")
  {
    this->Add(new StringValue("Name", ""), "Name");
    this->Add(new StringValue("Type", ""), "Type");
    this->Add(new StringValue("OutputDir", "na"), "OutputDir");
    this->Add(new ImageList, "ImageList");
    //    this->Add(new LandmarkList,"LandmarkList");
    this->Add(new MaskList, "MaskList");
    this->Add(new RegistrationList, "RegistrationList");
  }

  //
  // get the filename out of the list specified by listname.
  // templated over the list type and element typee.
  //
  template <class ListType, class ElementType>
  std::string GetFilenameByType(const char *listName,  const char *type) const
  {
    const ListType *   list = this->Get<ListType>(listName);
    const ElementType *element    = list->GetMatching<ElementType>("Type", type);

    return element == 0 ? std::string() :
           element->GetAttribute<StringValue>("Filename");
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

#if 0
  //
  // get the landmark with the type 'type'
  const std::string GetAtlasFilenameByType(const char *type) const
  {
    return this->GetFilenameByType<LandmarkList, LandmarkType>("LandmarkList",
                                                               type);
  }

  const std::string GetAtlasFilenameByType(const std::string & type) const
  {
    return GetAtlasFilenameByType( type.c_str() );
  }

#endif

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
  const TypeVector ImageTypes() const
  {
    const ImageList *imList = this->Get<ImageList>("ImageList");

    return imList->CollectAttValues<ImageDescription>("Type");
  }

#if 0
  //
  // get all potential landmark types
  const TypeVector LandmarkTypes() const
  {
    LandmarkList *landmarkList = this->Get<LandmarkList>("LandmarkList");

    return landmarkList->CollectAttValues<LandmarkType>("Type");
  }

#endif
  //
  // get all mask types
  const TypeVector MaskTypes() const
  {
    const MaskList *maskList = this->Get<MaskList>("MaskList");

    return maskList->CollectAttValues<MaskType>("Type");
  }
};

class DataSetList : public CompoundObjectBase
{
public:
  DataSetList() : CompoundObjectBase("DataSetList")
  {
  }
};

#endif
