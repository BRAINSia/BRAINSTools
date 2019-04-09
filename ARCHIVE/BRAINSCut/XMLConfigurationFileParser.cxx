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
#include "XMLConfigurationFileParser.h"
#include "ApplyModel.h"
#include "TrainingVectorConfigurationType.h"
#include "TrainingPrameters.h"
#include <ElementContainer.h>
#include "BRAINSCutExceptionStringHandler.h"

#define CheckDataSet()                              \
  if( dataSet == 0 )                               \
    {                                               \
    std::string s(name);                            \
    s += " occurs in XML file outside any dataset"; \
    BRAINSCutExceptionStringHandler ex( s.c_str() );         \
    throw ex;                                       \
    }

void
XMLConfigurationFileParser::StartElement(void *userData,
                                         const XML_Char *name,
                                         const XML_Char * *atts)
{
  // collect attributes

  StringMap attribMap;

  for( unsigned i = 0; atts[i] != nullptr; i += 2 )
    {
    attribMap[std::string(atts[i])] = std::string(atts[i + 1]);
    }

  // in case of sub-attribute
  std::list<ElementContainer *> *stack =
    static_cast<std::list<ElementContainer *> *>( userData );

  ElementContainer *current = *( stack->begin() );

  DataSet * dataSet = static_cast<DataSet *>( current );

  // name
  const std::string Name(name);

  if( Name == "AutoSegProcessDescription" )
    {
    // nothing to do, top level object is on the top of stack
    return;
    }
  else if( Name == "DataSet"  )
    {
    DataSet *currentDataSet = new DataSet;
    try
      {
      std::string currentDataSetName( attribMap.Get( Name.c_str(), "Name") );
      currentDataSet->SetAttribute<StringValue, std::string>("Name", currentDataSetName);
      currentDataSet->SetAttribute<StringValue, std::string>( "Type",   attribMap.Get( Name.c_str(), "Type") );
      if( currentDataSet->GetAttribute<StringValue>("Type") == "Apply" )
        {
        currentDataSet->SetAttribute<StringValue, std::string>( "OutputDir",
                                                                attribMap.Get( Name.c_str(), "OutputDir") );
        }

      myConfiguration->AddDataSet(currentDataSet);

      stack->push_front(currentDataSet);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "ProbabilityMap" )
    {
    try
      {
      ProbabilityMapList *mapList =
        myConfiguration->Get<ProbabilityMapList>("ProbabilityMapList");
      ProbabilityMapParser *map = new ProbabilityMapParser;
      std::string           structureID( attribMap.Get("ProbabilityMap",
                                                       "StructureID") );
      map->SetAttribute<StringValue>("StructureID",
                                     structureID);
      map->SetAttribute<StringValue>( "Filename",
                                      attribMap.Get("ProbabilityMap",
                                                    "Filename") );
      map->SetAttribute<FloatValue>( "Gaussian",
                                     attribMap.Get("ProbabilityMap",
                                                   "Gaussian") );
      map->SetAttribute<StringValue>( "GenerateVector",
                                      attribMap.Get("GenerateVector",
                                                    "GenerateVector") );
      mapList->Add(map, structureID);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "Registration" )
    {
    try
      {
      CheckDataSet();

      RegistrationList *regList =  dataSet->Get<RegistrationList>(
          "RegistrationList");
      RegistrationType *registration =  new RegistrationType;
      registration->SetAttribute<StringValue>(
        "SubjToAtlasRegistrationFilename",
        attribMap.Get("Registration",
                      "SubjToAtlasRegistrationFilename") );
      registration->SetAttribute<StringValue>(
        "AtlasToSubjRegistrationFilename",
        attribMap.Get("Registration",
                      "AtlasToSubjRegistrationFilename") );
      std::string id( attribMap.Get("Registration", "ID") );
      registration->SetAttribute<StringValue>("ID", id);
      regList->Add(registration, id);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "RegistrationConfiguration" )
    {
    try
      {
      RegistrationConfigurationParser *params =
        myConfiguration->Get<RegistrationConfigurationParser>("RegistrationConfiguration");
      params->SetAttribute<StringValue>( "ImageTypeToUse",
                                         attribMap.Get("RegistrationConfiguration",
                                                       "ImageTypeToUse") );
      params->SetAttribute<StringValue>( "ID",
                                         attribMap.Get("RegistrationConfiguration",
                                                       "ID") );
      params->SetAttribute<IntValue>( "BRAINSROIAutoDilateSize",
                                      attribMap.GetIfExist("RegistrationConfiguration",
                                                           "BRAINSROIAutoDilateSize") );
      params->SetAttribute<BooleanValue>( "ProbabilityMapRegistrationToSubject",
                                      attribMap.GetIfExist("RegistrationConfiguration",
                                                           "ProbabilityMapRegistrationToSubject") );
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "Mask" )
    {
    try
      {
      CheckDataSet();
      MaskType *  mask = new MaskType;
      std::string type( attribMap.Get("Mask", "Type") );
      mask->SetAttribute<StringValue>("Type", type);
      mask->SetAttribute<StringValue>( "Filename",
                                       attribMap.Get("Mask", "Filename") );

      MaskList *maskList = dataSet->Get<MaskList>("MaskList");
      maskList->Add(mask, type);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "Image" )
    {
    try
      {
      CheckDataSet();
      ImageDescription *image = new ImageDescription;
      std::string       type( attribMap.Get("Image", "Type") );
      image->SetAttribute<StringValue>("Type", type);
      image->SetAttribute<StringValue>( "Filename",
                                        attribMap.Get("Image", "Filename") );

      ImageList *imageList = dataSet->Get<ImageList>("ImageList");
      imageList->Add(image, type);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "SpatialLocation" )
    {
    try
      {
      CheckDataSet();
      SpatialLocationType *image = new SpatialLocationType;
      std::string          type( attribMap.Get("SpatialLocation", "Type") );
      image->SetAttribute<StringValue>("Type", type);
      image->SetAttribute<StringValue>( "Filename",
                                        attribMap.Get("SpatialLocation", "Filename") );

      SpatialLocationList *imageList = dataSet->Get<SpatialLocationList>("SpatialLocationList");
      imageList->Add(image, type);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "TrainingVectorConfiguration" || Name == "NeuralNetParams" ) // TODO Change NeuralNet Param
    {
    try
      {
      TrainingVectorConfigurationType *np = new TrainingVectorConfigurationType;
      np->SetAttribute<FloatValue>( "MaskSmoothingValue",
                                    attribMap.Get( Name.c_str(),
                                                   "MaskSmoothingValue") );
      np->SetAttribute<IntValue>( "GradientProfileSize",
                                  attribMap.Get( Name.c_str(),
                                                 "GradientProfileSize") );
      np->SetAttribute<StringValue>( "TrainingVectorFilename",
                                     attribMap.Get( Name.c_str(),
                                                    "TrainingVectorFilename") );
      np->SetAttribute<StringValue>( "TestVectorFilename",
                                     attribMap.Get( Name.c_str(),
                                                    "TestVectorFilename") );
      np->SetAttribute<StringValue>( "TrainingModelFilename",
                                     attribMap.Get( Name.c_str(),
                                                    "TrainingModelFilename") );
      np->SetAttribute<StringValue>( "Normalization",
                                     attribMap.Get( Name.c_str(),
                                                    "Normalization") );
      myConfiguration->Add(np, "TrainingVectorConfiguration");
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "RandomForestParameters" )
    {
    try
      {
      TrainingParameters *ap = new TrainingParameters("RandomForestParameters");
      ap->SetAttribute<IntValue>( "MaxDepth",
                                  attribMap.Get("RandomForestParameters",
                                                "MaxDepth") );
      ap->SetAttribute<IntValue>( "MinSampleCount",
                                  attribMap.Get("RandomForestParameters",
                                                "MinSampleCount") );
      ap->SetAttribute<BooleanValue>( "UseSurrogates",
                                      attribMap.Get("RandomForestParameters",
                                                    "UseSurrogates") );
      ap->SetAttribute<BooleanValue>( "CalcVarImportance",
                                      attribMap.Get("RandomForestParameters",
                                                    "CalcVarImportance") );
      ap->SetAttribute<IntValue>( "MaxTreeCount",
                                  attribMap.Get("RandomForestParameters",
                                                "MaxTreeCount") );
      myConfiguration->Add(ap, Name);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "ANNParameters" )
    {
    try
      {
      TrainingParameters *ap = new TrainingParameters("ANNParameters");
      // ap->SetAttribute<IntValue>( "VectorSize",    attribMap.Get("TrainingParameters",
      //          "VectorSize") );
      ap->SetAttribute<IntValue>( "Iterations",
                                  attribMap.Get("ANNParameters",
                                                "Iterations") );
      ap->SetAttribute<IntValue>( "MaximumVectorsPerEpoch",
                                  attribMap.Get("ANNParameters",
                                                "MaximumVectorsPerEpoch") );
      ap->SetAttribute<IntValue>( "EpochIterations",
                                  attribMap.Get("ANNParameters",
                                                "EpochIterations") );
      ap->SetAttribute<IntValue>( "ErrorInterval",
                                  attribMap.Get("ANNParameters",
                                                "ErrorInterval") );
      ap->SetAttribute<FloatValue>( "ActivationSlope",
                                    attribMap.Get("ANNParameters",
                                                  "ActivationSlope") );
      ap->SetAttribute<FloatValue>( "ActivationMinMax",
                                    attribMap.Get("ANNParameters",
                                                  "ActivationMinMax") );
      ap->SetAttribute<FloatValue>( "DesiredError",
                                    attribMap.Get("ANNParameters",
                                                  "DesiredError") );
      ap->SetAttribute<IntValue>( "NumberOfHiddenNodes",
                                  attribMap.Get("ANNParameters",
                                                "NumberOfHiddenNodes") );
      myConfiguration->Add(ap, Name);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else if( Name == "ApplyModel" )
    {
    try
      {
      ApplyModelType *am = new ApplyModelType;
      am->SetAttribute<FloatValue>( "MaskThresh",
                                    attribMap.Get("ApplyModel",
                                                  "MaskThresh") );
      am->SetAttribute<FloatValue>( "GaussianSmoothingSigma",
                                    attribMap.GetIfExist("ApplyModel",
                                                         "GaussianSmoothingSigma") );

      myConfiguration->Add(am, Name);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw;
      }
    }
  else
    {
    std::string message = "The Element name of \"";
    message += Name;
    message += "\" does not exist. Please check if the xml file is well-formed\n";
    throw BRAINSCutExceptionStringHandler(message);
    }
}

void
XMLConfigurationFileParser::EndElement(void *userData,
                                       const XML_Char *name)
{
  std::list<ElementContainer *> *stack =
    static_cast<std::list<ElementContainer *> *>( userData );

  if( std::string(name) == "DataSet" )
    {
    stack->pop_front();
    }
}

BRAINSCutConfiguration *
XMLConfigurationFileParser::GetConfiguration()
{
  return myConfiguration;
}

/*
void
XMLConfigurationFileParser::ReadXML()
{
  std::list<ElementContainer *> myConfigurationBuffer;
  myConfigurationBuffer.push_front( myConfiguration );

  SetUserData( &myConfigurationBuffer);
  Parse();
}
*/

/**
 * Validation function merged into this class from GenerateProbability class
 */
void
XMLConfigurationFileParser::ValidateDataSets()
{
  // HACK:  Needed to speed up testing.
  // std::list<DataSet *> dataSets = myConfiguration->GetTrainDataSets();

  std::cout << " ***************************************************" << std::endl
            << " Validation has not been implimented yet" << std::endl
            << " ***************************************************" << std::endl;
  // return true;

  /*
   * TODO:: change validation part to check the simple file existance checking
   *
   *
  for( std::list<DataSet *>::iterator it = dataSets.begin();
       it != dataSets.end(); ++it )
    {
    DataSet::StringVectorType allImageTypes( ( *it )->GetImageTypes() );
    const std::string   FirstImageName(
      ( *it )->GetImageFilenameByType(allImageTypes[0]) );
    InternalImageType::Pointer FirstImage =
      itkUtil::ReadImage<InternalImageType>(FirstImageName);
    FirstImage =
      itkUtil::ScaleAndCast<InternalImageType, InternalImageType>(
        FirstImage,
        0,
        HUNDRED_PERCENT_VALUE);
    for( unsigned int imindex = 1; imindex < allImageTypes.size(); imindex++ )
      {
      const std::string CurrentImageName( ( *it )->GetImageFilenameByType(
                                            allImageTypes[imindex]) );
      if( !itksys::SystemTools::FileExists( CurrentImageName.c_str() ) )
        {
        std::string errmsg(CurrentImageName);
        errmsg += " does not exist";
        throw  BRAINSCutExceptionStringHandler(errmsg);
        }
      InternalImageType::Pointer CurrentImage =
        itkUtil::ReadImage<InternalImageType>(CurrentImageName);
      CurrentImage =
        itkUtil::ScaleAndCast<InternalImageType,
                              InternalImageType>(
          CurrentImage,
          0,
          HUNDRED_PERCENT_VALUE);
      if( !itkUtil::ImagePhysicalDimensionsAreIdentical<
            InternalImageType, InternalImageType>(FirstImage,
                                                              CurrentImage) )
        {
        std::string errmsg(CurrentImageName);
        errmsg += " and ";
        errmsg += FirstImageName;
        errmsg += " differ";
        std::cout << "============" << FirstImageName << "===============\n"
                  << FirstImage << std::endl;
        std::cout << "============" << CurrentImageName
                  << "===============\n" << CurrentImage << std::endl;
        throw  BRAINSCutExceptionStringHandler(errmsg);
        }
      } // Each Image

    DataSet::StringVectorType allMaskTypes( ( *it )->GetMaskTypes() );
    for( unsigned maskindex = 0; maskindex < allMaskTypes.size(); maskindex++ )
      {
      const std::string CurrentMaskName( ( *it )->GetMaskFilenameByType(
                                           allMaskTypes[maskindex]) );
      if( !itksys::SystemTools::FileExists( CurrentMaskName.c_str() ) )
        {
        std::string errmsg(CurrentMaskName);
        errmsg += " does not exist";
        throw  BRAINSCutExceptionStringHandler(errmsg);
        }
      InternalImageType::Pointer CurrentMask =
        itkUtil::ReadImage<InternalImageType>(CurrentMaskName);
      CurrentMask =
        itkUtil::ScaleAndCast<InternalImageType,
                              InternalImageType>(
          CurrentMask,
          0,
          HUNDRED_PERCENT_VALUE);
      if( !itkUtil::ImagePhysicalDimensionsAreIdentical<
            InternalImageType, InternalImageType>(FirstImage,
                                                              CurrentMask) )
        {
        std::string errmsg(CurrentMaskName);
        errmsg += " and ";
        errmsg += FirstImageName;
        errmsg += " differ";
        std::cout << "============" << FirstImageName << "===============\n"
                  << FirstImage << std::endl;
        std::cout << "============" << CurrentMaskName
                  << "===============\n" << CurrentMask << std::endl;
        throw  BRAINSCutExceptionStringHandler(errmsg);
        }
      } // Each Mask
    }   // Each DataSet
  return true;
  */
}
