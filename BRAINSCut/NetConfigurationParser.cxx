#include "NetConfigurationParser.h"
#include "ApplyModel.h"
#include "NeuralParams.h"
#include "ANNParams.h"
#include <XMLElementContainer.h>
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
NetConfigurationParser::StartElement(void *userData,
                                     const XML_Char *name,
                                     const XML_Char * *atts)
{
  // collect attributes

  // Use GetIfExist( ) instead of Get() for optional attributes
  StringMap attribMap;

  for( unsigned i = 0; atts[i] != 0; i += 2 )
    {
    attribMap[std::string(atts[i])] = std::string(atts[i + 1]);
    }

  std::string Name(name);

  std::list<XMLElementContainer *> *stack =
    static_cast<std::list<XMLElementContainer *> *>( userData );

  XMLElementContainer *current = *( stack->begin() );

  // only one of these two dynamic casts will succeed, but
  // do them here to avoid duplication below.
  DataSet *         dataSet = dynamic_cast<DataSet *>( current );
  NetConfiguration *Local_netConfiguration = dynamic_cast<NetConfiguration *>( current );

  if( Name == "AutoSegProcessDescription" )
    {
    // nothing to do, top level object is on the top of stack
    return;
    }
  else if( Name == "DataSet" )
    {
    DataSet *currentDataSet = new DataSet;
    try
      {
      std::string currentDataSetName( attribMap.Get("DataSet", "Name") );
      currentDataSet->SetAttribute<StringValue, std::string>("Name", currentDataSetName);
      currentDataSet->SetAttribute<StringValue, std::string>( "Type",   attribMap.Get("DataSet", "Type") );
      if( currentDataSet->GetAttribute<StringValue>("Type") == "Apply" )
        {
        currentDataSet->SetAttribute<StringValue, std::string>( "OutputDir",   attribMap.Get("DataSet", "OutputDir") );
        }

      Local_netConfiguration->AddDataSet(currentDataSet);
      stack->push_front(currentDataSet);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw ex;
      }
    }
  else if( Name == "ProbabilityMap" )
    {
    try
      {
      ProbabilityMapList *mapList =
        Local_netConfiguration->Get<ProbabilityMapList>("ProbabilityMapList");
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
      throw ex;
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
      throw ex;
      }
    }
  else if( Name == "RegistrationConfiguration" )
    {
    try
      {
      RegistrationConfigurationParser *params =
        Local_netConfiguration->Get<RegistrationConfigurationParser>("RegistrationConfiguration");
      params->SetAttribute<StringValue>( "ImageTypeToUse",
                                         attribMap.Get("RegistrationConfiguration",
                                                       "ImageTypeToUse") );
      params->SetAttribute<StringValue>( "ID",
                                         attribMap.Get("RegistrationConfiguration",
                                                       "ID") );
      params->SetAttribute<IntValue>( "BRAINSROIAutoDilateSize",
                                      attribMap.GetIfExist("RegistrationConfiguration",
                                                           "BRAINSROIAutoDilateSize") );
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw ex;
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
      throw ex;
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
      throw ex;
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
      throw ex;
      }
    }
  else if( Name == "NeuralNetParams" )
    {
    try
      {
      NeuralParams *np = new NeuralParams;
      np->SetAttribute<FloatValue>( "MaskSmoothingValue",
                                    attribMap.Get("NeuralNetParams",
                                                  "MaskSmoothingValue") );
      // np->SetAttribute<IntValue>("GaussianSize",
      // attribMap.Get("NeuralNetParams",
      // "GaussianSize"));
      np->SetAttribute<IntValue>( "GradientProfileSize",
                                  attribMap.Get("NeuralNetParams",
                                                "GradientProfileSize") );
      // np->SetAttribute<IntValue>("IrisSize",
      //                       attribMap.Get("NeuralNetParams",
      //                                   "IrisSize"));
      np->SetAttribute<StringValue>( "TrainingVectorFilename",
                                     attribMap.Get("NeuralNetParams",
                                                   "TrainingVectorFilename") );
      np->SetAttribute<StringValue>( "TestVectorFilename",
                                     attribMap.Get("NeuralNetParams",
                                                   "TestVectorFilename") );
      np->SetAttribute<StringValue>( "TrainingModelFilename",
                                     attribMap.Get("NeuralNetParams",
                                                   "TrainingModelFilename") );
      np->SetAttribute<StringValue>( "Normalization",
                                     attribMap.Get("NeuralNetParams",
                                                   "Normalization") );
      Local_netConfiguration->Add(np, Name);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw ex;
      }
    }
  else if( Name == "ANNParams" )
    {
    try
      {
      ANNParams *ap = new ANNParams;
      // ap->SetAttribute<IntValue>( "VectorSize",    attribMap.Get("ANNParams",
      //          "VectorSize") );
      ap->SetAttribute<IntValue>( "Iterations",
                                  attribMap.Get("ANNParams",
                                                "Iterations") );
      ap->SetAttribute<IntValue>( "MaximumVectorsPerEpoch",
                                  attribMap.Get("ANNParams",
                                                "MaximumVectorsPerEpoch") );
      ap->SetAttribute<IntValue>( "EpochIterations",
                                  attribMap.Get("ANNParams",
                                                "EpochIterations") );
      ap->SetAttribute<IntValue>( "ErrorInterval",
                                  attribMap.Get("ANNParams",
                                                "ErrorInterval") );
      ap->SetAttribute<FloatValue>( "ActivationSlope",
                                    attribMap.Get("ANNParams",
                                                  "ActivationSlope") );
      ap->SetAttribute<FloatValue>( "ActivationMinMax",
                                    attribMap.Get("ANNParams",
                                                  "ActivationMinMax") );
      ap->SetAttribute<FloatValue>( "DesiredError",
                                    attribMap.Get("ANNParams",
                                                  "DesiredError") );
      ap->SetAttribute<IntValue>( "NumberOfHiddenNodes",
                                  attribMap.Get("ANNParams",
                                                "NumberOfHiddenNodes") );
      Local_netConfiguration->Add(ap, Name);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw ex;
      }
    }
  else if( Name == "ApplyModel" )
    {
    try
      {
      ApplyModelType *am = new ApplyModelType;
      am->SetAttribute<FloatValue>( "CutOutThresh",
                                    attribMap.Get("ApplyModel",
                                                  "CutOutThresh") );
      am->SetAttribute<FloatValue>( "MaskThresh",
                                    attribMap.Get("ApplyModel",
                                                  "MaskThresh") );
      Local_netConfiguration->Add(am, Name);
      }
    catch( BRAINSCutExceptionStringHandler& ex )
      {
      std::cerr << ex << std::endl;
      throw ex;
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
NetConfigurationParser::EndElement(void *userData,
                                   const XML_Char *name)
{
  std::list<XMLElementContainer *> *stack =
    static_cast<std::list<XMLElementContainer *> *>( userData );

  if( std::string(name) == "DataSet" )
    {
    stack->pop_front();
    }
}

NetConfiguration *
NetConfigurationParser::GetNetConfiguration()
{
  return netConfiguration;
}

void
NetConfigurationParser::ReadXML()
{
  std::list<XMLElementContainer *> netConfigurationBuffer;

  netConfigurationBuffer.push_front( netConfiguration );

  SetUserData( &netConfigurationBuffer);
  Parse();
}

/**
 * Validation function merged into this class from GenerateProbability class
 */
void
NetConfigurationParser::ValidateDataSets()
{
  // HACK:  Needed to speed up testing.
  // std::list<DataSet *> dataSets = netConfiguration->GetTrainDataSets();

  std::cout << " ***************************************************" << std::endl
            << " Validation has not been implimented yet" << std::endl
            << " ***************************************************" << std::endl;
  // return true;

  /*
   * chage validation part to check the simple file existance checking
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
