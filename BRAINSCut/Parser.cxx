#include "Parser.h"
#include "NeuralParams.h"
#include "ANNParams.h"
#include "SVMParams.h"
#include "ApplyModel.h"
#include <ProcessObjectBase.h>

#define CheckDataSet()                              \
  if( dataSet == 0 )                               \
    {                                               \
    std::string s(name);                            \
    s += " occurs in XML file outside any dataset"; \
    ProcessObjectException ex( s.c_str() );         \
    throw ex;                                       \
    }

void
Parser::StartElement(void *userData,
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

  std::string                     Name(name);
  std::list<ProcessObjectBase *> *stack =
    static_cast<std::list<ProcessObjectBase *> *>( userData );
  ProcessObjectBase *current = *( stack->begin() );
  // only one of these two dynamic casts will succeed, but
  // do them here to avoid duplication below.
  DataSet *           dataSet = dynamic_cast<DataSet *>( current );
  ProcessDescription *proc = dynamic_cast<ProcessDescription *>( current );
  if( Name == "AutoSegProcessDescription" )
    {
    // nothing to do, top level object is on the top of stack
    return;
    }
  else if( Name == "DataSet" )
    {
    DataSet *_dataSet = new DataSet;
    try
      {
      std::string _Name( attribMap.Get("DataSet", "Name") );
      _dataSet->SetAttribute<StringValue, std::string>("Name", _Name);
      _dataSet->SetAttribute<StringValue, std::string>( "Type",   attribMap.Get("DataSet", "Type") );
      if( _dataSet->GetAttribute<StringValue>("Type") == "Apply" )
        {
        _dataSet->SetAttribute<StringValue, std::string>( "OutputDir",   attribMap.Get("DataSet", "OutputDir") );
        }

      proc->AddDataSet(_dataSet);
      stack->push_front(_dataSet);
      }
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
      throw;
      }
    }
  else if( Name == "ProbabilityMap" )
    {
    try
      {
      ProbabilityMapList *mapList =
        proc->Get<ProbabilityMapList>("ProbabilityMapList");
      ProbabilityMap *map = new ProbabilityMap;
      std::string     structureID( attribMap.Get("ProbabilityMap",
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
      /*REGINA::
        * Add xml attribute for spherical coordinate system
        */
      map->SetAttribute<StringValue>( "rho",
                                      attribMap.Get("ProbabilityMap",
                                                    "rho") );
      map->SetAttribute<StringValue>( "phi",
                                      attribMap.Get("ProbabilityMap",
                                                    "phi") );
      map->SetAttribute<StringValue>( "theta",
                                      attribMap.Get("ProbabilityMap",
                                                    "theta") );
      mapList->Add(map, structureID);
      }
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
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
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
      throw;
      }
    }
  else if( Name == "RegistrationParams" )
    {
    try
      {
      RegistrationParams *params =
        proc->Get<RegistrationParams>("RegistrationParams");
      params->SetAttribute<StringValue>( "ImageTypeToUse",
                                         attribMap.Get("RegistrationParams",
                                                       "ImageTypeToUse") );
      params->SetAttribute<StringValue>( "ID",
                                         attribMap.Get("RegistrationParams",
                                                       "ID") );
      }
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
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
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
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
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
      throw;
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
      proc->Add(np, Name);
      }
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg
        +=
          "\" is not proper, Please check if the xml file is well-formed\n";
      msg
        +=
          " Check if Iris Size is included, we are not using this parameter anymore. If it does, simply get rid of that part.\n";
      msg +=     "*****************************************************\n";
      std::cerr
        << "********************  ERROR!!! **********************\n " << msg;
      throw;
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
      // ap->SetAttribute<FloatValue>("LearningRate",
      //                           attribMap.Get("ANNParams",
      //                                         "LearningRate"));
      // ap->SetAttribute<FloatValue>("MomentumRate",
      //                           attribMap.Get("ANNParams",
      //                                         "MomentumRate"));
      proc->Add(ap, Name);
      }
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
      throw;
      }
    }
  else if( Name == "SVMParams" )
    {
    try
      {
      SVMParams *sp = new SVMParams;
      sp->SetAttribute<IntValue>( "VectorSize",
                                  attribMap.Get("SVMParams",
                                                "VectorSize") );
      sp->SetAttribute<FloatValue>( "GaussianSize",
                                    attribMap.Get("SVMParams",
                                                  "GaussianSize") );
      proc->Add(sp, Name);
      }
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
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
      proc->Add(am, Name);
      }
    catch( ProcessObjectException& ex )
      {
      std::string msg = " ERROR :: One of attribute name of the \"";
      msg += Name;
      msg += "\" is not proper, Please check if the xml file is well-formed\n";
      std::cerr << "****ERROR!!!\n " << msg;
      throw;
      }
    }
  else
    {
    std::string message = "The Element name of \"";
    message += Name;
    message += "\" is not exist. Please check if the xml file is well-formed\n";
    std::cerr << "****ERROR!!!\n " << message;
    ProcessObjectException exception(message);
    throw exception;
    }
}

void
Parser::EndElement(void *userData,
                   const XML_Char *name)
{
  std::list<ProcessObjectBase *> *stack =
    static_cast<std::list<ProcessObjectBase *> *>( userData );

  if( std::string(name) == "DataSet" )
    {
    stack->pop_front();
    }
}

bool ReadXML(const char *filename, ProcessDescription & prob)
{
  std::list<ProcessObjectBase *> stack;

  stack.push_front(&prob);
  Parser parser(filename);
  parser.SetUserData(&stack);
  return parser.Parse();
}
