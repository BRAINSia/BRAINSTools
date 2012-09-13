#include "BRAINSCutDataHandler.h"
#include "XMLConfigurationFileParser.h"
#include "GenericTransformImage.h"

/** constructors */
BRAINSCutDataHandler
::BRAINSCutDataHandler( std::string modelConfigurationFilename )
{
  try
    {
    SetNetConfigurationFilename( modelConfigurationFilename );
    SetNetConfiguration();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_FAILURE);
    }
}

void
BRAINSCutDataHandler
::SetNetConfiguration()
{
  try
    {
    // read in xml images
    std::list<ElementContainer *> elementList;

    elementList.push_front( myConfiguration );

    XMLConfigurationFileParser BRAINSCutXMLConfigurationFileParser =
      XMLConfigurationFileParser( myConfigurationFilename );
    BRAINSCutXMLConfigurationFileParser.SetUserData( &elementList );
    BRAINSCutXMLConfigurationFileParser.Parse();

    myConfiguration = BRAINSCutXMLConfigurationFileParser.GetConfiguration();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_FAILURE);
    }
}

void
BRAINSCutDataHandler
::SetAtlasDataSet()
{
  atlasDataSet = myConfiguration->GetAtlasDataSet();
  std::cout << "registrationImageTypeToUse :: " << registrationImageTypeToUse << std::endl;

  if( registrationImageTypeToUse.empty() )
    {
    std::cout << "registrationImageTypeToUse is empty." << std::endl;
    exit(EXIT_FAILURE);
    }
  atlasFilename = atlasDataSet->GetImageFilenameByType( registrationImageTypeToUse);
  atlasBinaryFilename = atlasDataSet->GetMaskFilenameByType( "RegistrationROI" );
  std::cout << atlasBinaryFilename << std::endl;
}

void
BRAINSCutDataHandler
::SetAtlasImage()
{
  atlasImage = ReadImageByFilename( atlasFilename );
}

void
BRAINSCutDataHandler
::SetRhoPhiTheta()
{
  rho = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("rho") );
  phi = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("phi") );
  theta = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("theta") );
}

void
BRAINSCutDataHandler
::SetNetConfigurationFilename(const std::string filename)
{
  myConfigurationFilename = filename;
}

std::string
BRAINSCutDataHandler
::GetNetConfigurationFilename()
{
  return myConfigurationFilename;
}

int
BRAINSCutDataHandler
::GetTrainIteration()
{
  int trainIteration = myConfiguration->Get<TrainingParameters>("ANNParameters")
    ->GetAttribute<IntValue>("Iterations");

  return trainIteration;
}

DataSet::StringVectorType
BRAINSCutDataHandler
::GetROIIDsInOrder() const
{
  return this->m_roiIDsInOrder;
}

void
BRAINSCutDataHandler
::SetRegionsOfInterest()
{
  roiDataList = myConfiguration->Get<ProbabilityMapList>("ProbabilityMapList");
  this->m_roiIDsInOrder = roiDataList->CollectAttValues<ProbabilityMapParser>("StructureID");

  std::sort( this->m_roiIDsInOrder.begin(), this->m_roiIDsInOrder.end() ); // get l_caudate, l_globus, .. , r_caudate,
                                                                           // r_globus..
  roiCount = roiDataList->size();
}

/** registration related */
void
BRAINSCutDataHandler
::SetRegistrationParameters()
{
  registrationParser =
    myConfiguration->Get<RegistrationConfigurationParser>("RegistrationConfiguration");

  SetRegistrationImageTypeToUse(
    std::string( registrationParser->GetAttribute<StringValue>( "ImageTypeToUse") ) );

  registrationID = std::string(
      registrationParser->GetAttribute<StringValue>("ID") );

  roiAutoDilateSize = registrationParser->GetAttribute<IntValue>("BRAINSROIAutoDilateSize");
}

void
BRAINSCutDataHandler
::SetRegistrationImageTypeToUse( std::string type )
{
  registrationImageTypeToUse = type;
}

std::string
BRAINSCutDataHandler
::GetRegistrationImageTypeToUse()
{
  if( registrationImageTypeToUse.empty() )
    {
    std::cout << "registrationImageTypeToUse is empty." << std::endl;
    exit(EXIT_FAILURE);
    }
  return registrationImageTypeToUse;
}

std::string
BRAINSCutDataHandler
::GetRegistrationID()
{
  if( registrationID.empty() )
    {
    std::cout << "The registrationID is empty" << std::endl;
    exit(EXIT_FAILURE);
    }
  return registrationID;
}

std::string
BRAINSCutDataHandler
::GetSubjectToAtlasRegistrationFilename( DataSet& subject)
{
  std::string filename = subject.GetRegistrationWithID( registrationID )
    ->GetAttribute<StringValue>("SubjToAtlasRegistrationFilename");

  return filename;
}

std::string
BRAINSCutDataHandler
::GetAtlasToSubjectRegistrationFilename( DataSet& subject)
{
  std::string filename = subject.GetRegistrationWithID( registrationID )
    ->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename");

  return filename;
}

void
BRAINSCutDataHandler
::GetDeformedSpatialLocationImages( std::map<std::string, WorkingImagePointer>& warpedSpatialLocationImages,
                                    DataSet& subject)
{
  std::string atlasSubjectRegistrationFilename = GetAtlasToSubjectRegistrationFilename( subject );

  DisplacementFieldType::Pointer deformation = GetDeformationField( atlasSubjectRegistrationFilename );
  GenericTransformType::Pointer  genericTransform = GetGenericTransform( atlasSubjectRegistrationFilename );

  std::string subjectFilenameToUse = subject.GetImageFilenameByType(registrationImageTypeToUse);

  if( !itksys::SystemTools::FileExists( subjectFilenameToUse.c_str(), false ) )
    {
    std::cout << " Subject image file of "
              << subjectFilenameToUse
              << ", type of " << registrationImageTypeToUse
              << ", does not exist. "
              << std::endl;
    exit(EXIT_FAILURE);
    }

  WorkingImagePointer referenceImage =
    ReadImageByFilename( subject.GetImageFilenameByType(registrationImageTypeToUse) );

  const std::string transoformationPixelType = "float";

  warpedSpatialLocationImages.insert( std::pair<std::string, WorkingImagePointer>
                                        ("rho", GenericTransformImage<WorkingImageType,
                                                                      WorkingImageType,
                                                                      DisplacementFieldType>
                                          ( rho, referenceImage, deformation, genericTransform,
                                          0.0, "Linear", transoformationPixelType == "binary") ) );
  warpedSpatialLocationImages.insert( std::pair<std::string, WorkingImagePointer>
                                        ("phi", GenericTransformImage<WorkingImageType,
                                                                      WorkingImageType,
                                                                      DisplacementFieldType>
                                          ( phi, referenceImage, deformation, genericTransform,
                                          0.0, "Linear", transoformationPixelType == "binary") ) );
  warpedSpatialLocationImages.insert( std::pair<std::string, WorkingImagePointer>
                                        ("theta", GenericTransformImage<WorkingImageType,
                                                                        WorkingImageType,
                                                                        DisplacementFieldType>
                                          ( theta, referenceImage, deformation, genericTransform,
                                          0.0, "Linear", transoformationPixelType == "binary") ) );
}

void
BRAINSCutDataHandler
::GetImagesOfSubjectInOrder( WorkingImageVectorType& subjectImageList, DataSet& subject)
{
  DataSet::StringVectorType imageListFromAtlas = atlasDataSet->GetImageTypes(); // T1, T2, SG, ...

  std::sort( imageListFromAtlas.begin(), imageListFromAtlas.end() );            // SG, T1, T2, ... ascending order

  for( DataSet::StringVectorType::iterator imgTyIt = imageListFromAtlas.begin();
       imgTyIt != imageListFromAtlas.end();
       ++imgTyIt ) // imgTyIt = image type iterator
    {
    std::cout << *imgTyIt << std::endl;
    WorkingImagePointer currentTypeImage = ReadImageByFilename( subject.GetImageFilenameByType( *imgTyIt ) );
    subjectImageList.push_back( currentTypeImage );
    }
}

void
BRAINSCutDataHandler
::GetDeformedROIs( std::map<std::string, WorkingImagePointer>& warpedROIs,
                   DataSet& subject)
{
  std::string atlasSubjectRegistrationFilename = GetAtlasToSubjectRegistrationFilename( subject );

  /** Get the transformation file
   * Note that only one of transformation type will be used. Either deformation or transformation
   * That determined based on the file name at the GetDeformationField
   */
  DisplacementFieldType::Pointer deformation = GetDeformationField( atlasSubjectRegistrationFilename );
  GenericTransformType::Pointer  genericTransform = GetGenericTransform( atlasSubjectRegistrationFilename );

  WorkingImagePointer referenceImage =
    ReadImageByFilename( subject.GetImageFilenameByType(registrationImageTypeToUse) );

  const std::string transformationPixelType = "float";

  for( DataSet::StringVectorType::iterator roiTyIt = this->m_roiIDsInOrder.begin();
       roiTyIt != this->m_roiIDsInOrder.end();
       ++roiTyIt )
    {
    std::string roiFilename = roiDataList->GetMatching<ProbabilityMapParser>( "StructureID", (*roiTyIt).c_str() )
      ->GetAttribute<StringValue>("Filename");
    WorkingImagePointer currentROI = ReadImageByFilename( roiFilename );

    warpedROIs.insert( std::pair<std::string, WorkingImagePointer>(
                         (*roiTyIt), GenericTransformImage<WorkingImageType,
                                                           WorkingImageType,
                                                           DisplacementFieldType>
                           ( currentROI, referenceImage, deformation,
                           genericTransform, 0.0, "Linear",
                           transformationPixelType == "binary") ) );
    }
}

void
BRAINSCutDataHandler
::SetTrainingVectorConfiguration()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  trainingVectorConfiguration = myConfiguration->Get<TrainingVectorConfigurationType>("TrainingVectorConfiguration");
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
}

void
BRAINSCutDataHandler
::SetGradientSize()
{
  gradientSize = trainingVectorConfiguration->GetAttribute<IntValue>("GradientProfileSize");
}

unsigned int
BRAINSCutDataHandler
::GetGradientSize()
{
  return gradientSize;
}

void
BRAINSCutDataHandler
::SetNormalization()
{
  std::string normalizationString;

  try
    {
    normalizationString = trainingVectorConfiguration->GetAttribute<StringValue>("Normalization");
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_FAILURE);
    }

  if( normalizationString == "true" )
    {
    normalization = true;
    }
  else
    {
    normalization =  false;
    }
  std::cout << "Get Normalization from XML file --> " << normalization << std::endl;
}

bool
BRAINSCutDataHandler
::GetNormalization()
{
  return normalization;
}

/** model file name */
std::string
BRAINSCutDataHandler
::GetModelBaseName()
{
  std::string basename;

  try
    {
    basename = trainingVectorConfiguration->GetAttribute<StringValue>("TrainingModelFilename");
    }
  catch( ... )
    {
    std::cout << "Fail to get the ann model file name:: " << basename << std::endl;
    throw BRAINSCutExceptionStringHandler("Fail to get the ann model file name");
    exit(EXIT_FAILURE);
    }
  return basename;
}

std::string
BRAINSCutDataHandler
::GetANNModelFilename()
{
  return ANNModelFilename;
}

std::string
BRAINSCutDataHandler
::GetANNModelFilenameAtIteration( const int iteration)
{
  SetANNModelFilenameAtIteration( iteration );
  return ANNModelFilename;
}

void
BRAINSCutDataHandler
::SetANNModelFilenameAtIteration( const int iteration)
{
  ANNModelFilename = GetModelBaseName();

  char temp[10];
  sprintf( temp, "%09d", iteration );
  ANNModelFilename += temp;
}

std::string
BRAINSCutDataHandler
::GetRFModelFilename( int depth,
                      int NTrees)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  std::string basename = GetModelBaseName();

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  char tempDepth[5];

  sprintf( tempDepth, "%04u", depth );
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  char tempNTrees[5];
  sprintf( tempNTrees, "%04u", NTrees );
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  std::string filename = basename + "D" + tempDepth + "NT" + tempNTrees + ".gz";

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  return filename;
}

std::string
BRAINSCutDataHandler
::GetAtlasFilename()
{
  return atlasFilename;
}

std::string
BRAINSCutDataHandler
::GetAtlasBinaryFilename()
{
  return atlasBinaryFilename;
}

int
BRAINSCutDataHandler
::GetROIAutoDilateSize()
{
  return roiAutoDilateSize;
}

unsigned int
BRAINSCutDataHandler
::GetROICount()
{
  return roiCount;
}

WorkingImagePointer
BRAINSCutDataHandler
::GetAtlasImage()
{
  return atlasImage;
}

ProbabilityMapList *
BRAINSCutDataHandler
::GetROIDataList()
{
  return roiDataList;
}

DataSet *
BRAINSCutDataHandler
::GetAtlasDataSet()
{
  return atlasDataSet;
}

scalarType
BRAINSCutDataHandler
::GetANNOutputThreshold()
{
  scalarType annOutputThreshold =
    myConfiguration->Get<ApplyModelType>("ApplyModel")
    ->GetAttribute<FloatValue>("MaskThresh");

  if( annOutputThreshold < 0.0F )
    {
    std::string msg = " ANNOutput Threshold cannot be less than zero. \n";
    throw BRAINSCutExceptionStringHandler( msg );
    }
  return annOutputThreshold;
}

//
// Apply related
//
void
BRAINSCutDataHandler
::SetRandomForestModelFilename(std::string name)
{
  RandomForestModelFilename = name;
}

void
BRAINSCutDataHandler
::SetRandomForestModelFilename(int depth, int nTree)
{
  if( depth < 0 && nTree < 0 )
    {
    nTree = myConfiguration->Get<TrainingParameters>("RandomForestParameters")
      ->GetAttribute<IntValue>("MaxDepth");
    depth = myConfiguration->Get<TrainingParameters>("RandomForestParameters")
      ->GetAttribute<IntValue>("MaxTreeCount");
    }
  RandomForestModelFilename =  GetRFModelFilename( depth, nTree );
}

std::string
BRAINSCutDataHandler
::GetRandomForestModelFilename()
{
  return RandomForestModelFilename;
}

void
BRAINSCutDataHandler
::SetANNTestingSSEFilename()
{
  ANNTestingSSEFilename = GetModelBaseName();
  ANNTestingSSEFilename += "ValidationSetSSE.txt";
}

std::string
BRAINSCutDataHandler
::GetANNTestingSSEFilename()
{
  return ANNTestingSSEFilename;
}

void
BRAINSCutDataHandler
::SetTrainVectorFilename()
{
  trainVectorFilename =
    trainingVectorConfiguration->GetAttribute<StringValue>("TrainingVectorFilename");
  trainVectorFilename += "ANN"; // TODO
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "vector file at " << trainVectorFilename << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
}

std::string
BRAINSCutDataHandler
::GetTrainVectorFilename()
{
  if( trainVectorFilename.empty() )
    {
    std::string msg = "The train vector file name is empty.\n";
    throw BRAINSCutExceptionStringHandler( msg );
    }
  return trainVectorFilename;
}

BRAINSCutConfiguration::TrainDataSetListType
BRAINSCutDataHandler
::GetTrainDataSet()
{
  BRAINSCutConfiguration::TrainDataSetListType trainDataSetList;

  try
    {
    trainDataSetList = myConfiguration->GetTrainDataSets();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_SUCCESS);
    }
  return trainDataSetList;
}

BRAINSCutConfiguration::ApplyDataSetListType
BRAINSCutDataHandler
::GetApplyDataSet()
{
  BRAINSCutConfiguration::ApplyDataSetListType applyDataSetList;

  try
    {
    applyDataSetList = myConfiguration->GetApplyDataSets();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_SUCCESS);
    }
  return applyDataSetList;
}

scalarType
BRAINSCutDataHandler
::GetGaussianSmoothingSigma()
{
  scalarType gaussianSmoothingSigma =
    myConfiguration->Get<ApplyModelType>("ApplyModel")
    ->GetAttribute<FloatValue>("GaussianSmoothingSigma");

  return gaussianSmoothingSigma;
}

void
BRAINSCutDataHandler
::SetTrainConfiguration( std::string trainParameterName )
{
  TrainConfiguration = myConfiguration->Get<TrainingParameters>( trainParameterName.c_str() );
}

unsigned int
BRAINSCutDataHandler
::GetEpochIteration()
{
  unsigned int trainEpochIteration = TrainConfiguration->GetAttribute<IntValue>("EpochIterations");

  return trainEpochIteration;
}

float
BRAINSCutDataHandler
::GetDesiredError()
{
  float trainDesiredError = TrainConfiguration->GetAttribute<FloatValue>("DesiredError");

  return trainDesiredError;
}

unsigned int
BRAINSCutDataHandler
::GetMaximumDataSize()
{
  unsigned int trainMaximumDataSize = TrainConfiguration->GetAttribute<IntValue>("MaximumVectorsPerEpoch");

  return trainMaximumDataSize;
}

int
BRAINSCutDataHandler
::GetANNHiddenNodesNumber()
{
  int ANNHiddenNodesNumber = TrainConfiguration->GetAttribute<IntValue>("NumberOfHiddenNodes");

  return ANNHiddenNodesNumber;
}

float
BRAINSCutDataHandler
::GetActivationFunctionSlope()
{
  float activationSlope = TrainConfiguration->GetAttribute<FloatValue>("ActivationSlope");

  return activationSlope;
}

float
BRAINSCutDataHandler
::GetActivationFunctionMinMax()
{
  float activationMinMax = TrainConfiguration->GetAttribute<FloatValue>("ActivationMinMax");

  return activationMinMax;
}

/** Get Random Tree */
int
BRAINSCutDataHandler
::GetMaxDepth()
{
  int trainMaxDepth = TrainConfiguration->GetAttribute<IntValue>("MaxDepth");

  return trainMaxDepth;
}

int
BRAINSCutDataHandler
::GetMinSampleCount()
{
  int trainMinSampleCount = TrainConfiguration->GetAttribute<IntValue>("MinSampleCount");

  return trainMinSampleCount;
}

bool
BRAINSCutDataHandler
::GetUseSurrogates()
{
  bool trainUseSurrogates = TrainConfiguration->GetAttribute<BooleanValue>("UseSurrogates");

  return trainUseSurrogates;
}

bool
BRAINSCutDataHandler
::GetCalcVarImportance()
{
  bool trainCalcVarImportance = TrainConfiguration->GetAttribute<BooleanValue>("CalcVarImportance");

  return trainCalcVarImportance;
}

int
BRAINSCutDataHandler
::GetMaxTreeCount()
{
  int trainMaxTreeCount = TrainConfiguration->GetAttribute<IntValue>("MaxTreeCount");

  return trainMaxTreeCount;
}
