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
#include "BRAINSCutDataHandler.h"
#include "XMLConfigurationFileParser.h"
#include "GenericTransformImage.h"
#include "itkDisplacementFieldTransform.h"

/** constructors */
BRAINSCutDataHandler ::BRAINSCutDataHandler()
  : trainingVectorConfiguration( nullptr )
  , TrainConfiguration( nullptr )
  , m_atlasDataSet( nullptr )
  , m_atlasFilename( "" )
  , m_atlasBinaryFilename( "" )
  , m_atlasImage( nullptr )
  , m_roiDataList( nullptr )
  , m_roiIDsInOrder()
  , roiCount( 0 )
  , ROIRegistrationToSubject( true )
  , registrationParser( nullptr )
  , registrationImageTypeToUse( "" )
  , registrationID( "" )
  , roiAutoDilateSize( 0 )
  , m_rho( nullptr )
  , m_phi( nullptr )
  , m_theta( nullptr )
  , m_gradientSize( 0 )
  , m_trainVectorFilename( "" )
  , m_normalization( "" )
  , ANNModelFilename( "" )
  , RandomForestModelFilename( "" )
  , ANNTestingSSEFilename( "" )
  , myConfigurationFilename( "" )
  , myConfiguration( nullptr )
{}

BRAINSCutDataHandler ::BRAINSCutDataHandler( const std::string & modelConfigurationFilename )
  : trainingVectorConfiguration( nullptr )
  , TrainConfiguration( nullptr )
  , m_atlasDataSet( nullptr )
  , m_atlasFilename( "" )
  , m_atlasBinaryFilename( "" )
  , m_atlasImage( nullptr )
  , m_roiDataList( nullptr )
  , m_roiIDsInOrder()
  , roiCount( 0 )
  , ROIRegistrationToSubject( true )
  , registrationParser( nullptr )
  , registrationImageTypeToUse( "" )
  , registrationID( "" )
  , roiAutoDilateSize( 0 )
  , m_rho( nullptr )
  , m_phi( nullptr )
  , m_theta( nullptr )
  , m_gradientSize( 0 )
  , m_trainVectorFilename( "" )
  , m_normalization( "" )
  , ANNModelFilename( "" )
  , RandomForestModelFilename( "" )
  , ANNTestingSSEFilename( "" )
  , myConfigurationFilename( "" )
  , myConfiguration( nullptr )
{
  try
  {
    SetNetConfigurationFilename( modelConfigurationFilename );
    SetNetConfiguration();
  }
  catch ( BRAINSCutExceptionStringHandler & e )
  {
    std::cout << e.Error() << std::endl;
    exit( EXIT_FAILURE );
  }
}

void
BRAINSCutDataHandler ::SetNetConfiguration()
{
  try
  {
    // read in xml images
    std::list< ElementContainer * > elementList;

    elementList.push_front( myConfiguration );

    XMLConfigurationFileParser BRAINSCutXMLConfigurationFileParser =
      XMLConfigurationFileParser( myConfigurationFilename );
    BRAINSCutXMLConfigurationFileParser.SetUserData( &elementList );
    BRAINSCutXMLConfigurationFileParser.Parse();

    myConfiguration = BRAINSCutXMLConfigurationFileParser.GetConfiguration();
  }
  catch ( BRAINSCutExceptionStringHandler & e )
  {
    std::cout << e.Error() << std::endl;
    exit( EXIT_FAILURE );
  }
}

void
BRAINSCutDataHandler ::SetAtlasDataSet()
{
  m_atlasDataSet = myConfiguration->GetAtlasDataSet();
  std::cout << "registrationImageTypeToUse :: " << registrationImageTypeToUse << std::endl;

  if ( registrationImageTypeToUse.empty() )
  {
    std::cout << "registrationImageTypeToUse is empty." << std::endl;
    exit( EXIT_FAILURE );
  }
  m_atlasFilename = m_atlasDataSet->GetImageFilenameByType( registrationImageTypeToUse );
  m_atlasBinaryFilename = m_atlasDataSet->GetMaskFilenameByType( "RegistrationROI" );
  std::cout << m_atlasBinaryFilename << std::endl;
}

void
BRAINSCutDataHandler ::SetAtlasImage()
{
  m_atlasImage = ReadImageByFilename( m_atlasFilename );
}

void
BRAINSCutDataHandler ::SetRhoPhiTheta()
{
  m_rho = ReadImageByFilename( m_atlasDataSet->GetSpatialLocationFilenameByType( "rho" ) );
  m_phi = ReadImageByFilename( m_atlasDataSet->GetSpatialLocationFilenameByType( "phi" ) );
  m_theta = ReadImageByFilename( m_atlasDataSet->GetSpatialLocationFilenameByType( "theta" ) );
}

void
BRAINSCutDataHandler ::SetNetConfigurationFilename( const std::string & filename )
{
  myConfigurationFilename = filename;
}

std::string
BRAINSCutDataHandler ::GetNetConfigurationFilename()
{
  return myConfigurationFilename;
}

int
BRAINSCutDataHandler ::GetTrainIteration()
{
  int trainIteration =
    myConfiguration->Get< TrainingParameters >( "ANNParameters" )->GetAttribute< IntValue >( "Iterations" );

  return trainIteration;
}

DataSet::StringVectorType
BRAINSCutDataHandler ::GetROIIDsInOrder() const
{
  return this->m_roiIDsInOrder;
}

void
BRAINSCutDataHandler ::SetRegionsOfInterest()
{
  m_roiDataList = myConfiguration->Get< ProbabilityMapList >( "ProbabilityMapList" );
  this->m_roiIDsInOrder = m_roiDataList->CollectAttValues< ProbabilityMapParser >( "StructureID" );

  std::sort( this->m_roiIDsInOrder.begin(), this->m_roiIDsInOrder.end() ); // get l_caudate, l_globus, .. , r_caudate,
                                                                           // r_globus..
  roiCount = m_roiDataList->size();
}

/** registration related */
void
BRAINSCutDataHandler ::SetRegistrationParameters()
{
  registrationParser = myConfiguration->Get< RegistrationConfigurationParser >( "RegistrationConfiguration" );

  SetRegistrationImageTypeToUse( std::string( registrationParser->GetAttribute< StringValue >( "ImageTypeToUse" ) ) );

  registrationID = std::string( registrationParser->GetAttribute< StringValue >( "ID" ) );

  ROIRegistrationToSubject =
    bool( registrationParser->GetAttribute< BooleanValue >( "ProbabilityMapRegistrationToSubject" ) );

  roiAutoDilateSize = registrationParser->GetAttribute< IntValue >( "BRAINSROIAutoDilateSize" );
}

void
BRAINSCutDataHandler ::SetRegistrationImageTypeToUse( std::string type )
{
  registrationImageTypeToUse = type;
}

std::string
BRAINSCutDataHandler ::GetRegistrationImageTypeToUse()
{
  if ( registrationImageTypeToUse.empty() )
  {
    std::cout << "registrationImageTypeToUse is empty." << std::endl;
    exit( EXIT_FAILURE );
  }
  return registrationImageTypeToUse;
}

std::string
BRAINSCutDataHandler ::GetRegistrationID()
{
  if ( registrationID.empty() )
  {
    std::cout << "The registrationID is empty" << std::endl;
    exit( EXIT_FAILURE );
  }
  return registrationID;
}

std::string
BRAINSCutDataHandler ::GetSubjectToAtlasRegistrationFilename( DataSet & subject )
{
  std::string filename =
    subject.GetRegistrationWithID( registrationID )->GetAttribute< StringValue >( "SubjToAtlasRegistrationFilename" );

  return filename;
}

std::string
BRAINSCutDataHandler ::GetAtlasToSubjectRegistrationFilename( DataSet & subject )
{
  std::string filename =
    subject.GetRegistrationWithID( registrationID )->GetAttribute< StringValue >( "AtlasToSubjRegistrationFilename" );

  return filename;
}

void
BRAINSCutDataHandler ::GetDeformedSpatialLocationImages(
  std::map< std::string, WorkingImagePointer > & warpedSpatialLocationImages, DataSet & subject )
{
  using GenericTransformType = itk::Transform< double, 3, 3 >;
  std::string atlasSubjectRegistrationFilename = GetAtlasToSubjectRegistrationFilename( subject );

  DisplacementFieldType::Pointer deformation = GetDeformationField( atlasSubjectRegistrationFilename );
  GenericTransformType::Pointer  genericTransform = GetGenericTransform( atlasSubjectRegistrationFilename );

  // GenericTransformImage's behavior was OR not and -- if both a
  // deformation and a genericTransform were non-null then the generic
  // transform was ignored.
  if ( deformation.IsNotNull() )
  {
    using DisplacementFieldTransformType =
      itk::DisplacementFieldTransform< DeformationScalarType, DisplacementFieldType::ImageDimension >;
    DisplacementFieldTransformType::Pointer dispXfrm = DisplacementFieldTransformType::New();
    dispXfrm->SetDisplacementField( deformation );
    genericTransform = dispXfrm.GetPointer();
  }
  std::string subjectFilenameToUse = subject.GetImageFilenameByType( registrationImageTypeToUse );

  if ( !itksys::SystemTools::FileExists( subjectFilenameToUse.c_str(), false ) )
  {
    std::cout << " Subject image file of " << subjectFilenameToUse << ", type of " << registrationImageTypeToUse
              << ", does not exist. " << std::endl;
    exit( EXIT_FAILURE );
  }

  WorkingImagePointer referenceImage =
    ReadImageByFilename( subject.GetImageFilenameByType( registrationImageTypeToUse ) );

  const std::string transoformationPixelType = "float";

  warpedSpatialLocationImages.insert( std::pair< std::string, WorkingImagePointer >(
    "rho",
    GenericTransformImage< WorkingImageType, WorkingImageType, DisplacementFieldType >(
      m_rho, referenceImage, genericTransform.GetPointer(), 0.0, "Linear", transoformationPixelType == "binary" ) ) );
  warpedSpatialLocationImages.insert( std::pair< std::string, WorkingImagePointer >(
    "phi",
    GenericTransformImage< WorkingImageType, WorkingImageType, DisplacementFieldType >(
      m_phi, referenceImage, genericTransform.GetPointer(), 0.0, "Linear", transoformationPixelType == "binary" ) ) );
  warpedSpatialLocationImages.insert( std::pair< std::string, WorkingImagePointer >(
    "theta",
    GenericTransformImage< WorkingImageType, WorkingImageType, DisplacementFieldType >(
      m_theta, referenceImage, genericTransform.GetPointer(), 0.0, "Linear", transoformationPixelType == "binary" ) ) );
}

WorkingImagePointer
BRAINSCutDataHandler ::GetCandidateRegion( DataSet & subject ) const
{
  const std::string candidateRegionFilename = subject.GetImageFilenameByType( "candidateRegion" );

  if ( candidateRegionFilename == std::string( "" ) )
  {
    std::cout << "* No candidate region is given! " << std::endl;
    return nullptr;
  }
  else if ( !itksys::SystemTools::FileExists( candidateRegionFilename.c_str() ) )
  {
    std::cout << "Requested candidateRegion file does not exists!" << std::endl
              << ":: " << candidateRegionFilename << std::endl;
    std::exit( EXIT_FAILURE );
  }
  WorkingImagePointer candidateRegion = ReadImageByFilename( candidateRegionFilename );
  return candidateRegion;
}

void
BRAINSCutDataHandler ::ReadImagesOfSubjectInOrder( WorkingImageVectorType & subjectImageList, DataSet & subject )
{
  DataSet::StringVectorType imageListFromAtlas = m_atlasDataSet->GetImageTypes(); // T1, T2, SG, ...

  std::sort( imageListFromAtlas.begin(), imageListFromAtlas.end() ); // SG, T1, T2, ... ascending order

  for ( DataSet::StringVectorType::iterator imgTyIt = imageListFromAtlas.begin(); imgTyIt != imageListFromAtlas.end();
        ++imgTyIt ) // imgTyIt = image type iterator
  {
    std::cout << *imgTyIt << std::endl;
    WorkingImagePointer currentTypeImage = ReadImageByFilename( subject.GetImageFilenameByType( *imgTyIt ) );
    subjectImageList.push_back( currentTypeImage );
  }
}

void
BRAINSCutDataHandler ::GetDeformedROIs( std::map< std::string, WorkingImagePointer > & warpedROIs, DataSet & subject )
{
  using GenericTransformType = itk::Transform< double, 3, 3 >;
  std::string atlasSubjectRegistrationFilename = GetAtlasToSubjectRegistrationFilename( subject );

  /** Get the transformation file
   * Note that only one of transformation type will be used. Either deformation or transformation
   * That determined based on the file name at the GetDeformationField
   */
  DisplacementFieldType::Pointer deformation = GetDeformationField( atlasSubjectRegistrationFilename );
  GenericTransformType::Pointer  genericTransform = GetGenericTransform( atlasSubjectRegistrationFilename );

  // GenericTransformImage's behavior was OR not and -- if both a
  // deformation and a genericTransform were non-null then the generic
  // transform was ignored.
  if ( deformation.IsNotNull() )
  {
    using DisplacementFieldTransformType =
      itk::DisplacementFieldTransform< DeformationScalarType, DisplacementFieldType::ImageDimension >;
    DisplacementFieldTransformType::Pointer dispXfrm = DisplacementFieldTransformType::New();
    dispXfrm->SetDisplacementField( deformation );
    genericTransform = dispXfrm.GetPointer();
  }

  if ( !ROIRegistrationToSubject )
  {
    typedef itk::IdentityTransform< DeformationScalarType, DisplacementFieldType::ImageDimension >
                                   IdentityTransformType;
    IdentityTransformType::Pointer identityXfrm = IdentityTransformType::New();
    genericTransform = identityXfrm.GetPointer();
  }

  WorkingImagePointer referenceImage =
    ReadImageByFilename( subject.GetImageFilenameByType( registrationImageTypeToUse ) );

  const std::string transformationPixelType = "float";

  for ( DataSet::StringVectorType::iterator roiTyIt = this->m_roiIDsInOrder.begin();
        roiTyIt != this->m_roiIDsInOrder.end();
        ++roiTyIt )
  {
    std::string roiFilename = m_roiDataList->GetMatching< ProbabilityMapParser >( "StructureID", ( *roiTyIt ).c_str() )
                                ->GetAttribute< StringValue >( "Filename" );
    WorkingImagePointer currentROI = ReadImageByFilename( roiFilename );

    warpedROIs.insert( std::pair< std::string, WorkingImagePointer >(
      ( *roiTyIt ),
      GenericTransformImage< WorkingImageType, WorkingImageType, DisplacementFieldType >( currentROI,
                                                                                          referenceImage,
                                                                                          genericTransform.GetPointer(),
                                                                                          0.0,
                                                                                          "Linear",
                                                                                          transformationPixelType ==
                                                                                            "binary" ) ) );
  }
}

void
BRAINSCutDataHandler ::SetTrainingVectorConfiguration()
{
  trainingVectorConfiguration =
    myConfiguration->Get< TrainingVectorConfigurationType >( "TrainingVectorConfiguration" );
}

void
BRAINSCutDataHandler ::SetGradientSize()
{
  m_gradientSize = trainingVectorConfiguration->GetAttribute< IntValue >( "GradientProfileSize" );
}

unsigned int
BRAINSCutDataHandler ::GetGradientSize()
{
  return m_gradientSize;
}

void
BRAINSCutDataHandler ::SetNormalization()
{
  std::string normalizationString;

  try
  {
    normalizationString = trainingVectorConfiguration->GetAttribute< StringValue >( "Normalization" );
  }
  catch ( BRAINSCutExceptionStringHandler & e )
  {
    std::cout << e.Error() << std::endl;
    exit( EXIT_FAILURE );
  }

  if ( normalizationString == "true" )
  {
    m_normalization = "Linear";
  }
  else if ( normalizationString == "false" )
  {
    m_normalization = "None";
  }
  else
  {
    m_normalization = normalizationString;
  }
  std::cout << "Get Normalization from XML file --> " << m_normalization << std::endl;
}

std::string
BRAINSCutDataHandler ::GetNormalizationMethod()
{
  return m_normalization;
}

/** model file name */
std::string
BRAINSCutDataHandler ::GetModelBaseName()
{
  std::string basename;

  try
  {
    basename = trainingVectorConfiguration->GetAttribute< StringValue >( "TrainingModelFilename" );
  }
  catch ( ... )
  {
    std::cout << "Fail to get the ann model file name:: " << basename << std::endl;
    throw BRAINSCutExceptionStringHandler( "Fail to get the ann model file name" );
  }
  return basename;
}

std::string
BRAINSCutDataHandler ::GetANNModelFilename()
{
  return ANNModelFilename;
}

std::string
BRAINSCutDataHandler ::GetANNModelFilenameAtIteration( const int iteration )
{
  SetANNModelFilenameAtIteration( iteration );
  return ANNModelFilename;
}

void
BRAINSCutDataHandler ::SetANNModelFilenameAtIteration( const int iteration )
{
  ANNModelFilename = GetModelBaseName();

  char temp[10];
  sprintf( temp, "%09d", iteration );
  ANNModelFilename += temp;
}

std::string
BRAINSCutDataHandler ::GetRFModelFilename( unsigned int depth, unsigned int NTrees )
{
  std::string basename = GetModelBaseName();

  char tempDepth[5];

  sprintf( tempDepth, "%04u", depth );

  char tempNTrees[5];
  sprintf( tempNTrees, "%04u", NTrees );

  std::string filename = basename + "D" + tempDepth + "NT" + tempNTrees + ".gz";

  return filename;
}

std::string
BRAINSCutDataHandler ::GetAtlasFilename()
{
  return m_atlasFilename;
}

std::string
BRAINSCutDataHandler ::GetAtlasBinaryFilename()
{
  return m_atlasBinaryFilename;
}

int
BRAINSCutDataHandler ::GetROIAutoDilateSize()
{
  return roiAutoDilateSize;
}

unsigned int
BRAINSCutDataHandler ::GetROICount()
{
  return roiCount;
}

WorkingImagePointer
BRAINSCutDataHandler ::GetAtlasImage()
{
  return m_atlasImage;
}

ProbabilityMapList *
BRAINSCutDataHandler ::GetROIDataList()
{
  return m_roiDataList;
}

DataSet *
BRAINSCutDataHandler ::GetAtlasDataSet()
{
  return m_atlasDataSet;
}

scalarType
BRAINSCutDataHandler ::GetANNOutputThreshold()
{
  scalarType annOutputThreshold =
    myConfiguration->Get< ApplyModelType >( "ApplyModel" )->GetAttribute< FloatValue >( "MaskThresh" );

  if ( annOutputThreshold < 0.0F )
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
BRAINSCutDataHandler ::SetRandomForestModelFilename( std::string name )
{
  RandomForestModelFilename = name;
}

void
BRAINSCutDataHandler ::SetRandomForestModelFilename( int depth, int nTree )
{
  if ( depth < 0 && nTree < 0 )
  {
    nTree =
      myConfiguration->Get< TrainingParameters >( "RandomForestParameters" )->GetAttribute< IntValue >( "MaxDepth" );
    depth = myConfiguration->Get< TrainingParameters >( "RandomForestParameters" )
              ->GetAttribute< IntValue >( "MaxTreeCount" );
  }
  RandomForestModelFilename = GetRFModelFilename( depth, nTree );
}

std::string
BRAINSCutDataHandler ::GetRandomForestModelFilename()
{
  return RandomForestModelFilename;
}

void
BRAINSCutDataHandler ::SetANNTestingSSEFilename()
{
  ANNTestingSSEFilename = GetModelBaseName();
  ANNTestingSSEFilename += "ValidationSetSSE.txt";
}

std::string
BRAINSCutDataHandler ::GetANNTestingSSEFilename()
{
  return ANNTestingSSEFilename;
}

void
BRAINSCutDataHandler ::SetTrainVectorFilename()
{
  m_trainVectorFilename = trainingVectorConfiguration->GetAttribute< StringValue >( "TrainingVectorFilename" );
  m_trainVectorFilename += "ANN"; // TODO
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "vector file at " << m_trainVectorFilename << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
}

std::string
BRAINSCutDataHandler ::GetTrainVectorFilename()
{
  if ( m_trainVectorFilename.empty() )
  {
    std::string msg = "The train vector file name is empty.\n";
    throw BRAINSCutExceptionStringHandler( msg );
  }
  return m_trainVectorFilename;
}

BRAINSCutConfiguration::TrainDataSetListType
BRAINSCutDataHandler ::GetTrainDataSet()
{
  BRAINSCutConfiguration::TrainDataSetListType trainDataSetList;

  try
  {
    trainDataSetList = myConfiguration->GetTrainDataSets();
  }
  catch ( BRAINSCutExceptionStringHandler & e )
  {
    std::cout << e.Error() << std::endl;
    exit( EXIT_SUCCESS );
  }
  return trainDataSetList;
}

BRAINSCutConfiguration::ApplyDataSetListType
BRAINSCutDataHandler ::GetApplyDataSet()
{
  BRAINSCutConfiguration::ApplyDataSetListType applyDataSetList;

  try
  {
    applyDataSetList = myConfiguration->GetApplyDataSets();
  }
  catch ( BRAINSCutExceptionStringHandler & e )
  {
    std::cout << e.Error() << std::endl;
    exit( EXIT_SUCCESS );
  }
  return applyDataSetList;
}

scalarType
BRAINSCutDataHandler ::GetGaussianSmoothingSigma()
{
  scalarType gaussianSmoothingSigma =
    myConfiguration->Get< ApplyModelType >( "ApplyModel" )->GetAttribute< FloatValue >( "GaussianSmoothingSigma" );

  return gaussianSmoothingSigma;
}

void
BRAINSCutDataHandler ::SetTrainConfiguration( std::string trainParameterName )
{
  TrainConfiguration = myConfiguration->Get< TrainingParameters >( trainParameterName.c_str() );
}

unsigned int
BRAINSCutDataHandler ::GetEpochIteration()
{
  unsigned int trainEpochIteration = TrainConfiguration->GetAttribute< IntValue >( "EpochIterations" );

  return trainEpochIteration;
}

float
BRAINSCutDataHandler ::GetDesiredError()
{
  float trainDesiredError = TrainConfiguration->GetAttribute< FloatValue >( "DesiredError" );

  return trainDesiredError;
}

unsigned int
BRAINSCutDataHandler ::GetMaximumDataSize()
{
  unsigned int trainMaximumDataSize = TrainConfiguration->GetAttribute< IntValue >( "MaximumVectorsPerEpoch" );

  return trainMaximumDataSize;
}

int
BRAINSCutDataHandler ::GetANNHiddenNodesNumber()
{
  int ANNHiddenNodesNumber = TrainConfiguration->GetAttribute< IntValue >( "NumberOfHiddenNodes" );

  return ANNHiddenNodesNumber;
}

float
BRAINSCutDataHandler ::GetActivationFunctionSlope()
{
  float activationSlope = TrainConfiguration->GetAttribute< FloatValue >( "ActivationSlope" );

  return activationSlope;
}

float
BRAINSCutDataHandler ::GetActivationFunctionMinMax()
{
  float activationMinMax = TrainConfiguration->GetAttribute< FloatValue >( "ActivationMinMax" );

  return activationMinMax;
}

/** Get Random Tree */
int
BRAINSCutDataHandler ::GetMaxDepth()
{
  int trainMaxDepth = TrainConfiguration->GetAttribute< IntValue >( "MaxDepth" );

  return trainMaxDepth;
}

int
BRAINSCutDataHandler ::GetMinSampleCount()
{
  int trainMinSampleCount = TrainConfiguration->GetAttribute< IntValue >( "MinSampleCount" );

  return trainMinSampleCount;
}

bool
BRAINSCutDataHandler ::GetUseSurrogates()
{
  bool trainUseSurrogates = TrainConfiguration->GetAttribute< BooleanValue >( "UseSurrogates" );

  return trainUseSurrogates;
}

bool
BRAINSCutDataHandler ::GetCalcVarImportance()
{
  bool trainCalcVarImportance = TrainConfiguration->GetAttribute< BooleanValue >( "CalcVarImportance" );

  return trainCalcVarImportance;
}

int
BRAINSCutDataHandler ::GetMaxTreeCount()
{
  int trainMaxTreeCount = TrainConfiguration->GetAttribute< IntValue >( "MaxTreeCount" );

  return trainMaxTreeCount;
}
