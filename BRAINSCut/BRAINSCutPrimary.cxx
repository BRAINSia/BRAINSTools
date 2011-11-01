#include "BRAINSCutPrimary.h"
#include "NetConfigurationParser.h"

#include "GenericTransformImage.h"

/** constructors */
BRAINSCutPrimary
::BRAINSCutPrimary( std::string netConfigurationFilename )
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  SetNetConfigurationFilename( netConfigurationFilename );
  SetNetConfiguration();
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
}

void
BRAINSCutPrimary
::SetNetConfiguration()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  std::list<XMLElementContainer *> elementList;

  elementList.push_front( &BRAINSCutNetConfiguration );

  NetConfigurationParser BRIANSCutNetConfigurationParser = NetConfigurationParser( NetConfigurationFilename );
  BRIANSCutNetConfigurationParser.SetUserData( &elementList );
  BRIANSCutNetConfigurationParser.Parse();
}
}

void
BRAINSCutPrimary
::SetAtlasDataSet()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  atlasDataSet = BRAINSCutNetConfiguration.GetAtlasDataSet();
}

void
BRAINSCutPrimary
::SetAtlasImage()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  atlasImage = ReadImageByFilename( atlasDataSet->GetImageFilenameByType( registrationImageTypeToUse) );
}

void
BRAINSCutPrimary
::SetRhoPhiThetaFromNetConfiguration()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  rho = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("rho") );
  phi = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("phi") );
  theta = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("theta") );
}

void
BRAINSCutPrimary
::SetNetConfigurationFilename(const std::string filename)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  NetConfigurationFilename = filename;
}

std::string
BRAINSCutPrimary
::GetNetConfigurationFilename()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  return NetConfigurationFilename;
}

DataSet::StringVectorType
BRAINSCutPrimary
::GetROIIDsInOrder()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  return roiIDsInOrder;
}

void
BRAINSCutPrimary
::SetRegionsOfInterestFromNetConfiguration()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  roiDataList = BRAINSCutNetConfiguration.Get<ProbabilityMapList>("ProbabilityMapList");
  roiIDsInOrder = roiDataList->CollectAttValues<ProbabilityMapParser>("StructureID");

  std::sort( roiIDsInOrder.begin(), roiIDsInOrder.end() ); // get l_caudate, l_globus, .. , r_caudate, r_globus..
  roiCount = roiDataList->size();
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
}

void
BRAINSCutPrimary
::SetRegistrationParametersFromNetConfiguration()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  registrationParser =
    BRAINSCutNetConfiguration.Get<RegistrationConfigurationParser>("RegistrationConfiguration");

  registrationImageTypeToUse =
    std::string( registrationParser->GetAttribute<StringValue>( "ImageTypeToUse") );
  registrationID = std::string(
      registrationParser->GetAttribute<StringValue>("ID") );
}

inline
std::string
BRAINSCutPrimary
::GetAtlasToSubjectRegistrationFilename( DataSet& subject)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  std::string filename = subject.GetRegistrationWithID( registrationID )
    ->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename");
  std::cout << __LINE__ << "::" << __FILE__ << "::" << filename << std::endl;

  return filename;
}

void
BRAINSCutPrimary
::GetDeformedSpatialLocationImages( std::map<std::string, WorkingImagePointer>& warpedSpatialLocationImages,
                                    DataSet& subject)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  std::string atlasSubjectRegistrationFilename = GetAtlasToSubjectRegistrationFilename( subject );

  DeformationFieldType::Pointer deformation = GetDeformationField( atlasSubjectRegistrationFilename );
  GenericTransformType::Pointer genericTransform = GetGenericTransform( atlasSubjectRegistrationFilename );

  WorkingImagePointer referenceImage =
    ReadImageByFilename( subject.GetImageFilenameByType(registrationImageTypeToUse) );
  const std::string transoformationPixelType = "float";

  warpedSpatialLocationImages.insert( pair<std::string, WorkingImagePointer>
                                        ("rho", GenericTransformImage<WorkingImageType,
                                                                      WorkingImageType,
                                                                      DeformationFieldType>
                                          ( rho, referenceImage, deformation, genericTransform,
                                          0.0, "Linear", transoformationPixelType == "binary") ) );
  warpedSpatialLocationImages.insert( pair<std::string, WorkingImagePointer>
                                        ("phi", GenericTransformImage<WorkingImageType,
                                                                      WorkingImageType,
                                                                      DeformationFieldType>
                                          ( phi, referenceImage, deformation, genericTransform,
                                          0.0, "Linear", transoformationPixelType == "binary") ) );
  warpedSpatialLocationImages.insert( pair<std::string, WorkingImagePointer>
                                        ("theta", GenericTransformImage<WorkingImageType,
                                                                        WorkingImageType,
                                                                        DeformationFieldType>
                                          ( theta, referenceImage, deformation, genericTransform,
                                          0.0, "Linear", transoformationPixelType == "binary") ) );
}

void
BRAINSCutPrimary
::GetImagesOfSubjectInOrder( WorkingImageVectorType& subjectImageList, DataSet& subject)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  DataSet::StringVectorType imageListFromAtlas = atlasDataSet->GetImageTypes(); // T1, T2, SG, ...
  std::sort( imageListFromAtlas.begin(), imageListFromAtlas.end() );            // SG, T1, T2, ... ascending order

  for( DataSet::StringVectorType::iterator imgTyIt = imageListFromAtlas.begin();
       imgTyIt != imageListFromAtlas.end();
       ++imgTyIt ) // imgTyIt = image type iterator
    {
    WorkingImagePointer currentTypeImage = ReadImageByFilename( subject.GetImageFilenameByType( *imgTyIt ) );
    subjectImageList.push_back( currentTypeImage );
    }
}

void
BRAINSCutPrimary
::GetDeformedROIs( std::map<std::string, WorkingImagePointer>& warpedROIs,
                   DataSet& subject)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  std::string atlasSubjectRegistrationFilename = GetAtlasToSubjectRegistrationFilename( subject );

  DeformationFieldType::Pointer deformation = GetDeformationField( atlasSubjectRegistrationFilename );
  GenericTransformType::Pointer genericTransform = GetGenericTransform( atlasSubjectRegistrationFilename );

  WorkingImagePointer referenceImage =
    ReadImageByFilename( subject.GetImageFilenameByType(registrationImageTypeToUse) );

  const std::string transformationPixelType = "float";

  for( DataSet::StringVectorType::iterator roiTyIt = roiIDsInOrder.begin();
       roiTyIt != roiIDsInOrder.end();
       ++roiTyIt )
    {
    std::string roiFilename = roiDataList->GetMatching<ProbabilityMapParser>( "StructureID", (*roiTyIt).c_str() )
      ->GetAttribute<StringValue>("Filename");
    WorkingImagePointer currentROI = ReadImageByFilename( roiFilename );

    std::cout << __LINE__ << "::" << __FILE__ << std::endl
              << "Read :: " << roiFilename << std::endl;

    warpedROIs.insert( pair<std::string, WorkingImagePointer>(
                         (*roiTyIt), GenericTransformImage<WorkingImageType,
                                                           WorkingImageType,
                                                           DeformationFieldType>
                           ( currentROI, referenceImage, deformation,
                           genericTransform, 0.0, "Linear",
                           transformationPixelType == "binary") ) );
    }
}

void
BRAINSCutPrimary
::SetANNModelConfiguration()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  annModelConfiguration = BRAINSCutNetConfiguration.Get<NeuralParams>("NeuralNetParams");
}

void
BRAINSCutPrimary
::SetGradientSizeFromNetConfiguration()
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  gradientSize = annModelConfiguration->GetAttribute<IntValue>("GradientProfileSize");
}

/* inline functions */

inline WorkingImagePointer
BRAINSCutPrimary
::ReadImageByFilename( const std::string  filename )
{
  std::cout << __LINE__ << "::" << __FILE__ << "::" << filename << std::endl;
  WorkingImagePointer readInImage;

  ReadInImagePointer inputImage = itkUtil::ReadImage<ReadInImageType>(filename.c_str() );

  readInImage = itkUtil::ScaleAndCast<ReadInImageType,
                                      WorkingImageType>(inputImage,
                                                        ZeroPercentValue,
                                                        HundreadPercentValue);
  return readInImage;
}

inline
DeformationFieldType::Pointer
BRAINSCutPrimary
::GetDeformationField( std::string filename)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  const bool useTransform( filename.find(".mat") != std::string::npos );

  if( useTransform )
    {
    return NULL;
    }
  typedef itk::ImageFileReader<DeformationFieldType> DeformationReaderType;
  DeformationReaderType::Pointer deformationReader = DeformationReaderType::New();
  deformationReader->SetFileName( filename );
  deformationReader->Update();

  return deformationReader->GetOutput();
}

inline
GenericTransformType::Pointer
BRAINSCutPrimary
::GetGenericTransform( std::string filename)
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  const bool useDeformation( filename.find(".mat") == std::string::npos );

  if( useDeformation )
    {
    return NULL;
    }
  return itk::ReadTransformFromDisk( filename );
}

bool
BRAINSCutPrimary
::GetNormalizationFromNetConfiguration()
{
  std::string normalizationString = annModelConfiguration->GetAttribute<StringValue>("Normalization");

  if( normalizationString == "true" )
    {
    return true;
    }
  else
    {
    return false;
    }
}
