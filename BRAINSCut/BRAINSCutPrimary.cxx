#include "BRAINSCutPrimary.h"
#include "NetConfigurationParser.h"

#include "GenericTransformImage.h"

#include <itkSmoothingRecursiveGaussianImageFilter.h>
/** constructors */
BRAINSCutPrimary
::BRAINSCutPrimary( std::string netConfigurationFilename )
{
  SetNetConfigurationFilename( netConfigurationFilename );
  SetNetConfiguration();
}

void
BRAINSCutPrimary
::SetNetConfiguration()
{
  std::list<XMLElementContainer *> elementList;

  elementList.push_front( &BRAINSCutNetConfiguration );

  NetConfigurationParser BRIANSCutNetConfigurationParser = NetConfigurationParser( NetConfigurationFilename );
  BRIANSCutNetConfigurationParser.SetUserData( &elementList );
  BRIANSCutNetConfigurationParser.Parse();
}

void
BRAINSCutPrimary
::SetAtlasDataSet()
{
  atlasDataSet = BRAINSCutNetConfiguration.GetAtlasDataSet();
}

void
BRAINSCutPrimary
::SetAtlasFilename()
{
  atlasFilename = atlasDataSet->GetImageFilenameByType( registrationImageTypeToUse);
  atlasBinaryFilename = atlasDataSet->GetMaskFilenameByType( "RegistrationROI" );
  std::cout << atlasBinaryFilename << std::endl;
}

void
BRAINSCutPrimary
::SetAtlasImage()
{
  atlasImage = ReadImageByFilename( atlasFilename );
}

void
BRAINSCutPrimary
::SetRhoPhiThetaFromNetConfiguration()
{
  rho = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("rho") );
  phi = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("phi") );
  theta = ReadImageByFilename( atlasDataSet->GetSpatialLocationFilenameByType("theta") );
}

void
BRAINSCutPrimary
::SetNetConfigurationFilename(const std::string filename)
{
  NetConfigurationFilename = filename;
}

std::string
BRAINSCutPrimary
::GetNetConfigurationFilename()
{
  return NetConfigurationFilename;
}

DataSet::StringVectorType
BRAINSCutPrimary
::GetROIIDsInOrder()
{
  return roiIDsInOrder;
}

void
BRAINSCutPrimary
::SetRegionsOfInterestFromNetConfiguration()
{
  roiDataList = BRAINSCutNetConfiguration.Get<ProbabilityMapList>("ProbabilityMapList");
  roiIDsInOrder = roiDataList->CollectAttValues<ProbabilityMapParser>("StructureID");

  std::sort( roiIDsInOrder.begin(), roiIDsInOrder.end() ); // get l_caudate, l_globus, .. , r_caudate, r_globus..
  roiCount = roiDataList->size();
}

/** registration related */
void
BRAINSCutPrimary
::SetRegistrationParametersFromNetConfiguration()
{
  registrationParser =
    BRAINSCutNetConfiguration.Get<RegistrationConfigurationParser>("RegistrationConfiguration");

  registrationImageTypeToUse =
    std::string( registrationParser->GetAttribute<StringValue>( "ImageTypeToUse") );

  registrationID = std::string(
      registrationParser->GetAttribute<StringValue>("ID") );

  roiAutoDilateSize = registrationParser->GetAttribute<IntValue>("BRAINSROIAutoDilateSize");
}

std::string
BRAINSCutPrimary
::GetSubjectToAtlasRegistrationFilename( DataSet& subject)
{
  std::string filename = subject.GetRegistrationWithID( registrationID )
    ->GetAttribute<StringValue>("SubjToAtlasRegistrationFilename");

  return filename;
}

std::string
BRAINSCutPrimary
::GetAtlasToSubjectRegistrationFilename( DataSet& subject)
{
  std::string filename = subject.GetRegistrationWithID( registrationID )
    ->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename");

  return filename;
}

void
BRAINSCutPrimary
::GetDeformedSpatialLocationImages( std::map<std::string, WorkingImagePointer>& warpedSpatialLocationImages,
                                    DataSet& subject)
{
  std::string atlasSubjectRegistrationFilename = GetAtlasToSubjectRegistrationFilename( subject );

  DisplacementFieldType::Pointer deformation = GetDeformationField( atlasSubjectRegistrationFilename );
  GenericTransformType::Pointer  genericTransform = GetGenericTransform( atlasSubjectRegistrationFilename );

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
BRAINSCutPrimary
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
BRAINSCutPrimary
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

  for( DataSet::StringVectorType::iterator roiTyIt = roiIDsInOrder.begin();
       roiTyIt != roiIDsInOrder.end();
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
BRAINSCutPrimary
::SetANNModelConfiguration()
{
  annModelConfiguration = BRAINSCutNetConfiguration.Get<NeuralParams>("NeuralNetParams");
}

void
BRAINSCutPrimary
::SetGradientSizeFromNetConfiguration()
{
  gradientSize = annModelConfiguration->GetAttribute<IntValue>("GradientProfileSize");
}

WorkingImagePointer
BRAINSCutPrimary
::ReadImageByFilename( const std::string  filename )
{
  WorkingImagePointer readInImage;

  ReadInImagePointer inputImage = itkUtil::ReadImage<ReadInImageType>(filename.c_str() );

  readInImage = itkUtil::ScaleAndCast<ReadInImageType,
                                      WorkingImageType>(inputImage,
                                                        ZeroPercentValue,
                                                        HundreadPercentValue);
  return readInImage;
}

/* inline functions */

inline
DisplacementFieldType::Pointer
BRAINSCutPrimary
::GetDeformationField( std::string filename)
{
  const bool useTransform( filename.find(".mat") != std::string::npos );

  if( useTransform )
    {
    return NULL;
    }
  typedef itk::ImageFileReader<DisplacementFieldType> DeformationReaderType;
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
  std::string normalizationString;

  try
    {
    normalizationString = annModelConfiguration->GetAttribute<StringValue>("Normalization");
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_FAILURE);
    }

  if( normalizationString == "true" )
    {
    return true;
    }
  else
    {
    return false;
    }
}

WorkingImagePointer
BRAINSCutPrimary
::SmoothImage( const WorkingImagePointer image, const float GaussianValue)
{
  if( GaussianValue < 0 + FLOAT_TOLERANCE )
    {
    std::cout << "Gaussian value is less than tolerance. "
              << "No smoothing occurs at this time"
              << std::endl;
    return image;
    }
  std::cout << "Smooth Image with Gaussian value of :: "
            << GaussianValue
            << std::endl;
  typedef itk::SmoothingRecursiveGaussianImageFilter<WorkingImageType, WorkingImageType> SmoothingFilterType;
  SmoothingFilterType::Pointer smoothingFilter = SmoothingFilterType::New();

  smoothingFilter->SetInput( image);
  smoothingFilter->SetSigma( GaussianValue );

  smoothingFilter->Update();

  return smoothingFilter->GetOutput();
}
