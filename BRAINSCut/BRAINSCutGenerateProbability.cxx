#include "BRAINSCutGenerateProbability.h"
#include "NetConfigurationParser.h"
#include "NetConfiguration.h"

#include "itkIO.h"

/** constructors */
BRAINSCutGenerateProbability
::BRAINSCutGenerateProbability( std::string netConfigurationFilename ) : BRAINSCutPrimary( netConfigurationFilename )
{
  SetRegistrationParametersFromNetConfiguration();

  SetAtlasDataSet();
  SetAtlasFilename();
  SetAtlasImage();

  SetRegionsOfInterestFromNetConfiguration();
  SetTrainingDataSetsList();
}

/** Get/Set Methods */

void
BRAINSCutGenerateProbability
::SetTrainingDataSetsList()
{
  trainingDataSetList = BRAINSCutNetConfiguration.GetTrainDataSets();
}

/*
 * generate probability maps
 */
inline WorkingImageType::Pointer
AddImageToAccumulator( WorkingImageType::Pointer & image, WorkingImageType::Pointer & accumulator )
{
  typedef itk::BinaryThresholdImageFilter<WorkingImageType, WorkingImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();

  threshold->SetInput( image );
  threshold->SetInsideValue( 1.0F);
  threshold->SetOutsideValue( 0.0F );
  threshold->SetLowerThreshold( 0.1F);
  threshold->Update();

  typedef itk::AddImageFilter<WorkingImageType, WorkingImageType, WorkingImageType> AddType;
  AddType::Pointer adder = AddType::New();

  adder->SetInput1( threshold->GetOutput() );
  adder->SetInput2(accumulator);
  adder->Update();

  return adder->GetOutput();
}

void
BRAINSCutGenerateProbability
::GenerateProbabilityMaps()
{
  // GenerateRegistrations(BRAINSCutNetConfiguration, false, false,  1);

  /** generating spherical coordinate image does not have to be here */
  GenerateSymmetricalSphericalCoordinateImage();
  /** iterate through the rois*/
  for( unsigned int currentROIAt = 0; currentROIAt < roiCount; currentROIAt++ )
    {
    WorkingImageType::Pointer currentAccumulatedImages;
    CreateNewFloatImageFromTemplate( currentAccumulatedImages, atlasImage );

    std::string  currentROIID(    roiIDsInOrder[currentROIAt] );
    unsigned int currentROISubjectsCounter = 0;
    /** iterate through subject */
    for( std::list<DataSet *>::iterator currentSubjectIt = trainingDataSetList.begin();
         currentSubjectIt != trainingDataSetList.end();
         currentSubjectIt++ )
      {
      currentROISubjectsCounter++;
      /** deform ROI to Atlas */

      std::string currentRegistrationFilename = GetSubjectToAtlasRegistrationFilename( *(*currentSubjectIt) );
      /*std::string currentRegistrationFilename =
        (*currentSubjectIt)->GetRegistrationWithID( registrationID )
        ->GetAttribute<StringValue>( "SubjToAtlasRegistrationFIlename");*/

      std::string currentROIFilename = ( *currentSubjectIt)
        ->GetMaskFilenameByType( currentROIID );

      WorkingImageType::Pointer currentDeformedROI = ImageWarper<WorkingImageType>( currentRegistrationFilename,
                                                                                    currentROIFilename,
                                                                                    atlasImage );
      /** add the deformed roi to accumulator */
      currentAccumulatedImages =
        AddImageToAccumulator( currentDeformedROI, currentAccumulatedImages );
      } /** end of iteration for subject */

    /** average the accumulator based on the counts */
    WorkingImagePointer currentProbabilityImage =
      ImageMultiplyConstant<WorkingImageType>( currentAccumulatedImages,
                                               1.0F / static_cast<float>( currentROISubjectsCounter ) );

    /** get roi object */

    ProbabilityMapParser *currentROISet = roiDataList->GetMatching<ProbabilityMapParser>(
        "StructureID", currentROIID.c_str() );
    /** smooth the accumulated image */
    float               GaussianSigma = currentROISet->GetAttribute<FloatValue>("Gaussian");
    WorkingImagePointer currentSmoothProbabilityImage = SmoothImage( currentProbabilityImage, GaussianSigma );

    /** get filename */
    std::string currentProbabilityMapFilename( currentROISet->GetAttribute<StringValue>("Filename") );

    /** check the directory */
    std::string path = itksys::SystemTools::GetFilenamePath( currentProbabilityMapFilename );
    if( !itksys::SystemTools::FileExists( path.c_str(), false ) )
      {
      std::cout << " Probability map directory does not exist. Create as following:: "
                << path.c_str()
                << std::endl;
      itksys::SystemTools::MakeDirectory( path.c_str() );
      }

    /** write image */
    itkUtil::WriteImage<WorkingImageType>( currentSmoothProbabilityImage, currentProbabilityMapFilename );
    } /** end of iteration for roi */
}

inline WorkingImageType::IndexType::IndexValueType
TruncatedHalf(const WorkingImageType::SizeType::SizeValueType & v)
{
  return static_cast<WorkingImageType::IndexType::IndexValueType>
         ( static_cast<double>( v ) * 0.5 );
}

void
BRAINSCutGenerateProbability
::GenerateSymmetricalSphericalCoordinateImage()
{
  const WorkingImageType::SizeType atlasImageSize( atlasImage->GetLargestPossibleRegion().GetSize() );

  WorkingImageType::IndexType centerOfAtlas;

  centerOfAtlas[0] = TruncatedHalf(atlasImageSize[0]);
  centerOfAtlas[1] = TruncatedHalf(atlasImageSize[1]);
  centerOfAtlas[2] = TruncatedHalf(atlasImageSize[2]);

  itk::Point<WorkingPixelType, DIMENSION> centerOfAtlasPhysicalSpace;
  atlasImage->TransformIndexToPhysicalPoint( centerOfAtlas, centerOfAtlasPhysicalSpace );

  WorkingImageType::Pointer rhoImage, phiImage, thetaImage;

  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  CreateNewFloatImageFromTemplate(rhoImage, atlasImage);
  CreateNewFloatImageFromTemplate(phiImage, atlasImage);
  CreateNewFloatImageFromTemplate(thetaImage, atlasImage);

  itk::ImageRegionIterator<WorkingImageType> it(
    atlasImage, atlasImage->GetLargestPossibleRegion() );
  it.GoToBegin();

  itk::ImageRegionIterator<WorkingImageType> rhoit(
    rhoImage, rhoImage->GetLargestPossibleRegion() );
  rhoit.GoToBegin();

  itk::ImageRegionIterator<WorkingImageType> phiit(
    phiImage, phiImage->GetLargestPossibleRegion() );
  phiit.GoToBegin();

  itk::ImageRegionIterator<WorkingImageType> thetait(
    thetaImage, thetaImage->GetLargestPossibleRegion() );
  thetait.GoToBegin();

  itk::Point<float, DIMENSION> currentLocationPhysicalSpace;
  itk::Point<float, DIMENSION> LocationWithRespectToCenterOfImageInMM;

  while( !it.IsAtEnd() )
    {
    const WorkingImageType::IndexType CurrentIndex = it.GetIndex();
    atlasImage->TransformIndexToPhysicalPoint(CurrentIndex, currentLocationPhysicalSpace);
    for( unsigned i = 0; i < DIMENSION; i++ )
      {
      LocationWithRespectToCenterOfImageInMM[i] =
        currentLocationPhysicalSpace[i] - centerOfAtlasPhysicalSpace[i];
      }

    if( CurrentIndex[0] == ( centerOfAtlas[0] ) && CurrentIndex[1] ==
        ( centerOfAtlas[1] ) && CurrentIndex[2] == ( centerOfAtlas[2] ) )
      {
      std::cout << "CENTER_MATH AT (" << CurrentIndex << "): "
                << LocationWithRespectToCenterOfImageInMM << " = "
                << currentLocationPhysicalSpace << " - " << centerOfAtlasPhysicalSpace
                << std::endl;
      }

      {
      float rhoValue, phiValue, thetaValue;
      XYZToSpherical(LocationWithRespectToCenterOfImageInMM, rhoValue, phiValue, thetaValue);
      rhoit.Set(rhoValue);
      phiit.Set(phiValue);
      thetait.Set(thetaValue);
      }

    ++it;
    ++rhoit;
    ++phiit;
    ++thetait;
    }

  std::string RhoMapName = atlasDataSet->GetSpatialLocationFilenameByType("rho");
  std::string PhiMapName = atlasDataSet->GetSpatialLocationFilenameByType("phi");
  std::string ThetaMapName = atlasDataSet->GetSpatialLocationFilenameByType("theta");

  std::string rhoPath = itksys::SystemTools::GetFilenamePath( RhoMapName );
  itksys::SystemTools::MakeDirectory( rhoPath.c_str() );

  std::string phiPath = itksys::SystemTools::GetFilenamePath( PhiMapName );
  itksys::SystemTools::MakeDirectory( phiPath.c_str() );

  std::string thetaPath = itksys::SystemTools::GetFilenamePath( ThetaMapName );
  itksys::SystemTools::MakeDirectory( thetaPath.c_str() );

  // Check if rho,phi and theta file exists.
  itkUtil::WriteImage<WorkingImageType>(rhoImage, RhoMapName);
  itkUtil::WriteImage<WorkingImageType>(phiImage, PhiMapName);
  itkUtil::WriteImage<WorkingImageType>(thetaImage, ThetaMapName);
}

void
BRAINSCutGenerateProbability
::CreateNewFloatImageFromTemplate(WorkingImageType::Pointer & PointerToOutputImage,
                                  const WorkingImageType::Pointer & PreInitializedImage)
{
  WorkingImageType::RegionType region;

  PointerToOutputImage = WorkingImageType::New();
  region.SetSize( PreInitializedImage->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( PreInitializedImage->GetLargestPossibleRegion().GetIndex() );
  PointerToOutputImage->SetLargestPossibleRegion(region);
  PointerToOutputImage->SetBufferedRegion(region);
  PointerToOutputImage->SetRequestedRegion(region);
  PointerToOutputImage->CopyInformation(PreInitializedImage);
  PointerToOutputImage->Allocate();
  PointerToOutputImage->FillBuffer(0.0);
  // NO LONGER NEEDED CHECK_CORONAL(PointerToOutputImage->GetDirection());
  PointerToOutputImage->SetDirection( PreInitializedImage->GetDirection() );
  PointerToOutputImage->SetMetaDataDictionary( PreInitializedImage->GetMetaDataDictionary() );
  itk::ImageRegionIterator<WorkingImageType> bbri( PointerToOutputImage,
                                                   PointerToOutputImage->GetLargestPossibleRegion() );
  bbri = bbri.Begin();
  while( !bbri.IsAtEnd() )
    {
    // Zeroing voxel signal intensity values
    bbri.Set(itk::NumericTraits<WorkingImageType::PixelType>::Zero);
    ++bbri;
    }
}

void
BRAINSCutGenerateProbability
::XYZToSpherical(const itk::Point<float, 3> & LocationWithOriginAtCenterOfImage,
                 float & rhoValue, float & phiValue, float & thetaValue)
{
  /*Rho*/
#define _SQR(a) ( ( a ) * ( a ) )
  rhoValue = static_cast<float>
    ( vcl_sqrt( _SQR(LocationWithOriginAtCenterOfImage[0])
                + _SQR(LocationWithOriginAtCenterOfImage[1])
                + _SQR(LocationWithOriginAtCenterOfImage[2]) ) );
#undef _SQR
  /*Phi*/
  phiValue = 0.0F;
  if( LocationWithOriginAtCenterOfImage[0] < 0 )
    {
    phiValue = vcl_atan2(-LocationWithOriginAtCenterOfImage[0], LocationWithOriginAtCenterOfImage[1]);
    }
  else
    {
    phiValue = vcl_atan2(LocationWithOriginAtCenterOfImage[0], LocationWithOriginAtCenterOfImage[1]);
    }
  /*Theta*/
  thetaValue = 0.0F;
  if( LocationWithOriginAtCenterOfImage[2] < 0 )
    {
    thetaValue = vcl_atan2(-LocationWithOriginAtCenterOfImage[2], LocationWithOriginAtCenterOfImage[1]);
    }
  else
    {
    thetaValue = vcl_atan2(LocationWithOriginAtCenterOfImage[2], LocationWithOriginAtCenterOfImage[1]);
    }

  //  thetaValue = vcl_acos(LocationWithOriginAtCenterOfImage[2]/rhoValue);

  rhoValue = rhoValue / 128.0F;  // The largest brain ever will always fit in a sphere
  // with radius of 128MM centered at the AC point
  phiValue = phiValue / (vnl_math::pi);
  thetaValue = thetaValue / (vnl_math::pi);
}
