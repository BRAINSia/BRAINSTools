#include "BRAINSCutGenerateProbability.h"
#include "NetConfigurationParser.h"
#include "NetConfiguration.h"
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include "Utilities.h"
#include "itkIO.h"

/** constructors */
BRAINSCutGenerateProbability
::BRAINSCutGenerateProbability( std::string netConfigurationFilename ) : BRAINSCutPrimary( netConfigurationFilename )
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  SetRegistrationParametersFromNetConfiguration();

  SetAtlasDataSet();
  SetAtlasImage();

  SetRegionsOfInterestFromNetConfiguration();
  SetTrainingDataSetsList();
}

/** Get/Set Methods */

void
BRAINSCutGenerateProbability
::SetTrainingDataSetsList()
{
  try
    {
    trainingDataSetList = BRAINSCutNetConfiguration.GetTrainDataSets();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_SUCCESS);
    }
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
  GenerateRegistrations(BRAINSCutNetConfiguration, false, false,  1);

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

      std::string currentRegistrationFilename =
        (*currentSubjectIt)->GetRegistrationWithID( registrationID )
        ->GetAttribute<StringValue>( "SubjToAtlasRegistrationFIlename");

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

    float rho, phi, theta;

    XYZToSpherical(LocationWithRespectToCenterOfImageInMM, rho, phi, theta);

    rhoit.Set(rho);
    phiit.Set(phi);
    thetait.Set(theta);

    ++it;    ++rhoit;    ++phiit;    ++thetait;
    }

  std::string RhoMapName = atlasDataSet->GetSpatialLocationFilenameByType("rho");
  std::string PhiMapName = atlasDataSet->GetSpatialLocationFilenameByType("phi");
  std::string ThetaMapName = atlasDataSet->GetSpatialLocationFilenameByType("theta");

  // Check if rho,phi and theta file exists.
  itkUtil::WriteImage<WorkingImageType>(rhoImage, RhoMapName);
  itkUtil::WriteImage<WorkingImageType>(phiImage, PhiMapName);
  itkUtil::WriteImage<WorkingImageType>(thetaImage, ThetaMapName);
}

WorkingImagePointer
BRAINSCutGenerateProbability
::SmoothImage( const WorkingImagePointer image, const float GaussianValue)
{
  if( GaussianValue < 0 + FLOAT_TOLERANCE )
    {
    std::cout << "Gaussian value is less than tolerance. "
              << "No smoothing occurs at this time"
              << std::endl;
    return image;
    }
  typedef itk::SmoothingRecursiveGaussianImageFilter<WorkingImageType, WorkingImageType> SmoothingFilterType;
  SmoothingFilterType::Pointer smoothingFilter = SmoothingFilterType::New();

  smoothingFilter->SetInput( image);
  smoothingFilter->SetSigma( GaussianValue );

  smoothingFilter->Update();

  return smoothingFilter->GetOutput();
}
