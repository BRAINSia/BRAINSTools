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
#include "BRAINSCutGenerateRegistrations.h"

#include "itkBRAINSROIAutoImageFilter.h"
#include "BRAINSFitHelper.h"

// ----------------------------------------------------- //
BRAINSCutGenerateRegistrations
::BRAINSCutGenerateRegistrations(  BRAINSCutDataHandler& dataHandler ) :
  myDataHandler(nullptr),
  atlasToSubjectRegistraionOn(false),
  subjectDataSets()
{
  myDataHandler =  &dataHandler;
  myDataHandler->SetRegistrationParameters();
  myDataHandler->SetAtlasDataSet();
}

// ----------------------------------------------------- //
void
BRAINSCutGenerateRegistrations
::SetAtlasToSubjectRegistrationOn( bool onOff)
{
  atlasToSubjectRegistraionOn = onOff;
}

// ----------------------------------------------------- //
void
BRAINSCutGenerateRegistrations
::SetDataSet( bool applyDataSet )
{
  // if applyDataSEt==false, then use training dataset
  if( applyDataSet )
    {
    subjectDataSets = myDataHandler->GetApplyDataSet();
    }
  else
    {
    subjectDataSets = myDataHandler->GetTrainDataSet();
    }
}

// ----------------------------------------------------- //
void
BRAINSCutGenerateRegistrations
::GenerateRegistrations()
{
  for( std::list<DataSet *>::iterator subjectIt = subjectDataSets.begin();
       subjectIt != subjectDataSets.end();
       ++subjectIt )
    {
    const std::string subjectFilename( (*subjectIt)->GetImageFilenameByType(
                                         myDataHandler->GetRegistrationImageTypeToUse() ) );

    const RegistrationType *subjectRegistration = (*subjectIt)->GetRegistrationWithID(
        myDataHandler->GetRegistrationID() );
    const std::string SubjectToAtlasRegistrationFilename
      ( subjectRegistration->GetAttribute<StringValue>("SubjToAtlasRegistrationFilename") );
    const std::string AtlasToSubjRegistrationFilename
      ( subjectRegistration->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename") );
    const std::string SubjectBinaryFilename
      ( (*subjectIt)->GetMaskFilenameByType( "RegistrationROI" ) );

    if( atlasToSubjectRegistraionOn &&
        (!itksys::SystemTools::FileExists( AtlasToSubjRegistrationFilename.c_str() ) ) )
      {
      // create directories
      std::string directory = itksys::SystemTools::GetParentDirectory(
          SubjectToAtlasRegistrationFilename.c_str() );

      if( !itksys::SystemTools::FileExists( directory.c_str() ) )
        {
        itksys::SystemTools::MakeDirectory( directory.c_str() );
        }
      CreateTransformFile(  myDataHandler->GetAtlasFilename(),        // moving image
                            subjectFilename,                          // fixed image
                            myDataHandler->GetAtlasBinaryFilename(),  // moving ROI
                            SubjectBinaryFilename,                    // fixed ROI
                            AtlasToSubjRegistrationFilename,
                            false );
      }
    else if( (!atlasToSubjectRegistraionOn) &&
             (!itksys::SystemTools::FileExists( SubjectToAtlasRegistrationFilename.c_str() ) ) )
      {
      // create directories
      std::string directory = itksys::SystemTools::GetParentDirectory(
          SubjectToAtlasRegistrationFilename.c_str() );

      if( !itksys::SystemTools::FileExists( directory.c_str() ) )
        {
        itksys::SystemTools::MakeDirectory( directory.c_str() );
        }
      CreateTransformFile(  subjectFilename,                         // moving image
                            myDataHandler->GetAtlasFilename(),       // fixed image
                            SubjectBinaryFilename,                   // moving ROI
                            myDataHandler->GetAtlasBinaryFilename(), // fixed ROI
                            SubjectToAtlasRegistrationFilename,
                            false );
      }
    }
}

void
BRAINSCutGenerateRegistrations
::CreateTransformFile(const std::string & MovingImageFilename,
                      const std::string & FixedImageFilename,
                      const std::string & MovingBinaryImageFilename,
                      const std::string & FixedBinaryImageFilename,
                      const std::string & OutputRegName,
                      bool verbose)
{
  // Create Helper Class of BRAINSFit for BSpline Registraion.
  // CAUTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Check if input images are all 'float'
  // BRAINSFit takes only float type images!
  // For now, BRAINSCut only uses float type, so it is safe!!!!!!

  typedef itk::BRAINSFitHelper BSplineRegistrationHelperType;
  BSplineRegistrationHelperType::Pointer BSplineRegistrationHelper =
    BSplineRegistrationHelperType::New();

  // Set BSpline Trainsformation Type
  std::vector<std::string> transformType( 4 );

  // transformType[0]="Rigid";
  transformType[0] = "ScaleVersor3D";
  transformType[1] = "ScaleSkewVersor3D";
  transformType[2] = "Affine";
  transformType[3] = "BSpline";

  BSplineRegistrationHelper->SetTransformType( transformType );

  // Set Displacement Field Size
  BSplineRegistrationHelper->SetMaxBSplineDisplacement(4);

  // BSpline Grid Size
  std::vector<int> splineGridSize(3);

  splineGridSize[0] = 28;
  splineGridSize[1] = 20;
  splineGridSize[2] = 24;

  BSplineRegistrationHelper->SetSplineGridSize( splineGridSize );

  // Number of Iterations
  std::vector<int> numberOfIterations(4);

  numberOfIterations[0] = 1500;
  numberOfIterations[1] = 1500;
  numberOfIterations[2] = 1500;
  numberOfIterations[3] = 1500;
  // numberOfIterations[4]=1500;

  BSplineRegistrationHelper->SetNumberOfIterations( numberOfIterations );

  // Set Minimum Step Length
  std::vector<double> minimumStepLength(4);

  minimumStepLength[0] = 0.005;
  minimumStepLength[1] = 0.005;
  minimumStepLength[2] = 0.005;
  minimumStepLength[3] = 0.005;
  // minimumStepLength[4]=0.005;

  BSplineRegistrationHelper->SetMinimumStepLength( minimumStepLength );

  // Set Histogram Matching Options
  // BSplineRegistrationHelper->SetNumberOfHistogramBins(250);
  // TODO: decide what to do here too much freedom?
  // BSplineRegistrationHelper->SetNumberOfMatchPoints(10);
  // TODO:: HISTOGRAMMATCHING OPTION is now Disabled!
  BSplineRegistrationHelper->SetHistogramMatch(false);

  // Set Fixed Volume
  WorkingImageType::Pointer fixedVolume = ReadImageByFilename( FixedImageFilename );

  BSplineRegistrationHelper->SetFixedVolume( fixedVolume );

  // Set Moving Volume
  WorkingImageType::Pointer movingVolume = ReadImageByFilename( MovingImageFilename );

  BSplineRegistrationHelper->SetMovingVolume( movingVolume );

  typedef itk::Image<unsigned char, 3> LocalBinaryImageType;

  // - Fixed Image Binary Mask

  if( FixedBinaryImageFilename == "NA" ||
      FixedBinaryImageFilename == "na" ||
      FixedBinaryImageFilename == "" )
    {
    typedef itk::BRAINSROIAutoImageFilter<WorkingImageType,
                                          LocalBinaryImageType> ROIAutoFilterType;

    ROIAutoFilterType::Pointer fixedVolumeROIFilter = ROIAutoFilterType::New();

    fixedVolumeROIFilter->SetInput( fixedVolume );
    fixedVolumeROIFilter->SetDilateSize(myDataHandler->GetROIAutoDilateSize() );
    fixedVolumeROIFilter->Update();

    BSplineRegistrationHelper->SetFixedBinaryVolume(
      fixedVolumeROIFilter->GetSpatialObjectROI() );
    // get brain region size
    }
  else
    {
    typedef itk::ImageFileReader<LocalBinaryImageType> BinaryImageReaderType;
    BinaryImageReaderType::Pointer binaryFixedImageReader = BinaryImageReaderType::New();
    binaryFixedImageReader->SetFileName( FixedBinaryImageFilename );
    binaryFixedImageReader->Update();

    typedef itk::ImageMaskSpatialObject<3> binarySpatialObjectType;
    binarySpatialObjectType::Pointer binaryFixedObject
      = binarySpatialObjectType::New();
    binaryFixedObject->SetImage( binaryFixedImageReader->GetOutput() );

    BSplineRegistrationHelper->SetFixedBinaryVolume( binaryFixedObject );
    }

  // - Moving Image Binary Mask

  if( MovingBinaryImageFilename == "NA" ||
      MovingBinaryImageFilename == "na" ||
      MovingBinaryImageFilename == "" )
    {
    typedef itk::BRAINSROIAutoImageFilter<WorkingImageType,
                                          LocalBinaryImageType> ROIAutoFilterType;

    ROIAutoFilterType::Pointer movingVolumeROIFilter = ROIAutoFilterType::New();

    movingVolumeROIFilter->SetInput( movingVolume );
    movingVolumeROIFilter->SetDilateSize(myDataHandler->GetROIAutoDilateSize() );
    movingVolumeROIFilter->Update();

    BSplineRegistrationHelper->SetMovingBinaryVolume(
      movingVolumeROIFilter->GetSpatialObjectROI() );
    }
  else
    {
    typedef itk::ImageFileReader<LocalBinaryImageType> BinaryImageReaderType;
    BinaryImageReaderType::Pointer binaryMovingImageReader =
      BinaryImageReaderType::New();
    binaryMovingImageReader->SetFileName( MovingBinaryImageFilename );
    binaryMovingImageReader->Update();

    typedef itk::ImageMaskSpatialObject<3> binarySpatialObjectType;
    binarySpatialObjectType::Pointer binaryMovingObject
      = binarySpatialObjectType::New();
    binaryMovingObject->SetImage( binaryMovingImageReader->GetOutput() );

    BSplineRegistrationHelper->SetMovingBinaryVolume( binaryMovingObject );
    }

  // Set Other Options

  // TODO :: Following has to be user input
  constexpr unsigned int maxBSplineDisplacement = 7;
  const double       maskInferiorCutOffFromCenter = 65.0;
  const double       translationScale = 1000.0;
  const double       reproportionalScale = 1.0;
  const double       skewScale = 1.0;

  BSplineRegistrationHelper->SetSamplingPercentage( 0.05 ); // Use 5% of image
  BSplineRegistrationHelper->SetMaxBSplineDisplacement( maxBSplineDisplacement );
  BSplineRegistrationHelper->SetInitializeTransformMode( "useCenterOfHeadAlign" );
  BSplineRegistrationHelper->SetMaskInferiorCutOffFromCenter(
    maskInferiorCutOffFromCenter );
  // BSplineRegistrationHelper->SetOutputTransform( OutputRegName );
  BSplineRegistrationHelper->SetTranslationScale( translationScale );
  BSplineRegistrationHelper->SetReproportionScale( reproportionalScale );
  BSplineRegistrationHelper->SetSkewScale( skewScale );

  //  Start Registration

  // TODO: is this line really print before start registration????

  if( verbose == true )
    {
    BSplineRegistrationHelper->PrintCommandLine(true, "BSplineRegistrationHelper");
    }

  BSplineRegistrationHelper->Update();

  if( verbose == true )
    {
    std::cout << " - Write deformation " << std::endl
              << " :: " << OutputRegName
              << std::endl;
    }
  itk::WriteTransformToDisk<double>( BSplineRegistrationHelper->GetCurrentGenericTransform(),
                        OutputRegName );
  // Write out Transformed Output As Well
  // - EX. from GenericTransformImage.hxx

  WorkingImageType::Pointer DeformedMovingImage;
  DeformedMovingImage = TransformResample<WorkingImageType,
                                          WorkingImageType>(
      movingVolume.GetPointer(),
      fixedVolume.GetPointer(),
      0.0F,
      GetInterpolatorFromString<WorkingImageType>("Linear").GetPointer(),
      BSplineRegistrationHelper->GetCurrentGenericTransform().GetPointer() );

  typedef itk::ImageFileWriter<WorkingImageType> DeformedVolumeWriterType;

  DeformedVolumeWriterType::Pointer deformedVolumeWriter = DeformedVolumeWriterType::New();

  deformedVolumeWriter->SetFileName( OutputRegName + "_output.nii.gz" );
  deformedVolumeWriter->SetInput( DeformedMovingImage );
  deformedVolumeWriter->Update();
}
