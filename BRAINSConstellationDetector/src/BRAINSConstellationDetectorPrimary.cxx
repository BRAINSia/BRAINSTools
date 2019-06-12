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
/*
 */
#include "BRAINSConstellationDetectorPrimary.h"

BRAINSConstellationDetectorPrimary::BRAINSConstellationDetectorPrimary()
{
  this->m_houghEyeDetectorMode = 1;
  this->m_mspQualityLevel = 2;
  this->m_writedebuggingImagesLevel = 0;
  this->m_numberOfThreads = -1;
  this->m_otsuPercentileThreshold = 0.01;
  this->m_acLowerBound = 1000.0;
  this->m_trimRescaledIntensities = 4.4172;

  this->m_radiusMPJ = -1;
  this->m_radiusAC = -1;
  this->m_radiusPC = -1;
  this->m_radiusVN4 = -1;

  this->m_cutOutHeadInOutputVolume = false;
  this->m_rescaleIntensities = false;
  this->m_forceHoughEyeDetectorReportFailure = false;
  this->m_debug = false;
  this->m_verbose = false;

  this->m_inputTemplateModel = itksys::SystemTools::GetProgramPath( pathOut.c_str() ) + "/" + "T1.mdl";
  this->m_llsModel = itksys::SystemTools::GetProgramPath( pathOut.c_str() ) + "/" + "LLSModel.h5";

  this->m_inputVolume = "";
  this->m_outputVolume = "";
  this->m_outputResampledVolume = "";
  this->m_outputTransform = "";
  this->m_outputLandmarksInInputSpace = "";
  this->m_outputLandmarksInACPCAlignedSpace = "";
  this->m_outputMRML = "";
  this->m_outputVerificationScript = "";
  this->m_outputUntransformedClippedVolume = "";
  this->m_inputLandmarksEMSP = "";
  this->m_writeBranded2DImage = "";
  this->m_backgroundFillValueString = "0";
  this->m_interpolationMode = "Linear";
  this->m_rescaleIntensitiesOutputRange.push_back( 40 );
  this->m_rescaleIntensitiesOutputRange.push_back( 4000 );
  this->m_forceACPoint.push_back( 0 );
  this->m_forcePCPoint.push_back( 0 );
  this->m_forceVN4Point.push_back( 0 );
  this->m_forceRPPoint.push_back( 0 );
  this->m_resultsDir = "./";
  this->m_atlasVolume = "";
  this->m_atlasLandmarks = "";
  this->m_atlasLandmarkWeights = "";

  this->m_outputLandmarksInInputSpaceMap.clear();
  this->m_outputLandmarksInACPCAlignedSpaceMap.clear();
}

BRAINSConstellationDetectorPrimary::ImagePointType
BRAINSConstellationDetectorPrimary ::localFindCenterHeadFunc(
  BRAINSConstellationDetectorPrimary::ImageType::ConstPointer img )
{
  // ------------------------------------
  // Find center of head mass
  std::cout << "\nFinding center of head mass..." << std::endl;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();

  findCenterFilter->SetAxis( 2 );
  findCenterFilter->SetOtsuPercentileThreshold( 0.01 );
  findCenterFilter->SetClosingSize( 7 );
  findCenterFilter->SetHeadSizeLimit( 700 );
  findCenterFilter->SetBackgroundValue( 0 );
  findCenterFilter->SetInput( img );
  findCenterFilter->Update();
  const ImagePointType centerOfHeadMass = findCenterFilter->GetCenterOfBrain();
  return centerOfHeadMass;
}

bool
BRAINSConstellationDetectorPrimary::Compute( void )
{
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder( this->m_numberOfThreads );

  // ------------------------------------
  // Read external files
  std::cout << "\nReading in external files..." << std::endl;

  // read in lls model file
  std::map< std::string, std::vector< double > >  llsMeans;
  std::map< std::string, LandmarkIO::MatrixType > llsMatrices;
  std::map< std::string, double >                 searchRadii;

  LLSModel theModel;

  theModel.SetFileName( this->m_llsModel );
  if ( theModel.Read() != 0 )
  {
    std::cerr << "Error reading LLS Model" << std::endl;
    return EXIT_FAILURE;
  }
  llsMeans = theModel.GetLLSMeans();
  llsMatrices = theModel.GetLLSMatrices();
  searchRadii = theModel.GetSearchRadii();

  // ------------------------------------
  // load image
  std::cout << "\nLoading image..." << std::endl;

  // Input image is read as a double image;
  // then it is rescaled to a specific dynamic range;
  // Finally it is cast to a Short type image.
  using AtlasReaderType = itk::ImageFileReader< DImageType3D >;
  AtlasReaderType::Pointer reader = AtlasReaderType::New();
  reader->SetFileName( this->m_inputVolume );
  try
  {
    reader->Update();
  }
  catch ( itk::ExceptionObject & err )
  {
    std::cerr << " Error while reading image file( s ) with ITK:\n " << err << std::endl;
  }
  std::cout << "Processing: " << this->m_inputVolume << std::endl;

  DImageType3D::Pointer rescaledInputVolume = StandardizeMaskIntensity< DImageType3D, ByteImageType >(
    reader->GetOutput(), nullptr, 0.0005, 1.0 - 0.0005, 1, 0.95 * MAX_IMAGE_OUTPUT_VALUE, 0, MAX_IMAGE_OUTPUT_VALUE );

  using CasterType = itk::CastImageFilter< DImageType3D, SImageType >;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( rescaledInputVolume );
  caster->Update();

  SImageType::Pointer origSpaceInputVolume = caster->GetOutput();

  /*
   * Look for existing manually identified landmark point files adjacent
   * to the original image.
   */
  if ( this->m_inputLandmarksEMSP.empty() ) // Only look for default files if not specified on command line.
  {
    const static std::string fcsv_extension{ ".fcsv" };
    const std::string        root_dir = itksys::SystemTools::GetParentDirectory( this->m_inputVolume );
    std::string              potentialLandmarkFileName =
      root_dir + "/" + itksys::SystemTools::GetFilenameWithoutExtension( this->m_inputVolume );

    // If the extension is .gz, then also remove the preceding extension.
    if ( itksys::SystemTools::GetFilenameExtension( this->m_inputVolume ) == ".gz" )
    {
      potentialLandmarkFileName =
        root_dir + "/" + itksys::SystemTools::GetFilenameWithoutExtension( potentialLandmarkFileName );
    }
    potentialLandmarkFileName += fcsv_extension;
    if ( itksys::SystemTools::FileExists( potentialLandmarkFileName, true ) )
    {
      std::cerr << "WARNING: Using the side-car landmark override file to pre-load landmarks!: "
                << "\n         " << potentialLandmarkFileName << std::endl;
      this->m_inputLandmarksEMSP = potentialLandmarkFileName;
    }
    else
    {
      std::cerr << "WARNING: Side-car landmark file not found to pre-load landmarks!"
                << "\n         " << potentialLandmarkFileName << std::endl;
      // Now looking for file encoded as meta data in header.
      const char * const        metaDataEMSP_FCSVName = "EMSP_FCSV_FILENAME";
      itk::MetaDataDictionary & dict = reader->GetOutput()->GetMetaDataDictionary();
      std::string               ImageMetaDataEMSPFileOverride = "";
      // if it exists and the string matches what we put in on the image to write, AOK.
      if ( itk::ExposeMetaData< std::string >( dict, metaDataEMSP_FCSVName, ImageMetaDataEMSPFileOverride ) != false )
      {
        std::string directoryName = itksys::SystemTools::GetParentDirectory( this->m_inputVolume.c_str() );
        if ( directoryName == "" )
        {
          directoryName = ".";
        }
        this->m_inputLandmarksEMSP = directoryName + "/" + ImageMetaDataEMSPFileOverride;

        std::cout << "STATUS:  Found meta-data for EMSP override with value: " << this->m_inputLandmarksEMSP
                  << std::endl;
      }
    }
    std::cerr << "\n\n\n\n\n";
  }

  // load corresponding landmarks in EMSP aligned space from file if possible
  LandmarksMapType origSpaceLandmarks;
  if ( !this->m_inputLandmarksEMSP.empty() )
  {
    origSpaceLandmarks = ReadSlicer3toITKLmk( this->m_inputLandmarksEMSP );
  }


  // ------------------------------------
  // Find eye centers with BRAINS Hough Eye Detector
  HoughEyeDetectorType::Pointer                  houghEyeDetector = HoughEyeDetectorType::New();
  itk::VersorRigid3DTransform< double >::Pointer orig2eyeFixedImageVersorTransform = nullptr;
  itk::VersorRigid3DTransform< double >::Pointer org2eyeFixedLandmarkVersorTransform = nullptr;


  LandmarksMapType   eyeFixedSpaceLandmarks;
  ImageType::Pointer eyeFixedResampledImage = nullptr;

  bool hasBothEyeLandmarksDefinedInOrigSpace = false;
  if ( ( origSpaceLandmarks.find( "LE" ) != origSpaceLandmarks.end() ) &&
       ( origSpaceLandmarks.find( "RE" ) != origSpaceLandmarks.end() ) )
  {
    hasBothEyeLandmarksDefinedInOrigSpace = true;
  }

  if ( hasBothEyeLandmarksDefinedInOrigSpace )
  {
    std::cout << "\nLoaded eye centers information for BRAINS Hough Eye Detector." << std::endl;
    std::cout << "Skip estimation steps for eye centers." << std::endl;
    const SImageType::PointType origSpaceLE = origSpaceLandmarks.find( "LE" )->second;
    const SImageType::PointType origSpaceRE = origSpaceLandmarks.find( "RE" )->second;

    orig2eyeFixedImageVersorTransform =
      itk::ResampleFromEyePoints< ImageType, ImageType >( origSpaceLE, origSpaceRE, origSpaceInputVolume );
    eyeFixedResampledImage = itk::RigidResampleInPlayByVersor3D< ImageType, ImageType >(
      origSpaceInputVolume, orig2eyeFixedImageVersorTransform );

    // Now transform landmarks based on the eye centerings.
    org2eyeFixedLandmarkVersorTransform = itk::VersorRigid3DTransform< double >::New();
    orig2eyeFixedImageVersorTransform->GetInverse( org2eyeFixedLandmarkVersorTransform );
    for ( auto & landmark_pair : origSpaceLandmarks )
    {
      eyeFixedSpaceLandmarks[landmark_pair.first] =
        org2eyeFixedLandmarkVersorTransform->TransformPoint( landmark_pair.second );
    }

    // If CM is not provided from original space and moved to eyeFixedSpace, then compute it here from the eyeFixedImage
    if ( eyeFixedSpaceLandmarks.find( "CM" ) == eyeFixedSpaceLandmarks.end() )
    {
      eyeFixedSpaceLandmarks["CM"] = localFindCenterHeadFunc( eyeFixedResampledImage );
    }
  }
  else
  {
    std::cout << "\nFinding eye centers with BRAINS Hough Eye Detector..." << std::endl;
    houghEyeDetector->SetInput( origSpaceInputVolume );
    houghEyeDetector->SetHoughEyeDetectorMode( this->m_houghEyeDetectorMode );
    houghEyeDetector->SetResultsDir( this->m_resultsDir ); // debug output dir
    houghEyeDetector->SetWritedebuggingImagesLevel( this->m_writedebuggingImagesLevel );

    origSpaceLandmarks["CM"] = localFindCenterHeadFunc( origSpaceInputVolume );
    houghEyeDetector->SetCenterOfHeadMass( origSpaceLandmarks.at( "CM" ) );

    try
    {
      houghEyeDetector->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      std::cerr << "Cannot find eye centers" << std::endl;
      std::cerr << excep << std::endl;
    }
    catch ( ... )
    {
      std::cout << "Failed to find eye centers exception occurred" << std::endl;
    }
    org2eyeFixedLandmarkVersorTransform = houghEyeDetector->GetModifiableVersorTransform();
    std::cout << "PRE COHM: " << origSpaceLandmarks.at( "CM" );

    eyeFixedSpaceLandmarks["CM"] =
      houghEyeDetector->GetInvVersorTransform()->TransformPoint( origSpaceLandmarks.at( "CM" ) );
    std::cout << "POST HOUGH INV COHM: " << eyeFixedSpaceLandmarks.at( "CM" );

    origSpaceLandmarks["LE"] = houghEyeDetector->GetLE();
    origSpaceLandmarks["RE"] = houghEyeDetector->GetRE();
    eyeFixedResampledImage = houghEyeDetector->GetOutput();
  }

  // ------------------------------------
  // Find MPJ, AC, PC, and VN4 points with BRAINS Constellation Detector
  std::cout << "\nFinding named points with BRAINS Constellation Detector..." << std::endl;
  itk::BRAINSConstellationDetector2< ImageType, ImageType >::Pointer constellation2 =
    itk::BRAINSConstellationDetector2< ImageType, ImageType >::New();

  //  if( !origSpaceLandmarks.empty() )
  //    {
  //      constellation2->SetLandmarksEMSP( origSpaceLandmarks );
  //    }

  constellation2->SetLEPoint( origSpaceLandmarks.at( "LE" ) );
  constellation2->SetREPoint( origSpaceLandmarks.at( "RE" ) );
  constellation2->SetCenterOfHeadMass( eyeFixedSpaceLandmarks.at( "CM" ) ); // This is likely wrong!
  constellation2->SetHoughEyeTransform( org2eyeFixedLandmarkVersorTransform );
  constellation2->SetInput( eyeFixedResampledImage );


  // HACK: --- REMOVE ME
  std::string dggpath = "/tmp/nolmkinit_";
  if ( hasBothEyeLandmarksDefinedInOrigSpace )
  {
    dggpath = "/tmp/manuallmkinit_";
  }
  std::cout << "\n\nHACK: REMOVE ME\n\n" << std::endl;
  std::cout << "\n\nHACK: REMOVE ME\n\n" << std::endl;
  std::cout << "\n\nHACK: REMOVE ME\n\n" << std::endl;
  std::cout << "\n\nHACK: REMOVE ME\n\n" << std::endl;
  std::cout << "\n\nHACK: REMOVE ME\n\n" << std::endl;
  std::cout << "\n\nHACK: REMOVE ME\n\n" << std::endl;
  std::cout << "\n\nHACK: REMOVE ME\n\n" << std::endl;

  WriteITKtoSlicer3Lmk( dggpath + "eye_orig_space.fcsv", origSpaceLandmarks );
  WriteITKtoSlicer3Lmk( dggpath + "eye_fixed_space.fcsv", eyeFixedSpaceLandmarks );

  WriterType::Pointer dggwriter = WriterType::New();
  dggwriter->SetFileName( dggpath + "eye_fixed.nii.gz" );
  dggwriter->SetInput( eyeFixedResampledImage );
  dggwriter->Update();
  dggwriter->SetFileName( dggpath + "eye_orig_space.nii.gz" );
  dggwriter->SetInput( origSpaceInputVolume );
  dggwriter->Update();


  // tell the constellation detector if Hough eye detector fails
  // HACK constellation2->SetHoughEyeFailure( houghEyeDetector->GetFailure() );
  constellation2->SetInputTemplateModel( this->m_inputTemplateModel );
  constellation2->SetMspQualityLevel( this->m_mspQualityLevel );
  constellation2->SetOtsuPercentileThreshold( this->m_otsuPercentileThreshold );
  constellation2->SetAcLowerBound( this->m_acLowerBound );
  constellation2->SetCutOutHeadInOutputVolume( this->m_cutOutHeadInOutputVolume );
  constellation2->SetRescaleIntensities( this->m_rescaleIntensities );
  constellation2->SetTrimRescaledIntensities( this->m_trimRescaledIntensities );
  constellation2->SetRescaleIntensitiesOutputRange( this->m_rescaleIntensitiesOutputRange );
  constellation2->SetBackgroundFillValueString( this->m_backgroundFillValueString );
  constellation2->SetInterpolationMode( this->m_interpolationMode );
  constellation2->SetForceACPoint( this->m_forceACPoint ); // In original space
  constellation2->SetForcePCPoint( this->m_forcePCPoint );
  constellation2->SetForceVN4Point( this->m_forceVN4Point );
  constellation2->SetForceRPPoint( this->m_forceRPPoint );
  constellation2->SetRadiusMPJ( this->m_radiusMPJ );
  constellation2->SetRadiusAC( this->m_radiusAC );
  constellation2->SetRadiusPC( this->m_radiusPC );
  constellation2->SetRadiusVN4( this->m_radiusVN4 );
  constellation2->SetDebug( this->m_debug );
  constellation2->SetVerbose( this->m_verbose );
  constellation2->SetWritedebuggingImagesLevel( this->m_writedebuggingImagesLevel );
  constellation2->SetWriteBranded2DImage( this->m_writeBranded2DImage );
  constellation2->SetResultsDir( this->m_resultsDir );
  constellation2->SetLlsMatrices( llsMatrices );
  constellation2->SetLlsMeans( llsMeans );
  constellation2->SetSearchRadii( searchRadii );
  constellation2->SetOriginalInputImage( origSpaceInputVolume );
  constellation2->SetatlasVolume( this->m_atlasVolume );
  constellation2->SetatlasLandmarks( this->m_atlasLandmarks );
  constellation2->SetatlasLandmarkWeights( this->m_atlasLandmarkWeights );
  constellation2->Update();


  // Save landmarks in input/output or original/aligned space
  this->m_outputLandmarksInACPCAlignedSpaceMap = constellation2->GetAlignedPoints();

  // TODO: Use PrepareOutputsLandmarks here.
  for ( LandmarksMapType::const_iterator lit = constellation2->GetAlignedPoints().begin();
        lit != constellation2->GetAlignedPoints().end();
        ++lit )
  {
    const VersorTransformType::OutputPointType transformedPoint =
      constellation2->GetOrigToACPCVersorTransform()->TransformPoint( lit->second );
    this->m_outputLandmarksInInputSpaceMap[lit->first] = transformedPoint;
  }

  // ----------------------
  // Write results to disk
  std::cout << "\nWriting results to files..." << std::endl;
  if ( this->m_outputTransform.compare( "" ) != 0 )
  {
    TransformWriterType::Pointer writer = TransformWriterType::New();
    writer->SetInput( constellation2->GetOrigToACPCVersorTransform() );
    writer->SetFileName( this->m_outputTransform );
    try
    {
      writer->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      std::cerr << "Cannot write the outputTransform file!" << std::endl;
      std::cerr << excep << std::endl;
    }
    std::cout << "The output rigid transform file is written." << std::endl;
  }

  std::string preferedOutputReferenceImage = "";
  if ( this->m_outputVolume.compare( "" ) != 0 )
  {
    preferedOutputReferenceImage = this->m_outputVolume;
    // This will be overwritten if outputResampledVolume is set
    // Write the aligned image to a file
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( this->m_outputVolume );
    writer->SetInput( constellation2->GetOutput() );
#if ITK_VERSION_MAJOR >= 5
    writer->SetUseCompression( true );
#endif
    try
    {
      writer->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      std::cerr << "Cannot write the outputVolume image!" << std::endl;
      std::cerr << excep << std::endl;
    }
    std::cout << "The output unresampled volume is written." << std::endl;
  }

  if ( this->m_outputResampledVolume.compare( "" ) != 0 )
  {
    preferedOutputReferenceImage = this->m_outputResampledVolume;
    // This will be overwritten if outputResampledVolume is set
    // Write the aligned image to a file
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( this->m_outputResampledVolume );
    // TODO: Encapsulate this for applying a transform resampling.
    writer->SetInput( constellation2->GetOutputResampledImage() );
#if ITK_VERSION_MAJOR >= 5
    writer->SetUseCompression( true );
#endif
    try
    {
      writer->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      std::cerr << "Cannot write the outputResampledVolume image!" << std::endl;
      std::cerr << excep << std::endl;
    }
    std::cout << "The output resampled output volume is written." << std::endl;
  }

  // ------------------------
  if ( this->m_outputLandmarksInInputSpace.compare( "" ) != 0 )
  {
    WriteITKtoSlicer3Lmk( this->m_outputLandmarksInInputSpace, this->m_outputLandmarksInInputSpaceMap );
    std::cout << "The output landmarks list file in the original space is written." << std::endl;
  }

  if ( this->m_outputLandmarksInACPCAlignedSpace.compare( "" ) != 0 )
  {
    WriteITKtoSlicer3Lmk( this->m_outputLandmarksInACPCAlignedSpace, this->m_outputLandmarksInACPCAlignedSpaceMap );
    std::cout << "The output landmarks list file in the output space is written." << std::endl;
    if ( preferedOutputReferenceImage.compare( "" ) == 0 )
    {
      std::cout << "WARNING no aligned output volume is requested." << std::endl;
    }
  }

  // ----------------------
  if ( this->m_outputMRML.compare( "" ) != 0 )
  {
    WriteMRMLFile( this->m_outputMRML,
                   this->m_outputLandmarksInInputSpace,
                   this->m_outputLandmarksInACPCAlignedSpace,
                   this->m_inputVolume,
                   preferedOutputReferenceImage,
                   this->m_outputTransform,
                   this->m_outputLandmarksInInputSpaceMap,
                   this->m_outputLandmarksInACPCAlignedSpaceMap,
                   constellation2->GetOrigToACPCVersorTransform() );
    std::cout << "The output mrml scene file is written." << std::endl;
  }

  if ( this->m_outputUntransformedClippedVolume.compare( "" ) != 0 )
  {
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( this->m_outputUntransformedClippedVolume );
    writer->SetInput( constellation2->GetOutputUntransformedClippedVolume() );
#if ITK_VERSION_MAJOR >= 5
    writer->SetUseCompression( true );
#endif
    try
    {
      writer->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      std::cerr << "Cannot write the outputUntransformedClippedVolume image!" << std::endl;
      std::cerr << excep << std::endl;
    }
    std::cout << "The output untransformed clipped volume is written." << std::endl;
  }

  if ( ( this->m_outputVerificationScript.compare( "" ) != 0 ) &&
       ( this->m_outputLandmarksInACPCAlignedSpace.compare( "" ) != 0 ) && ( this->m_outputVolume.compare( "" ) != 0 ) )
  {
    writeVerificationScript(
      this->m_outputVerificationScript, this->m_outputVolume, this->m_outputLandmarksInACPCAlignedSpace );
    std::cout << "The verification script is written." << std::endl;
  }
  return EXIT_SUCCESS;
}
