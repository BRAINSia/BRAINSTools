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
  this->m_rescaleIntensitiesOutputRange.push_back(40);
  this->m_rescaleIntensitiesOutputRange.push_back(4000);
  this->m_forceACPoint.push_back(0);
  this->m_forcePCPoint.push_back(0);
  this->m_forceVN4Point.push_back(0);
  this->m_forceRPPoint.push_back(0);
  this->m_resultsDir = "./";
  this->m_atlasVolume = "";
  this->m_atlasLandmarks = "";
  this->m_atlasLandmarkWeights = "";

  this->m_outputLandmarksInInputSpaceMap.clear();
  this->m_outputLandmarksInACPCAlignedSpaceMap.clear();
}

bool BRAINSConstellationDetectorPrimary::Compute( void )
{
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(this->m_numberOfThreads);

  // ------------------------------------
  // Read external files
  std::cout << "\nReading in external files..." << std::endl;

  // read in lls model file
  std::map<std::string, std::vector<double> > llsMeans;
  std::map<std::string, MatrixType>           llsMatrices;
  std::map<std::string, double>               searchRadii;

  LLSModel theModel;

  theModel.SetFileName( this->m_llsModel );
  if( theModel.Read() != 0 )
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
  typedef itk::ImageFileReader<DImageType3D> AtlasReaderType;
  AtlasReaderType::Pointer reader = AtlasReaderType::New();
  reader->SetFileName( this->m_inputVolume );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << " Error while reading image file( s ) with ITK:\n "
      << err << std::endl;
    }
  std::cout << "Processing: " << this->m_inputVolume << std::endl;

  DImageType3D::Pointer rescaledInputVolume =
    StandardizeMaskIntensity<DImageType3D, ByteImageType>(reader->GetOutput(),
                                                          ITK_NULLPTR,
                                                          0.0005, 1.0 - 0.0005,
                                                          1, 0.95 * MAX_IMAGE_OUTPUT_VALUE,
                                                          0, MAX_IMAGE_OUTPUT_VALUE);

  typedef itk::CastImageFilter<DImageType3D, SImageType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( rescaledInputVolume );
  caster->Update();

  SImageType::Pointer inputVolume = caster->GetOutput();

    {
    const char * const        metaDataEMSP_FCSVName = "EMSP_FCSV_FILENAME";
    itk::MetaDataDictionary & dict = reader->GetOutput()->GetMetaDataDictionary();
    std::string               ImageMetaDataEMSPFileOverride =  "";
    // if it exists and the string matches what we put in on the image to write, AOK.
    if( itk::ExposeMetaData<std::string>(dict, metaDataEMSP_FCSVName, ImageMetaDataEMSPFileOverride) != false )
      {
      std::string directoryName = itksys::SystemTools::GetParentDirectory(this->m_inputVolume.c_str() );
      if( directoryName == "" )
        {
        directoryName = ".";
        }
      this->m_inputLandmarksEMSP = directoryName + "/" + ImageMetaDataEMSPFileOverride;

      std::cout << "STATUS:  Found meta-data for EMSP override with value: " << this->m_inputLandmarksEMSP << std::endl;
      }
    }

  // load corresponding landmarks in EMSP aligned space from file if possible
  LandmarksMapType landmarksEMSP;
  if( this->m_inputLandmarksEMSP.compare( "" ) != 0 )
    {
    landmarksEMSP = ReadSlicer3toITKLmk( this->m_inputLandmarksEMSP );
    }

  // ------------------------------------
  // Find center of head mass
  std::cout << "\nFinding center of head mass..." << std::endl;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
  findCenterFilter->SetInput( inputVolume );
  findCenterFilter->SetAxis( 2 );
  findCenterFilter->SetOtsuPercentileThreshold( 0.01 );
  findCenterFilter->SetClosingSize( 7 );
  findCenterFilter->SetHeadSizeLimit( 700 );
  findCenterFilter->SetBackgroundValue( 0 );
  findCenterFilter->Update();
  ImagePointType centerOfHeadMass = findCenterFilter->GetCenterOfBrain();

  // ------------------------------------
  // Find eye centers with BRAINS Hough Eye Detector
  HoughEyeDetectorType::Pointer houghEyeDetector = HoughEyeDetectorType::New();

  if( ( landmarksEMSP.find( "LE" ) != landmarksEMSP.end() )
    && ( landmarksEMSP.find( "RE" ) != landmarksEMSP.end() ) )
    {
    std::cout << "\nLoaded eye centers information for BRAINS Hough Eye Detector." << std::endl;
    std::cout << "Skip estimation steps for eye centers." << std::endl;
    }
  else if( this->m_forceHoughEyeDetectorReportFailure == true )
    {
    houghEyeDetector->SetFailure( true );
    std::cout << "\nThe Hough eye detector is doomed to failure as notified." << std::endl;
    std::cout << "Skip estimation steps for eye centers." << std::endl;
    }
  else
    {
    std::cout << "\nFinding eye centers with BRAINS Hough Eye Detector..." << std::endl;
    houghEyeDetector->SetInput( inputVolume );
    houghEyeDetector->SetHoughEyeDetectorMode( this->m_houghEyeDetectorMode );
    houghEyeDetector->SetResultsDir( this->m_resultsDir );           // debug output dir
    houghEyeDetector->SetWritedebuggingImagesLevel( this->m_writedebuggingImagesLevel );
    houghEyeDetector->SetCenterOfHeadMass( centerOfHeadMass );
    try
      {
      houghEyeDetector->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Cannot find eye centers" << std::endl;
      std::cerr << excep << std::endl;
      }
    catch( ... )
      {
      std::cout << "Failed to find eye centers exception occured" << std::endl;
      }
    }

  // ------------------------------------
  // Find MPJ, AC, PC, and VN4 points with BRAINS Constellation Detector
  std::cout << "\nFinding named points with BRAINS Constellation Detector..." << std::endl;

  itk::BRAINSConstellationDetector2<ImageType, ImageType>::Pointer constellation2 = itk::BRAINSConstellationDetector2<ImageType, ImageType>::New();

  if( ( landmarksEMSP.find( "LE" ) != landmarksEMSP.end() )
    && ( landmarksEMSP.find( "RE" ) != landmarksEMSP.end() ) )
    {
    constellation2->SetInput( inputVolume );
    constellation2->SetLandmarksEMSP( landmarksEMSP );
    constellation2->SetCenterOfHeadMass( centerOfHeadMass );
    }
  else
    {
    if( !landmarksEMSP.empty() )
      {
      constellation2->SetLandmarksEMSP( landmarksEMSP );
      }

    if( !houghEyeDetector->GetFailure() )
      {
      ImagePointType houghTransformedCOHM =
        houghEyeDetector->GetInvVersorTransform()->TransformPoint( centerOfHeadMass );

      constellation2->SetLEPoint( houghEyeDetector->GetLE() );
      constellation2->SetREPoint( houghEyeDetector->GetRE() );
      constellation2->SetInput( houghEyeDetector->GetOutput() );
      constellation2->SetHoughEyeTransform( houghEyeDetector->GetModifiableVersorTransform() );
      constellation2->SetCenterOfHeadMass( houghTransformedCOHM );
      }
    else
      {
      constellation2->SetInput( inputVolume );
      constellation2->SetCenterOfHeadMass( centerOfHeadMass );
      }
    }

  // tell the constellation detector whehter Hough eye detector fails
  constellation2->SetHoughEyeFailure( houghEyeDetector->GetFailure() );
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
  constellation2->SetForceACPoint( this->m_forceACPoint );     // In original space
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
  constellation2->SetOriginalInputImage( inputVolume );
  constellation2->SetatlasVolume( this->m_atlasVolume );
  constellation2->SetatlasLandmarks( this->m_atlasLandmarks );
  constellation2->SetatlasLandmarkWeights( this->m_atlasLandmarkWeights );
  constellation2->Update();


  // Save landmarks in input/output or original/aligned space
  this->m_outputLandmarksInACPCAlignedSpaceMap = constellation2->GetAlignedPoints();

  // TODO: Use PrepareOutputsLandmarks here.
  for( LandmarksMapType::const_iterator lit = constellation2->GetAlignedPoints().begin(); lit != constellation2->GetAlignedPoints().end(); ++lit )
    {
    const VersorTransformType::OutputPointType transformedPoint = constellation2->GetOrigToACPCVersorTransform()->TransformPoint( lit->second );
    this->m_outputLandmarksInInputSpaceMap[lit->first] = transformedPoint;
    }

  // ----------------------
  // Write results to disk
  std::cout << "\nWriting results to files..." << std::endl;
  if( this->m_outputTransform.compare( "" ) != 0 )
    {
    TransformWriterType::Pointer writer = TransformWriterType::New();
    writer->SetInput( constellation2->GetOrigToACPCVersorTransform() );
    writer->SetFileName( this->m_outputTransform );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Cannot write the outputTransform file!" << std::endl;
      std::cerr << excep << std::endl;
      }
    std::cout << "The output rigid transform file is written." << std::endl;
    }

  std::string preferedOutputReferenceImage = "";
  if( this->m_outputVolume.compare( "" ) != 0 )
    {
    preferedOutputReferenceImage = this->m_outputVolume;
    // This will be overwritten if outputResampledVolume is set
    // Write the aligned image to a file
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( this->m_outputVolume );
    writer->SetInput( constellation2->GetOutput() );
    writer->SetUseCompression( true );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Cannot write the outputVolume image!" << std::endl;
      std::cerr << excep << std::endl;
      }
    std::cout << "The output unresampled volume is written." << std::endl;
    }

  if( this->m_outputResampledVolume.compare( "" ) != 0 )
    {
    preferedOutputReferenceImage = this->m_outputResampledVolume;
    // This will be overwritten if outputResampledVolume is set
    // Write the aligned image to a file
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( this->m_outputResampledVolume );
    // TODO: Encapsulate this for applying a transform resampling.
    writer->SetInput( constellation2->GetOutputResampledImage() );
    writer->SetUseCompression( true );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Cannot write the outputResampledVolume image!" << std::endl;
      std::cerr << excep << std::endl;
      }
    std::cout << "The output resampled output volume is written." << std::endl;
    }

  // ------------------------
  if( this->m_outputLandmarksInInputSpace.compare( "" ) != 0 )
    {
    WriteITKtoSlicer3Lmk( this->m_outputLandmarksInInputSpace,
      this->m_outputLandmarksInInputSpaceMap );
    std::cout << "The output landmarks list file in the original space is written." << std::endl;
    }

  if( this->m_outputLandmarksInACPCAlignedSpace.compare( "" ) != 0 )
    {
    WriteITKtoSlicer3Lmk( this->m_outputLandmarksInACPCAlignedSpace,
      this->m_outputLandmarksInACPCAlignedSpaceMap );
    std::cout << "The output landmarks list file in the output space is written." << std::endl;
    if( preferedOutputReferenceImage.compare( "" ) == 0 )
      {
      std::cout << "WARNING no aligned output volume is requested." << std::endl;
      }
    }

  // ----------------------
  if( this->m_outputMRML.compare( "" ) != 0 )
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

  if( this->m_outputUntransformedClippedVolume.compare( "" ) != 0 )
    {
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( this->m_outputUntransformedClippedVolume );
    writer->SetInput( constellation2->GetOutputUntransformedClippedVolume() );
    writer->SetUseCompression( true );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Cannot write the outputUntransformedClippedVolume image!" << std::endl;
      std::cerr << excep << std::endl;
      }
    std::cout << "The output untransformed clipped volume is written." << std::endl;
    }

  if( ( this->m_outputVerificationScript.compare( "" ) != 0 )
    && ( this->m_outputLandmarksInACPCAlignedSpace.compare( "" ) != 0 )
    && ( this->m_outputVolume.compare( "" ) != 0 ) )
    {
    writeVerificationScript( this->m_outputVerificationScript, this->m_outputVolume,
      this->m_outputLandmarksInACPCAlignedSpace );
    std::cout << "The verification script is written." << std::endl;
    }
  return EXIT_SUCCESS;
}
