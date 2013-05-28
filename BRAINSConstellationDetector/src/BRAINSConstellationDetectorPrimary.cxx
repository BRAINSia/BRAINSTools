/*
 */
#include "BRAINSConstellationDetectorPrimary.h"

BRAINSConstellationDetectorPrimary::BRAINSConstellationDetectorPrimary()
{
  this->m_houghEyeDetectorMode = 1;
  this->m_mspQualityLevel = 2;
  this->m_writedebuggingImagesLevel = 0;
  this->m_numberOfThreads = 1;
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
  this->m_outputLandmarkWeights = "";   // SHOULD BE DELETED
  this->m_inputLandmarksPaired = "";    // SHOULD BE DELETED
  this->m_outputLandmarksPaired = "";   // SHOULD BE DELETED
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

  ReaderType::Pointer reader = ReaderType::New();
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

  // std::cout << "original input file : " << reader->GetOutput() << std::endl;  //Added by Ali

  // ------------------------------------
  // Find center of head mass
  std::cout << "\nFinding center of head mass..." << std::endl;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
  findCenterFilter->SetInput( reader->GetOutput() );
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
    houghEyeDetector->SetInput( reader->GetOutput() );
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

  Constellation2Type::Pointer constellation2 = Constellation2Type::New();

  if( ( landmarksEMSP.find( "LE" ) != landmarksEMSP.end() )
      && ( landmarksEMSP.find( "RE" ) != landmarksEMSP.end() ) )
    {
    constellation2->SetInput( reader->GetOutput() );
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
      constellation2->SetInput( reader->GetOutput() );
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
  constellation2->SetOriginalInputImage( reader->GetOutput() );
  constellation2->Update();

  VersorTransformType::Pointer acpcTransformFromConstellation2 = VersorTransformType::New();
  acpcTransformFromConstellation2->SetFixedParameters( constellation2->GetVersorTransform()->GetFixedParameters() );
  acpcTransformFromConstellation2->SetParameters( constellation2->GetVersorTransform()->GetParameters() );

  // Get the final transform
  VersorTransformType::Pointer finalTransform = VersorTransformType::New();
  VersorTransformType::Pointer invFinalTransform = VersorTransformType::New();
  if( this->m_atlasVolume.empty() )
    {
    finalTransform->SetFixedParameters( constellation2->GetVersorTransform()->GetFixedParameters() );
    finalTransform->SetParameters( constellation2->GetVersorTransform()->GetParameters() );
    finalTransform->GetInverse( invFinalTransform );
    }
  else
    {
    ReaderType::Pointer atlasReader = ReaderType::New();
    atlasReader->SetFileName( this->m_atlasVolume );
    try
      {
      atlasReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Error while reading atlasVolume file:\n "
                << err << std::endl;
      }

    std::cout << "read atlas" << std::endl;
    // TODO: prob needs a try-catch
    LandmarksMapType referenceAtlasLandmarks = ReadSlicer3toITKLmk( this->m_atlasLandmarks );
    std::cout << "read atlas landmarks " << std::endl;

    LandmarksMapType acpcLandmarks = constellation2->GetAlignedPoints();

    // Now take invFinalTransform, and  IGNORE THE VERSION from line 253create a better version of invFinalTransform
    // using BRAINSFit.
    // take the the subjects landmarks in original space, and  landmarks from a reference Atlas, and compute an initial
    // affine transform
    // ( using logic from BRAINSLandmarkInitializer) and create initToAtlasAffineTransform.

    typedef itk::AffineTransform<double, Dimension> AffineTransformType;
    AffineTransformType::Pointer initToAtlasAffineTransform = AffineTransformType::New();

    typedef itk::LandmarkBasedTransformInitializer<AffineTransformType, ImageType,
                                                   ImageType> LandmarkBasedInitializerType;
    LandmarkBasedInitializerType::Pointer landmarkBasedInitializer = LandmarkBasedInitializerType::New();

    typedef std::map<std::string, float> WeightType;
    WeightType landmarkWeights;
    if( this->m_atlasLandmarkWeights != "" )
      {
      std::cout << "setting weights. " << std::endl;
      std::ifstream weightFileStream( this->m_atlasLandmarkWeights.c_str() );

      if( !weightFileStream.is_open() )
        {
        std::cerr << "Fail to open weight file " << std::endl;
        exit(EXIT_FAILURE);
        }

      std::string line;
      while( getline( weightFileStream, line ) )
        {
        const size_t      firstComma = line.find(',', 0);
        const std::string landmark = line.substr( 0, firstComma );
        const float       weight   = atof( (line.substr( firstComma + 1, line.length() - 1 ) ).c_str() );
        landmarkWeights[landmark] = weight;
        }
      }

    typedef LandmarkBasedInitializerType::LandmarkPointContainer LandmarkContainerType;
    LandmarkContainerType fixedLmks;
    LandmarkContainerType movingLmks;
    typedef LandmarksMapType::const_iterator LandmarkConstIterator;
    LandmarkBasedInitializerType::LandmarkWeightType landmarkWgts;
    for( LandmarkConstIterator fixedIt = referenceAtlasLandmarks.begin(); fixedIt != referenceAtlasLandmarks.end();
         ++fixedIt )
      {
      LandmarkConstIterator movingIt = acpcLandmarks.find( fixedIt->first );
      if( movingIt != acpcLandmarks.end() )
        {
        fixedLmks.push_back( fixedIt->second);
        movingLmks.push_back( movingIt->second);
        if( !this->m_atlasLandmarkWeights.empty() )
          {
          if( landmarkWeights.find( fixedIt->first ) != landmarkWeights.end() )
            {
            landmarkWgts.push_back( landmarkWeights[fixedIt->first] );
            }
          else
            {
            std::cout << "Landmark for " << fixedIt->first << " does not exist. "
                      << "Set the weight to 0.5 "
                      << std::endl;
            landmarkWgts.push_back( 0.5F );
            }
          }
        }
      else
        {
        std::cout << "i shouldnt be here" << std::endl;
        }
      }

    if( !this->m_atlasLandmarkWeights.empty() )
      {
      landmarkBasedInitializer->SetLandmarkWeight( landmarkWgts );
      }

    landmarkBasedInitializer->SetFixedLandmarks( fixedLmks );
    landmarkBasedInitializer->SetMovingLandmarks( movingLmks );
    landmarkBasedInitializer->SetTransform( initToAtlasAffineTransform );
    landmarkBasedInitializer->InitializeTransform();

    initToAtlasAffineTransform->Compose( constellation2->GetVersorTransform(), true );
    typedef itk::BRAINSFitHelper HelperType;
    HelperType::Pointer brainsFitHelper = HelperType::New();

    // Now Run BRAINSFitHelper class initialized with initToAtlasAffineTransform, original image, and atlas image
    // adapted from BRAINSABC/brainseg/AtlasRegistrationMethod.hxx - do I need to change any of these parameters?
    brainsFitHelper->SetNumberOfSamples(500000);
    brainsFitHelper->SetNumberOfHistogramBins(50);
    std::vector<int> numberOfIterations(1);
    numberOfIterations[0] = 1500;
    brainsFitHelper->SetNumberOfIterations(numberOfIterations);
    brainsFitHelper->SetTranslationScale(1000);
    brainsFitHelper->SetReproportionScale(1.0);
    brainsFitHelper->SetSkewScale(1.0);

    typedef itk::Image<float, 3>                            FloatImageType;
    typedef itk::CastImageFilter<ImageType, FloatImageType> CastFilterType;

    CastFilterType::Pointer fixedCastFilter = CastFilterType::New();
    fixedCastFilter->SetInput( atlasReader->GetOutput() );
    fixedCastFilter->Update();
    brainsFitHelper->SetFixedVolume( fixedCastFilter->GetOutput() );

    CastFilterType::Pointer movingCastFilter = CastFilterType::New();
    movingCastFilter->SetInput( reader->GetOutput() );
    movingCastFilter->Update();
    brainsFitHelper->SetMovingVolume( movingCastFilter->GetOutput() );

    std::vector<double> minimumStepSize(1);
    minimumStepSize[0] = 0.005;
    brainsFitHelper->SetMinimumStepLength(minimumStepSize);
    std::vector<std::string> transformType(1);
    transformType[0] = "Affine";
    brainsFitHelper->SetTransformType(transformType);

    brainsFitHelper->SetCurrentGenericTransform( initToAtlasAffineTransform.GetPointer() );
    brainsFitHelper->Update();
    // TODO: REMOVE ioguz WriteTransformToDisk( brainsFitHelper->GetCurrentGenericTransform()
    // ,"/Users/ioguz/forAli/brainsFitTransform.txt" );

    typedef VersorRigid3DTransformType VersorRigidTransformType;
    VersorRigidTransformType::Pointer brainsFitExtractedVersorRigid =
      itk::ComputeRigidTransformFromGeneric( brainsFitHelper->GetCurrentGenericTransform().GetPointer() );
    if( brainsFitExtractedVersorRigid.IsNull() )
      {
      // Fail if something weird happens.  TODO: This should throw an exception.
      std::cout << "brainsFitExtractedVersorRigid is null. It means we're not registering to the atlas, after all."
                << std::endl;
      std::cout << "FAILIING" << std::endl;
      exit(-1);
      }

    // as a final step, translate the AC back to the origin.
      {
      LandmarkConstIterator                           acIter = acpcLandmarks.find( "AC" );
      const VersorRigidTransformType::OutputPointType acOrigPoint =
        acpcTransformFromConstellation2->TransformPoint( acIter->second );

      finalTransform->SetFixedParameters( brainsFitExtractedVersorRigid->GetFixedParameters() );
      finalTransform->SetParameters( brainsFitExtractedVersorRigid->GetParameters() );
      finalTransform->GetInverse( invFinalTransform );

      // TODO:  CHECK if this can be less convoluted. Too many inverses used.  translate the forward by positive
      // rather than inverse by negative.
      //
      VersorRigidTransformType::OutputPointType acPoint = invFinalTransform->TransformPoint( acOrigPoint );
        {
        VersorRigidTransformType::OffsetType translation;
        translation[0] = -acPoint[0];
        translation[1] = -acPoint[1];
        translation[2] = -acPoint[2];
        invFinalTransform->Translate( translation, true );
        }
      invFinalTransform->GetInverse( finalTransform );

      // TODO: Remove VersorRigidTransformType::OutputPointType acFinalPoint =  invFinalTransform->TransformPoint (
      // acOrigPoint );
      }
    }

  // Landmark weights. All of them equal 1 wright now, but they should be calculated later
  LandmarksWeightMapType LandmarksWeightMap;
  double                 weights = 1;
  // Save landmarks in input/output or original/aligned space
  for( LandmarksMapType::const_iterator lit = constellation2->GetAlignedPoints().begin();
       lit != constellation2->GetAlignedPoints().end(); ++lit )
    {
    VersorTransformType::OutputPointType transformedPoint =
      acpcTransformFromConstellation2->TransformPoint( lit->second );
    this->m_outputLandmarksInInputSpaceMap[lit->first] = transformedPoint;
    if( this->m_atlasVolume.empty() )
      {
      this->m_outputLandmarksInACPCAlignedSpaceMap[lit->first] = lit->second;
      }
    else
      {
      // not final transform any more - rather, this is just the original acpc transform XXX
      this->m_outputLandmarksInACPCAlignedSpaceMap[lit->first] = invFinalTransform->TransformPoint( transformedPoint );
      }
    // or something like constellation2->GetOriginalPoints()[lit->first];
    LandmarksWeightMap[lit->first] = weights;
    }

  // ----------------------
  // Write results to disk
  std::cout << "\nWriting results to files..." << std::endl;
  if( this->m_outputTransform.compare( "" ) != 0 )
    {
    TransformWriterType::Pointer writer = TransformWriterType::New();
    writer->SetInput( finalTransform );
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

    short BackgroundFillValue;
    if( this->m_backgroundFillValueString == std::string("BIGNEG") )
      {
      BackgroundFillValue = -32768;
      }
    else
      {
      BackgroundFillValue = atoi( this->m_backgroundFillValueString.c_str() );
      }
    // the input image may have rescaled intensities, etc
    writer->SetInput(
      TransformResample<SImageType, SImageType>(
        constellation2->GetImageToBeResampled(),
        MakeIsoTropicReferenceImage(),
        BackgroundFillValue,
        GetInterpolatorFromString<SImageType>(this->m_interpolationMode),
        finalTransform.GetPointer() ) );
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

  // TODO: Due to writing a new Aligner, these flags are not needed anymore and they should be deleted
  if( ( this->m_outputLandmarksPaired.compare( "" ) != 0 ) &&
      ( this->m_inputLandmarksPaired.compare( "" ) != 0 ) )
    {
    LandmarksMapType                 origLandmarks = ReadSlicer3toITKLmk( this->m_inputLandmarksPaired );
    LandmarksMapType                 alignedLandmarks;
    LandmarksMapType::const_iterator olit = origLandmarks.begin();
    for( ; olit != origLandmarks.end(); ++olit )
      {
      alignedLandmarks[olit->first] = invFinalTransform->TransformPoint( olit->second );
      }
    WriteITKtoSlicer3Lmk( this->m_outputLandmarksPaired, alignedLandmarks );
    std::cout << "The acpc-aligned landmark list file is written." << std::endl;
    }
  else if( ( ( this->m_outputLandmarksPaired.compare( "" ) != 0 ) &&
             ( this->m_inputLandmarksPaired.compare( "" ) == 0 ) )
           ||
           ( ( this->m_outputLandmarksPaired.compare( "" ) == 0 ) &&
             ( this->m_inputLandmarksPaired.compare( "" ) != 0 ) ) )
    {
    std::cerr << "The outputLandmark parameter should be paired with"
              << "the inputLandmark parameter." << std::endl;
    std::cerr << "No output acpc-aligned landmark list file is generated" << std::endl;
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

  // TODO: Due to writing a new weight generator, this flag are not needed anymore and they should be deleted
  if( this->m_outputLandmarkWeights.compare( "" ) != 0 )
    {
    WriteITKtoSlicer3LmkWts( this->m_outputLandmarkWeights,
                             LandmarksWeightMap );
    std::cout << "The output landmark weights list file is written." << std::endl;
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
                   finalTransform );
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
