/*
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "itksys/SystemTools.hxx"
#include "BRAINSThreadControl.h"

#include "BRAINSConstellationDetectorCLP.h"
#include "itkFindCenterOfBrainFilter.h"
#include "BRAINSHoughEyeDetector.h"
#include "BRAINSConstellationDetector2.h"
#include "Slicer3LandmarkIO.h"
#include "Slicer3LandmarkWeightIO.h"
#include "LLSModel.h"

#include "itkCommand.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"
#include "itkVersorRigid3DTransform.h"

#include "stdio.h"
#include <iostream>
#include <fstream>
#include "string.h"
#include <map>

// Image, filter, transform typedef
const unsigned int LocalImageDimension = 3;
typedef short PixelType;

typedef itk::Image<PixelType, LocalImageDimension> ImageType;
typedef ImageType::Pointer                         ImagePointerType;
typedef ImageType::PointType                       ImagePointType;
typedef ImageType::SpacingType                     ImageSpacingType;
typedef ImageType::SizeType                        ImageSizeType;
typedef ImageType::DirectionType                   ImageDirectionType;
typedef ImageType::IndexType                       ImageIndexType;
typedef std::map<std::string, ImagePointType>      LandmarksMapType;
typedef std::map<std::string, double>              LandmarksWeightMapType; // -- Added by Ali

typedef itk::ImageFileReader<ImageType>         ReaderType;
typedef itk::ImageFileWriter<ImageType>         WriterType;
typedef itk::FindCenterOfBrainFilter<ImageType> FindCenterFilter;
typedef itk::BRAINSHoughEyeDetector<
    ImageType, ImageType>                             HoughEyeDetectorType;
typedef itk::BRAINSConstellationDetector2<
    ImageType, ImageType>                             Constellation2Type;
typedef itk::TransformFileWriter            TransformWriterType;
typedef itk::VersorRigid3DTransform<double> VersorTransformType;

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  BRAINSUtils::SetThreadCount(numberOfThreads);

  // ------------------------------------
  // Verify input parameters
  std::cout << "Verifying input parameters..." << std::endl;
  if( inputVolume.compare( "" ) == 0 )
    {
    std::cerr << "To run the program please specify the input volume filename." << std::endl;
    std::cerr << "Type " << argv[0] << " -h for more help." << std::endl;
    return EXIT_FAILURE;
    }

  // Get a warning if none of the main output filename is specified
  if( ( outputVolume.compare( "" ) == 0 )
      && ( outputResampledVolume.compare( "" ) == 0 )
      && ( outputTransform.compare( "" ) == 0 )
      && ( outputLandmarksInACPCAlignedSpace.compare( "" ) == 0 )
      && ( outputLandmarksInInputSpace.compare( "" ) == 0 )
      && ( outputUntransformedClippedVolume.compare( "" ) == 0 )
      && ( ( inputLandmarksPaired.compare( "" ) == 0 ) ||
           ( outputLandmarksPaired.compare( "" ) == 0 ) ) )
    {
    std::cout << "WARNING: None of the main output filename is specified!" << std::endl;
    std::cout << "Try to specify at least one of the following output filenames:" << std::endl;
    std::cout << "outputVolume" << std::endl;
    std::cout << "outputResampledVolume" << std::endl;
    std::cout << "outputTransform" << std::endl;
    std::cout << "outputLandmarksInACPCAlignedSpace" << std::endl;
    std::cout << "outputLandmarksInInputSpace" << std::endl;
    std::cout << "outputUntransformedClippedVolume\n" << std::endl;
    }

  // set to default
  if( inputTemplateModel.compare( "T1.mdl" ) == 0 )
    {
    std::string pathOut;
    std::string errorMsg;

    // itksys_stl::string errorMessage;
    if( !itksys::SystemTools::FindProgramPath( argv[0], pathOut, errorMsg) )

      {
      std::cerr << "Error: Input Model File not found" << std::endl;
      std::cerr << errorMsg << std::endl;

      return 1;
      }
    // std::cerr << "ERROR: Could not find " << argv[0] << ": " << errorMessage
    // << std::endl;

    inputTemplateModel = itksys::SystemTools::GetProgramPath( pathOut.c_str() ) + "/" + inputTemplateModel;
    std::cout << "Set inputTemplateModel to default: " << std::endl;
    std::cout << inputTemplateModel << std::endl;
    }

  // ------------------------------------
  // Read external files
  std::cout << "\nReading in external files..." << std::endl;

  // load corresponding landmarks in EMSP aligned space from file if possible
  LandmarksMapType landmarksEMSP;
  if( inputLandmarksEMSP.compare( "" ) != 0 )
    {
    landmarksEMSP = ReadSlicer3toITKLmk( inputLandmarksEMSP );
    }

  // read in lls model file
  std::map<std::string, std::vector<double> > llsMeans;
  std::map<std::string, MatrixType>           llsMatrices;
  std::map<std::string, double>               searchRadii;

  if( llsModel == "" )
    {
    std::cerr << "Missing LLSModel file" << std::endl;
    return EXIT_FAILURE;
    }
  LLSModel theModel;
  theModel.SetFileName(llsModel);
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
  reader->SetFileName( inputVolume );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << " Error while reading image file( s ) with ITK:\n "
              << err << std::endl;
    }
  std::cout << "Processing: " << inputVolume << std::endl;

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
  else if( forceHoughEyeDetectorReportFailure == true )
    {
    houghEyeDetector->SetFailure( true );
    std::cout << "\nThe Hough eye detector is doomed to failure as notified." << std::endl;
    std::cout << "Skip estimation steps for eye centers." << std::endl;
    }
  else
    {
    std::cout << "\nFinding eye centers with BRAINS Hough Eye Detector..." << std::endl;
    houghEyeDetector->SetInput( reader->GetOutput() );
    houghEyeDetector->SetHoughEyeDetectorMode( houghEyeDetectorMode );
    houghEyeDetector->SetResultsDir( resultsDir );       // debug output dir
    houghEyeDetector->SetWritedebuggingImagesLevel( writedebuggingImagesLevel );
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
    if( landmarksEMSP.size() > 0 )
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
      constellation2->SetHoughEyeTransform( houghEyeDetector->GetVersorTransform() );
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
  constellation2->SetInputTemplateModel( inputTemplateModel );
  constellation2->SetMspQualityLevel( mspQualityLevel );
  constellation2->SetOtsuPercentileThreshold( otsuPercentileThreshold );
  constellation2->SetAcLowerBound( acLowerBound );
  constellation2->SetCutOutHeadInOutputVolume( cutOutHeadInOutputVolume );
  constellation2->SetRescaleIntensities( rescaleIntensities );
  constellation2->SetTrimRescaledIntensities( trimRescaledIntensities );
  constellation2->SetRescaleIntensitiesOutputRange( rescaleIntensitiesOutputRange );
  constellation2->SetBackgroundFillValueString( backgroundFillValueString );
  constellation2->SetInterpolationMode( interpolationMode );
  constellation2->SetForceACPoint( forceACPoint );   // In original space
  constellation2->SetForcePCPoint( forcePCPoint );
  constellation2->SetForceVN4Point( forceVN4Point );
  constellation2->SetForceRPPoint( forceRPPoint );
  constellation2->SetRadiusMPJ( radiusMPJ );
  constellation2->SetRadiusAC( radiusAC );
  constellation2->SetRadiusPC( radiusPC );
  constellation2->SetRadiusVN4( radiusVN4 );
  constellation2->SetDebug( debug );
  constellation2->SetVerbose( verbose );
  constellation2->SetWritedebuggingImagesLevel( writedebuggingImagesLevel );
  constellation2->SetWriteBranded2DImage( writeBranded2DImage );
  constellation2->SetResultsDir( resultsDir );
  constellation2->SetLlsMatrices( llsMatrices );
  constellation2->SetLlsMeans( llsMeans );
  constellation2->SetSearchRadii( searchRadii );
  constellation2->SetOriginalInputImage( reader->GetOutput() );
  constellation2->Update();

  // Get the final transform
  VersorTransformType::Pointer finalTransform = VersorTransformType::New();
  VersorTransformType::Pointer invFinalTransform = VersorTransformType::New();
  finalTransform->SetParameters( constellation2->GetVersorTransform()->GetParameters() );
  finalTransform->GetInverse( invFinalTransform );

  // Landmark weights. All of them equal 1 wright now, but they should be calculated later
  LandmarksWeightMapType LandmarksWeightMap;
  double                 weights = 1;

  // Save landmarks in input/output or original/aligned space
  LandmarksMapType                 outputLandmarksInInputSpaceMap;
  LandmarksMapType                 outputLandmarksInACPCAlignedSpaceMap;
  LandmarksMapType::const_iterator lit =
    constellation2->GetAlignedPoints().begin();
  for( ; lit != constellation2->GetAlignedPoints().end(); ++lit )
    {
    outputLandmarksInACPCAlignedSpaceMap[lit->first] = lit->second;
    outputLandmarksInInputSpaceMap[lit->first] =
      finalTransform->TransformPoint( lit->second );
    // or something like constellation2->GetOriginalPoints()[lit->first];

    LandmarksWeightMap[lit->first] = weights;
    }

  // ----------------------
  // Write results to disk
  std::cout << "\nWriting results to files..." << std::endl;
  if( outputTransform.compare( "" ) != 0 )
    {
    TransformWriterType::Pointer writer = TransformWriterType::New();
    writer->SetInput( finalTransform );
    writer->SetFileName( outputTransform );
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
  if( outputVolume.compare( "" ) != 0 )
    {
    preferedOutputReferenceImage = outputVolume;
    // This will be overwritten if outputResampledVolume is set
    // Write the aligned image to a file
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputVolume );
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

  if( outputResampledVolume.compare( "" ) != 0 )
    {
    preferedOutputReferenceImage = outputResampledVolume;
    // This will be overwritten if outputResampledVolume is set
    // Write the aligned image to a file
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputResampledVolume );
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

  if( ( outputLandmarksPaired.compare( "" ) != 0 ) &&
      ( inputLandmarksPaired.compare( "" ) != 0 ) )
    {
    LandmarksMapType                 origLandmarks = ReadSlicer3toITKLmk( inputLandmarksPaired );
    LandmarksMapType                 alignedLandmarks;
    LandmarksMapType::const_iterator olit = origLandmarks.begin();
    for( ; olit != origLandmarks.end(); ++olit )
      {
      alignedLandmarks[olit->first] = invFinalTransform->TransformPoint( olit->second );
      }
    WriteITKtoSlicer3Lmk( outputLandmarksPaired, alignedLandmarks );
    std::cout << "The acpc-aligned landmark list file is written." << std::endl;
    }
  else if( ( ( outputLandmarksPaired.compare( "" ) != 0 ) &&
             ( inputLandmarksPaired.compare( "" ) == 0 ) )
           ||
           ( ( outputLandmarksPaired.compare( "" ) == 0 ) &&
             ( inputLandmarksPaired.compare( "" ) != 0 ) ) )
    {
    std::cerr << "The outputLandmark parameter should be paired with"
              << "the inputLandmark parameter." << std::endl;
    std::cerr << "No output acpc-aligned landmark list file is generated" << std::endl;
    std::cerr << "Type " << argv[0] << " -h for more help." << std::endl;
    }

  if( outputLandmarksInInputSpace.compare( "" ) != 0 )
    {
    WriteITKtoSlicer3Lmk( outputLandmarksInInputSpace,
                          outputLandmarksInInputSpaceMap );
    std::cout << "The output landmarks list file in the original space is written." << std::endl;
    }

  if( outputLandmarksInACPCAlignedSpace.compare( "" ) != 0 )
    {
    WriteITKtoSlicer3Lmk( outputLandmarksInACPCAlignedSpace,
                          outputLandmarksInACPCAlignedSpaceMap );
    std::cout << "The output landmarks list file in the output space is written." << std::endl;
    if( preferedOutputReferenceImage.compare( "" ) == 0 )
      {
      std::cout << "WARNING no aligned output volume is requested." << std::endl;
      }
    }

  if( outputLandmarkWeights.compare( "" ) != 0 )
    {
    WriteITKtoSlicer3LmkWts( outputLandmarkWeights,
                             LandmarksWeightMap );
    std::cout << "The output landmark weights list file is written." << std::endl;
    }

  if( outputMRML.compare( "" ) != 0 )
    {
    WriteMRMLFile( outputMRML,
                   outputLandmarksInInputSpace,
                   outputLandmarksInACPCAlignedSpace,
                   inputVolume,
                   preferedOutputReferenceImage,
                   outputTransform,
                   outputLandmarksInInputSpaceMap,
                   outputLandmarksInACPCAlignedSpaceMap,
                   finalTransform );
    std::cout << "The output mrml scene file is written." << std::endl;
    }

  if( outputUntransformedClippedVolume.compare( "" ) != 0 )
    {
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputUntransformedClippedVolume );
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

  if( ( outputVerificationScript.compare( "" ) != 0 )
      && ( outputLandmarksInACPCAlignedSpace.compare( "" ) != 0 )
      && ( outputVolume.compare( "" ) != 0 ) )
    {
    writeVerificationScript( outputVerificationScript, outputVolume, outputLandmarksInACPCAlignedSpace );
    std::cout << "The verification script is written." << std::endl;
    }

  return 0;
}
