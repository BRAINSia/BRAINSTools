/*
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "BRAINSConstellationDetectorPrimary.h"
#include "BRAINSConstellationDetectorCLP.h"

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

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

  // set the template model to default
  if( inputTemplateModel.compare( "" ) == 0 )
    {
    std::string pathOut;
    std::string errorMsg;

    if( !itksys::SystemTools::FindProgramPath( argv[0], pathOut, errorMsg) )

      {
      std::cerr << "Error: Input Template Model File not found" << std::endl;
      std::cerr << errorMsg << std::endl;

      return 1;
      }

    inputTemplateModel = itksys::SystemTools::GetProgramPath( pathOut.c_str() ) + "/" + "T1.mdl";
    std::cout << "Set inputTemplateModel to default: " << std::endl;
    std::cout << inputTemplateModel << std::endl;
    }

  // set the llsModel to default
  if( llsModel == "" )
    {
    std::string pathOut;
    std::string errorMsg;

    if( !itksys::SystemTools::FindProgramPath( argv[0], pathOut, errorMsg) )

      {
      std::cerr << "Error: Input LLSModel File not found" << std::endl;
      std::cerr << errorMsg << std::endl;

      return 1;
      }

    llsModel = itksys::SystemTools::GetProgramPath( pathOut.c_str() ) + "/" + "LLSModel.h5";
    std::cout << "Set LLSModel to default: " << std::endl;
    std::cout << llsModel << std::endl;
    }

  BRAINSConstellationDetectorPrimary BCD;

  BCD.SetNumberOfThreads( numberOfThreads );
  BCD.SetInputLandmarksEMSP( inputLandmarksEMSP );
  BCD.SetLLSModel( llsModel );
  BCD.SetInputVolume( inputVolume );
  BCD.SetForceHoughEyeDetectorReportFailure( forceHoughEyeDetectorReportFailure );
  BCD.SetHoughEyeDetectorMode( houghEyeDetectorMode );
  BCD.SetResultsDir( resultsDir );
  BCD.SetWritedebuggingImagesLevel( writedebuggingImagesLevel );
  BCD.SetInputTemplateModel( inputTemplateModel );
  BCD.SetMspQualityLevel( mspQualityLevel );
  BCD.SetOtsuPercentileThreshold( otsuPercentileThreshold );
  BCD.SetAcLowerBound( acLowerBound );
  BCD.SetCutOutHeadInOutputVolume( cutOutHeadInOutputVolume );
  BCD.SetRescaleIntensities( rescaleIntensities );
  BCD.SetTrimRescaledIntensities( trimRescaledIntensities );
  BCD.SetRescaleIntensitiesOutputRange( rescaleIntensitiesOutputRange );
  BCD.SetBackgroundFillValueString( backgroundFillValueString );
  BCD.SetInterpolationMode( interpolationMode );
  BCD.SetForceACPoint( forceACPoint );
  BCD.SetForcePCPoint( forcePCPoint );
  BCD.SetForceVN4Point( forceVN4Point );
  BCD.SetForceRPPoint( forceRPPoint );
  BCD.SetRadiusMPJ( radiusMPJ );
  BCD.SetRadiusAC( radiusAC );
  BCD.SetRadiusPC( radiusPC );
  BCD.SetRadiusVN4( radiusVN4 );
  BCD.SetDebug( debug );
  BCD.SetVerbose( verbose );
  BCD.SetWriteBranded2DImage( writeBranded2DImage );
  BCD.SetOutputTransform( outputTransform );
  BCD.SetOutputVolume( outputVolume );
  BCD.SetOutputResampledVolume( outputResampledVolume );
  BCD.SetInputLandmarksPaired( inputLandmarksPaired );
  BCD.SetOutputLandmarksPaired( outputLandmarksPaired );
  BCD.SetOutputLandmarksInInputSpace( outputLandmarksInInputSpace );
  BCD.SetOutputLandmarksInACPCAlignedSpace( outputLandmarksInACPCAlignedSpace );
  BCD.SetOutputLandmarkWeights( outputLandmarkWeights );
  BCD.SetOutputMRML( outputMRML );
  BCD.SetOutputVerificationScript( outputVerificationScript );
  BCD.SetOutputUntransformedClippedVolume( outputUntransformedClippedVolume );

  try
    {
    BCD.Compute();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception Object caught:\n" << err << std::endl;
    return EXIT_FAILURE;
    }

  return 0;
}
