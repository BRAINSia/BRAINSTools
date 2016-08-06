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
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "BRAINSConstellationDetectorPrimary.h"
#include "BRAINSConstellationDetectorCLP.h"

#include "BRAINSConstellationDetectorVersion.h"

int main( int argc, char *argv[] )
{
  PARSE_ARGS;
  FFTWInit(""); //Initialize for FFTW in order to improve performance of subsequent runs
  BRAINSRegisterAlternateIO();

  const std::string Version(BCDVersionString);
  std::cout << "Run BRAINSConstellationDetector Version: " << Version << std::endl;
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
      && ( outputUntransformedClippedVolume.compare( "" ) == 0 ) )
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
  BCD.SetAtlasVolume( atlasVolume );
  BCD.SetAtlasLandmarks( atlasLandmarks );
  BCD.SetAtlasLandmarkWeights( atlasLandmarkWeights );
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
  BCD.SetOutputLandmarksInInputSpace( outputLandmarksInInputSpace );
  BCD.SetOutputLandmarksInACPCAlignedSpace( outputLandmarksInACPCAlignedSpace );
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

    // Write a log script to the disk to report failure
    std::stringstream failureLogFileStream;
    failureLogFileStream << err;
    std::ofstream failureLogScript;
    failureLogScript.open("BCD_FAILED.txt");
    if( !failureLogScript.is_open() )
      {
      std::cerr << "Error: Can't write failure log file: BCD_FAILED.txt " << std::endl;
      std::cerr.flush();
      }
    failureLogScript << failureLogFileStream.str();
    failureLogScript.close();

    return EXIT_FAILURE;
    }

  return 0;
}
