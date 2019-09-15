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
#include "BRAINSDemonWarpCommonLibWin32Header.h"
#include "BRAINSCommonLibWin32Header.h"
#include "BRAINSThreadControl.h"
#include "BRAINSDemonWarpCLP.h"

#include "BRAINSDemonWarpTemplates.h"

#include "DebugImageWrite.h"
DEFINE_DEBUG_IMAGE_COUNTER;

#ifdef USE_DebugImageViewer
/*************************
 * Have a global variable to
 * add debugging information.
 */
DebugImageViewerClient DebugImageDisplaySender;
#endif

int
main( int argc, char * argv[] )
{
  struct BRAINSDemonWarpAppParameters command;

  {
    PARSE_ARGS;
    BRAINSRegisterAlternateIO();
    const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder( numberOfThreads );

#ifdef USE_DebugImageViewer
    DebugImageDisplaySender.SetEnabled( UseDebugImageViewer );
    DebugImageDisplaySender.SetPromptUser( PromptAfterImageSend );
#endif

    command.registrationFilterType = registrationFilterType;
    command.movingVolume = movingVolume;
    command.fixedVolume = fixedVolume;
    command.outputVolume = outputVolume;
    command.outputDisplacementFieldVolume = outputDisplacementFieldVolume;
    command.inputPixelType = inputPixelType;
    command.outputPixelType = outputPixelType;
    command.outputDisplacementFieldPrefix = outputDisplacementFieldPrefix;
    command.outputCheckerboardVolume = outputCheckerboardVolume;
    command.outputNormalized = outputNormalized;
    command.outputDebug = outputDebug;
    command.maskProcessingMode = maskProcessingMode;
    command.fixedBinaryVolume = fixedBinaryVolume;
    command.movingBinaryVolume = movingBinaryVolume;
    command.lowerThresholdForBOBF = lowerThresholdForBOBF;
    command.upperThresholdForBOBF = upperThresholdForBOBF;
    command.backgroundFillValue = backgroundFillValue;
    // Not yet implemented.
    // command.forceCoronalZeroOrigin = forceCoronalZeroOrigin;
    //    command.movingLandmarks = movingLandmarks;
    //    command.fixedLandmarks = fixedLandmarks;
    //    command.initializeWithFourier = initializeWithFourier;
    command.initializeWithDisplacementField = initializeWithDisplacementField;
    command.initializeWithTransform = initializeWithTransform;

    command.histogramMatch = histogramMatch;
    command.numberOfHistogramLevels = numberOfHistogramBins;
    command.numberOfMatchPoints = numberOfMatchPoints;
    command.numberOfLevels = numberOfPyramidLevels;
    command.numberOfIterations.SetSize( numberOfPyramidLevels );

    command.maxStepLength = maxStepLength;
    command.gradientType = gradientType;
    command.smoothDisplacementFieldSigma = smoothDisplacementFieldSigma;
    command.smoothingUp = smoothingUp;
    command.numberOfBCHApproximationTerms = numberOfBCHApproximationTerms;
    command.interpolationMode = interpolationMode;
    for ( int i = 0; i < numberOfPyramidLevels; i++ )
    {
      command.numberOfIterations[i] = arrayOfPyramidLevelIterations[i];
    }
    for ( int i = 0; i < 3; i++ )
    {
      command.theMovingImageShrinkFactors[i] = minimumMovingPyramid[i];
      command.theFixedImageShrinkFactors[i] = minimumFixedPyramid[i];
    }
    for ( int i = 0; i < 3; i++ )
    {
      command.checkerboardPatternSubdivisions[i] = checkerboardPatternSubdivisions[i];
      command.seedForBOBF[i] = seedForBOBF[i];
      command.neighborhoodForBOBF[i] = neighborhoodForBOBF[i];
      command.medianFilterSize[i] = medianFilterSize[i];
    }
  }

  //  bool debug=true;
  if ( command.outputDebug )
  {
    std::cout << "                   movingVolume: " << command.movingVolume << std::endl
              << "                    fixedVolume: " << command.fixedVolume << std::endl
              << "                   outputVolume: " << command.outputVolume << std::endl
              << "   outputDisplacementFieldVolume: " << command.outputDisplacementFieldVolume << std::endl
              << "                 inputPixelType: " << command.inputPixelType << std::endl
              << "                outputPixelType: " << command.outputPixelType << std::endl
              << "  outputDisplacementFieldPrefix: " << command.outputDisplacementFieldPrefix << std::endl
              << "       outputCheckerboardVolume: " << command.outputCheckerboardVolume << std::endl
              << "checkerboardPatternSubdivisions: " << command.checkerboardPatternSubdivisions << std::endl
              << "               outputNormalized: " << command.outputNormalized << std::endl
              << "                    outputDebug: " << command.outputDebug << std::endl
              << "             maskProcessingMode: " << command.maskProcessingMode << std::endl
              << "              fixedBinaryVolume: " << command.fixedBinaryVolume << std::endl
              << "             movingBinaryVolume: " << command.movingBinaryVolume << std::endl
              << "          lowerThresholdForBOBF: " << command.lowerThresholdForBOBF << std::endl
              << "          upperThresholdForBOBF: " << command.upperThresholdForBOBF << std::endl
              << "            backgroundFillValue: " << command.backgroundFillValue << std::endl
              << "                    seedForBOBF: " << command.seedForBOBF << std::endl
              << "            neighborhoodForBOBF: " << command.neighborhoodForBOBF << std::endl
              << "               medianFilterSize: " << command.medianFilterSize
              << std::endl
              /** NOT YET IMPLEMENTED
               *  << "        movingLandmarks: " << command.movingLandmarks << std::endl
               *  << "         fixedLandmarks: " << command.fixedLandmarks << std::endl
               *  << "     initializeWithFourier: " << command.initializeWithFourier
               *  << std::endl
               */
              << "  initializeWithDisplacementField: " << command.initializeWithDisplacementField << std::endl
              << "       initializeWithTransform: " << command.initializeWithTransform << std::endl
              << "                    gradientType: " << command.gradientType << std::endl
              << "                   maxStepLength: " << command.maxStepLength << std::endl
              << "     smoothDisplacementFieldSigma: " << command.smoothDisplacementFieldSigma << std::endl
              << "                     smoothingUp: " << command.smoothingUp << std::endl
              << "                   histogramMatch: " << command.histogramMatch << std::endl
              << "                histogram levels: " << command.numberOfHistogramLevels << std::endl
              << "                 matching points: " << command.numberOfMatchPoints << std::endl
              << "   numberOfBCHApproximationTerms: " << command.numberOfBCHApproximationTerms << std::endl;
  }

  bool violated = false;
  if ( command.movingVolume.size() == 0 )
  {
    violated = true;
    std::cout << "  --movingVolume Required! " << std::endl;
  }
  if ( command.fixedVolume.size() == 0 )
  {
    violated = true;
    std::cout << "  --fixedVolume Required! " << std::endl;
  }

  if ( ( command.checkerboardPatternSubdivisions[0] == 0 ) || ( command.checkerboardPatternSubdivisions[1] == 0 ) ||
       ( command.checkerboardPatternSubdivisions[2] == 0 ) )
  {
    std::cout << "Invalid Patameters. The value of checkboardPatternSubdivisions should not be zero!" << std::endl;
    return EXIT_FAILURE;
  }

  if ( violated )
  {
    return EXIT_FAILURE;
  }

  // Test if the input data type is valid
  if ( command.inputPixelType != "" )
  {
    // check to see if valid type
    if ( ( CompareNoCase( command.inputPixelType, std::string( "uchar" ) ) ) &&
         ( CompareNoCase( command.inputPixelType, std::string( "short" ) ) ) &&
         ( CompareNoCase( command.inputPixelType, std::string( "ushort" ) ) ) &&
         ( CompareNoCase( command.inputPixelType, std::string( "int" ) ) ) &&
         ( CompareNoCase( command.inputPixelType, std::string( "float" ) ) )
#ifdef _USE_UNCOMMON_TYPES // This is commented out because it causes too many
                           // segments in one object file for the intel
                           // compiler
         && ( CompareNoCase( command.inputPixelType, std::string( "uint" ) ) ) &&
         ( CompareNoCase( command.inputPixelType, std::string( "double" ) ) )
#endif
    )
    {
      std::cout << "Error. Invalid data type string specified with --inputPixelType!" << std::endl;
      std::cout << "Use one of the following:" << std::endl;
      PrintDataTypeStrings();
      return EXIT_FAILURE;
    }
  }

  if ( command.outputPixelType != "" )
  {
    // check to see if valid type
    if ( ( CompareNoCase( command.outputPixelType, std::string( "uchar" ) ) ) &&
         ( CompareNoCase( command.outputPixelType, std::string( "SHORT" ) ) ) &&
         ( CompareNoCase( command.outputPixelType, std::string( "ushort" ) ) ) &&
         ( CompareNoCase( command.outputPixelType, std::string( "int" ) ) ) &&
         ( CompareNoCase( command.outputPixelType, std::string( "float" ) ) )
#ifdef _USE_UNCOMMON_TYPES // This is commented out because it causes too many
                           // segments in one object file for the intel
                           // compiler
         && ( CompareNoCase( command.outputPixelType, std::string( "uint" ) ) ) &&
         ( CompareNoCase( command.outputPixelType, std::string( "double" ) ) )
#endif
    )
    {
      std::cout << "Error. Invalid data type string specified with --outputPixelType!" << std::endl;
      std::cout << "Use one of the following:" << std::endl;
      PrintDataTypeStrings();
      return EXIT_FAILURE;
    }
  }

  // Call the process output data type function based on the input data type.
  if ( CompareNoCase( command.inputPixelType, std::string( "uchar" ) ) == 0 )
  {
    ProcessOutputType_uchar( command );
  }
  else if ( CompareNoCase( command.inputPixelType, std::string( "short" ) ) == 0 )
  {
    ProcessOutputType_short( command );
  }
  else if ( CompareNoCase( command.inputPixelType, std::string( "ushort" ) ) == 0 )
  {
    ProcessOutputType_ushort( command );
  }
  else if ( CompareNoCase( command.inputPixelType, std::string( "int" ) ) == 0 )
  {
    ProcessOutputType_int( command );
  }
  else if ( CompareNoCase( command.inputPixelType, std::string( "float" ) ) == 0 )
  {
    ProcessOutputType_float( command );
  }
#ifdef _USE_UNCOMMON_TYPES // This is commented out because it causes too many
                           // segments in one object file for the intel
                           // compiler
  else if ( CompareNoCase( command.inputPixelType, std::string( "uint" ) ) == 0 )
  {
    ProcessOutputType_uint( command );
  }
  else if ( CompareNoCase( command.inputPixelType, std::string( "double" ) ) == 0 )
  {
    ProcessOutputType_double( command );
  }
#endif
  else
  {
    std::cout << "Error. Invalid data type for --inputPixelType!  Use one of these:" << std::endl;
    PrintDataTypeStrings();
    return EXIT_FAILURE;
  }
  return 0;
}
