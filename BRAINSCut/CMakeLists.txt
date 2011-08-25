project(BRAINSCut)
set(LOCAL_PROJECT_NAME BRAINSCut)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

enable_testing()
#include(Dart)
include(CTest)
include(CPack)
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBuildMacros.cmake)
include(${BRAINS_CMAKE_HELPER_DIR}/PreventInSourceBuilds.cmake)
set(BRAINSLogo ${BRAINS_CMAKE_HELPER_DIR}/BRAINSLogo.h CACHE
FILEPATH "Path to BRAINS Logo")


find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})
message(STATUS "BRAINSCut OpenCV_DIR = ${OpenCV_DIR}")
find_package( OpenCV REQUIRED)
#include( ${OpenCV_USE_FILE} )

#-----------------------------------------------------------------------------
# Output directories.
#
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBRAINS3BuildMacros.cmake)
#SETOPTIONALDEBUGIMAGEVIEWER()


# No need to add with nothing to compile. add_subdirectory(BRAINSFit_Common)
include(${BRAINS_CMAKE_HELPER_DIR}/IJMacros.txt)

###
include_directories(${BRAINSCut_SOURCE_DIR}/ProcessDescription)
add_subdirectory(ProcessDescription)
add_subdirectory(TrainedModels)
add_subdirectory(FeatureCreators)
add_subdirectory(Tools)
add_subdirectory(RegistrationScripts)
file(GLOB BRAINSCut_HDRS *.h *.txx)

add_library(BRAINSCutCommon
            Utilities.cxx
            GenerateRegistrations.cxx
            XMLParser.cxx
            Parser.cxx
            Svm.cxx
            ShuffleVectors.cxx
            ${BRAINSCut_HDRS}
            )
### OPENCV
target_link_libraries( BRAINSCutCommon
                       ProcessDescription
                       BRAINSCommonLib
                       ${ITK_LIBRARIES}
                       ${OpenCV_LIBS}
                       )


set(CLP BRAINSCut)
set(${CLP}_SOURCE BRAINSCut.cxx BRAINSCut.cxx
  GenerateProbability.cxx
  CreateVectors.cxx
  TrainModel.cxx
  ApplyModel.cxx
  )
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSCut.xml "${${CLP}_SOURCE}" "" "main" "BRAINSCutCommon")
target_link_libraries(BRAINSCut BRAINSCutCommon BRAINSCommonLib)

set(CLP GenerateSummedGradientImage)
set(${CLP}_SOURCE GenerateSummedGradientImage.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} GenerateSummedGradientImage.xml "${${CLP}_SOURCE}" "" "main" "BRAINSCutCommon")
target_link_libraries(GenerateSummedGradientImage BRAINSCutCommon)

set(CLP GenerateBrainClippedImage)
set(${CLP}_SOURCE GenerateBrainClippedImage.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} GenerateBrainClippedImage.xml "${${CLP}_SOURCE}" "" "main" "BRAINSCutCommon")
target_link_libraries(GenerateBrainClippedImage BRAINSCutCommon)

add_executable(ShuffleVectorsModule ShuffleVectorsModule.cxx)
target_link_libraries(ShuffleVectorsModule BRAINSCutCommon)

add_executable(DumpBinTrainingVectors DumpBinTrainingVectors.cxx)
target_link_libraries(DumpBinTrainingVectors BRAINSCutCommon)

add_executable(GenerateCsfClippedFromClassifiedImage GenerateCsfClippedFromClassifiedImage.cxx)
target_link_libraries(GenerateCsfClippedFromClassifiedImage
  BRAINSCutCommon)

add_executable(WriteFeaturesforROI WriteFeaturesforROI.cxx)
target_link_libraries(WriteFeaturesforROI BRAINSCutCommon)

install(TARGETS BRAINSCut
                GenerateCsfClippedFromClassifiedImage
                GenerateSummedGradientImage
                GenerateBrainClippedImage
                DumpBinTrainingVectors
                ShuffleVectorsModule
                BRAINSCutCommon
                WriteFeaturesforROI
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib/static)


# --------------------------------------------------------- #
# BRAINSCut Testing
# --------------------------------------------------------- #
##
#if ( 0 ) ## HACK:  Just silencing failing test until they can be properly addressed

enable_testing()

## BRAINSCut NetConfiguration(XML) file for testing
##
set( InputNetConfigurationForTesting ${BRAINSCut_SOURCE_DIR}/TestSuite/NetConfigurations/BRAINSCutTest.xml.in)
set( OutputNetConfigurationForTesting ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCutTest.xml )

configure_file( ${InputNetConfigurationForTesting} ${OutputNetConfigurationForTesting}  IMMEDIATE @ONLY)

##
## BRAINSCut Test Run
add_test( GenerateProbabilityTest $<TARGET_FILE:BRAINSCut>
                                  --generateProbability --netConfiguration ${OutputNetConfigurationForTesting})
add_test( CreateVectorsTest $<TARGET_FILE:BRAINSCut>
                                  --createVectors --netConfiguration ${OutputNetConfigurationForTesting})
add_test( TrainModelTest $<TARGET_FILE:BRAINSCut>
                                  --trainModel --netConfiguration ${OutputNetConfigurationForTesting})
add_test( ApplyModelTest $<TARGET_FILE:BRAINSCut>
                                  --applyModel --netConfiguration ${OutputNetConfigurationForTesting})
add_test( AllInOneTest $<TARGET_FILE:BRAINSCut>
                                  --generateProbability  --createVectors  --trainModel  --applyModel
                                  --netConfiguration ${OutputNetConfigurationForTesting})
#endif()









