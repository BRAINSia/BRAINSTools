
##- project(BRAINSCut)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil( BRAINSCut_ITK
  ITKAnisotropicSmoothing
  ITKBinaryMathematicalMorphology
  ITKConnectedComponents
  ITKDistanceMap
  ITKFastMarching
  ITKImageCompare
  ITKImageFeature
  ITKImageGradient
  ITKImageLabel
  ITKLabelVoting
  ITKLevelSets
  ITKMathematicalMorphology
  ITKMetricsv4
  ITKOptimizersv4
  ITKRegionGrowing
  ITKRegistrationCommon
  ITKRegistrationMethodsv4
  ITKReview
  ITKSmoothing
  ITKTestKernel
  ITKThresholding
)

#
# OpenCV
#
message(STATUS "BRAINSCut OpenCV_DIR = ${OpenCV_DIR}")
find_package( OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

set(MODULE_FOLDER "Module-BRAINSCut")

###
add_subdirectory(BRAINSFeatureCreators)

include_directories(${BRAINSTools_SOURCE_DIR}/BRAINSCut/BRAINSCutConfiguration)
include_directories(${BRAINSTools_SOURCE_DIR}/BRAINSCut)
add_subdirectory(BRAINSCutConfiguration)


file(GLOB BRAINSCut_HDRS *.h *.hxx)

set(BRAINSCutCOMMONLIB_SRCS
  XMLParser.cxx
  XMLConfigurationFileParser.cxx
  ShuffleVectors.cxx
  BRAINSCutDataHandler.cxx
  BRAINSCutGenerateRegistrations.cxx
  BRAINSCutGenerateProbability.cxx
  BRAINSCutCreateVector.cxx
  FeatureInputVector.cxx
  BRAINSCutApplyModel.cxx
  BRAINSCutTrainModel.cxx
  BRAINSCutVectorTrainingSet.cxx
  BRAINSCutUtilities.cxx
  ${BRAINSCut_HDRS}
  )

#TODO: DELETE CreateVectors.cxx
#TODO:  Rename BRAINSCutDataHandler to something more meaningful

## Build BRAINSCutCOMMONLIB library
##
add_library(BRAINSCutCOMMONLIB STATIC ${BRAINSCutCOMMONLIB_SRCS})
set_target_properties(BRAINSCutCOMMONLIB PROPERTIES FOLDER ${MODULE_FOLDER})

#message(STATUS "${OpenCV_LIBS}")      ## Recent versions of OpenCV changed CMake variable name.
#message(STATUS "${OpenCV_LIBRARIES}")
target_link_libraries(BRAINSCutCOMMONLIB
  NetConfigurationCOMMONLIB
  ${BRAINSCut_ITK_LIBRARIES}
  ${OpenCV_LIBS}
  )

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(BRAINSCutCOMMONLIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(BRAINSCutCOMMONLIB PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

## Build all the programs
##
set(ALL_PROGS_LIST
  BRAINSCut
  ShuffleVectorsModule  ## TODO: Make into SEM program
  SimilarityIndex
  )
foreach( prog ${ALL_PROGS_LIST} )
  StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSCutCOMMONLIB )
endforeach()

add_subdirectory(cnpy) ## https://github.com/rogersce/cnpy Read/write numpy formatted arrays.
add_executable(ImageFeatureExtractor ImageFeatureExtractor.cpp)
target_link_libraries(ImageFeatureExtractor BRAINSCommonLib cnpy)
set_target_properties(ImageFeatureExtractor PROPERTIES FOLDER ${MODULE_FOLDER})

if(1) ## HACK Regina many tests are failing and needs investigation,
      ## HACK Regina some tests are segmentation faulting
  add_subdirectory(TestSuite)
endif()
