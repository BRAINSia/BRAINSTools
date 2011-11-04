project(BRAINSCut)

#message(STATUS "BRAINSCut OpenCV_DIR = ${OpenCV_DIR}")
find_package( OpenCV REQUIRED)

###
add_subdirectory(BRAINSFeatureCreators)

include_directories(${BRAINSCut_SOURCE_DIR}/NetConfiguration)
add_subdirectory(NetConfiguration)


file(GLOB BRAINSCut_HDRS *.h *.txx)

set(BRAINSCutCOMMONLIB_SRCS
  Utilities.cxx
  GenerateRegistrations.cxx
  XMLParser.cxx
  NetConfigurationParser.cxx
  ShuffleVectors.cxx
  BRAINSCutGenerateProbability.cxx
  BRAINSCutCreateVector.cxx
  FeatureInputVector.cxx
  BRAINSCutApplyModel.cxx
  BRAINSCutTrainModel.cxx
  BRAINSCutVectorTrainingSet.cxx
  GenerateProbability.cxx
  TrainModel.cxx
  ApplyModel.cxx
  BRAINSCutPrimary.cxx
  ${BRAINSCut_HDRS}
  )

#TODO: DELETE CreateVectors.cxx
#TODO:  Rename BRAINSCutPrimary to something more meaningful

## Build BRAINSCutCOMMONLIB library
##
add_library(BRAINSCutCOMMONLIB STATIC ${BRAINSCutCOMMONLIB_SRCS})
target_link_libraries(BRAINSCutCOMMONLIB
  NetConfigurationCOMMONLIB
  ${ITK_LIBRARIES}
  ${OpenCV_LIBS}
  )


## Build all the programs
##
set(ALL_PROGS_LIST
  BRAINSCut
  ShuffleVectorsModule  ## TODO: Make into SEM program
  WriteFeaturesforROI   ## TODO: Make into SEM program
  SimilarityIndex 
  )
foreach( prog ${ALL_PROGS_LIST} )
  StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSCutCOMMONLIB )
endforeach()

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
