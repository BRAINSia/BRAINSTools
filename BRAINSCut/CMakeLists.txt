project(BRAINSCut)

message(STATUS "BRAINSCut OpenCV_DIR = ${OpenCV_DIR}")
find_package( OpenCV REQUIRED)

###
add_subdirectory(BRAINSFeatureCreators)

include_directories(${BRAINSCut_SOURCE_DIR}/BRAINSCutConfiguration)
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
#message(FATAL_ERROR "${OpenCV_LIBS}")
target_link_libraries(BRAINSCutCOMMONLIB
  NetConfigurationCOMMONLIB
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
  ${VTK_LIBRARIES}
  )


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

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
