project(BRAINSConstellationDetectorTool)
set(LOCAL_PROJECT_NAME BRAINSConstellationDetectorTool)
cmake_minimum_required(VERSION 2.8)


## Build landmarksConstellationCOMMONLIB library
##
add_library(landmarksConstellationCOMMONLIB STATIC
  landmarksConstellationCommon.cxx landmarkIO.cxx
  landmarksConstellationDetector.cxx
  TrimForegroundInDirection.cxx
  LLSModel.cxx)
target_link_libraries(landmarksConstellationCOMMONLIB BRAINSCommonLib ${ITK_LIBRARIES})

## Build all the programs
##
set(ALL_PROGS_LIST
  BRAINSConstellationModeler
  BRAINSLinearModelerEPCA
  BRAINSConstellationDetector
  BRAINSAlignMSP
  BRAINSClipInferior
  BRAINSTrimForegroundInDirection
  BRAINSLmkTransform
  TransformFromFiducials
  BRAINSEyeDetector
  insertMidACPCpoint
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES landmarksConstellationCOMMONLIB )
endforeach()

