project(BRAINSConstellationDetectorTool)

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
  BRAINSTransformFromFiducials
  BRAINSEyeDetector
  insertMidACPCpoint
  fcsv_to_hdf5
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES landmarksConstellationCOMMONLIB )
endforeach()
