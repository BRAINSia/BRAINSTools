project(BRAINSConstellationDetectorProgs)
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
  )
foreach(prog ${ALL_PROGS_LIST})

  if(0) # Build against Slicer
    ## Include the Slicer macro for setting up default locations!
    SlicerMacroBuildCLI(
      NAME ${prog}
      LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
      TARGET_LIBRARIES landmarksConstellationCOMMONLIB ${OPTIONAL_DEBUG_LINK_LIBRARIES}
      CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
      VERBOSE
    )
  else()
    SEMMacroBuildCLI(
      NAME ${prog}
      LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
      TARGET_LIBRARIES landmarksConstellationCOMMONLIB ${OPTIONAL_DEBUG_LINK_LIBRARIES}
      CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
      VERBOSE
      )
  endif()
endforeach()

