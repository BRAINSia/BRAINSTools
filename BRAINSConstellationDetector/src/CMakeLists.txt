## Build landmarksConstellationCOMMONLIB library
##
add_library(landmarksConstellationCOMMONLIB STATIC
  landmarksConstellationCommon.cxx landmarkIO.cxx
  landmarksConstellationDetector.cxx
  TrimForegroundInDirection.cxx
  LLSModel.cxx
  PrepareOutputImages.cxx
  BRAINSConstellationDetectorPrimary.cxx)

target_link_libraries(landmarksConstellationCOMMONLIB BRAINSCommonLib ${BRAINSConstellationDetector_ITK_LIBRARIES} ${VTK_LIBRARIES})


#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(landmarksConstellationCOMMONLIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(landmarksConstellationCOMMONLIB PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

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
  landmarksConstellationAligner
  landmarksConstellationWeights
  BinaryMaskEditorBasedOnLandmarks
  BRAINSConstellationLandmarksTransform
  # ComputeReflectiveCorrelationMetric # --A debugging program, should be compiled when needed.
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES landmarksConstellationCOMMONLIB ${VTK_LIBRARIES} )

endforeach()
