
##- project(BRAINSFit)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(BRAINSFit_ITK
  ITKMetricsv4
  ITKOptimizersv4
  ITKSmoothing
  ITKConnectedComponents
  ITKThresholding
  ITKImageCompare
  ITKMathematicalMorphology
  ITKBinaryMathematicalMorphology
  ITKImageLabel
  ITKRegistrationMethodsv4
  ITKDistanceMap
  ITKRegionGrowing
  ITKRegistrationCommon
  ITKTestKernel
)

#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSFitCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSFitCommonLib.h
  )

set(ALL_PROGS_LIST
  BRAINSFit
  PerformMetricTest
  )

set(BRAINSFitLibraries BRAINSCommonLib ${ANTS_LIBS} ${BRAINSFit_ITK_LIBRARIES})
DebugImageViewerLibAdditions(BRAINSFitLibraries)
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES ${BRAINSFitLibraries})
endforeach()

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()
