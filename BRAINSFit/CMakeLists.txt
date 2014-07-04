
##- project(BRAINSFit)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  ITKIOTransformBase
  ITKImageFilterBase
  ITKImageFunction
  ITKImageGrid
  ITKMetricsv4
  ITKOptimizers
  ITKSmoothing
  ITKSpatialObjects
  ITKTransform
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSFitCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSFitCommonLib.h
  )

set(ALL_PROGS_LIST
  BRAINSFit
  )

set(BRAINSFitLibraries BRAINSCommonLib ${ANTS_LIBS})
DebugImageViewerLibAdditions(BRAINSFitLibraries)
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES ${BRAINSFitLibraries})
endforeach()

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
    add_subdirectory(TestSuite)
endif()
