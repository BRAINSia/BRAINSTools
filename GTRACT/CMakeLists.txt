
project( GTRACT )

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil( GTRACT_ITK
  ITKOptimizersv4
  ITKDistanceMap
  ITKMathematicalMorphology
  ITKRegionGrowing
  ITKImageCompare
  ITKMetricsv4
  ITKRegistrationMethodsv4
  ITKDiffusionTensorImage
  ITKFastMarching
  ITKImageAdaptors
  ITKImageCompose
  ITKImageFilterBase
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKOptimizers
  ITKRegistrationCommon
  ITKSmoothing
  ITKSpatialObjects
  ITKStatistics
  ITKThresholding
  ITKTransform
  ITKNrrdIO
  )

#-----------------------------------------------------------------------------
# Output directories.
#

option(DEBUGGING_PRINT_IMAGES "Use the DEBUG_IMAGE_VIEWER for debugging" ON)

###
# On Visual Studio 8 MS deprecated C. This removes all 1.276E1265 security
# warnings
# MSVC deprecation warnings suppressed via GTRACTCommon target (Common/CMakeLists.txt)
# VTK linked via target_link_libraries on GTRACTCommon (Common/CMakeLists.txt)
#-----------------------------------------------------------------------------
# Configure files with settings for use by the build.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gtractConfigure.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/gtractConfigure.h)

# GTRACT include dirs are set via target_include_directories on GTRACTCommon (Common/CMakeLists.txt)

set(MODULE_FOLDER "Module-GTRACT")

#-----------------------------------------------------------------------------
# Build Generate CLP
#add_subdirectory( Utils )

#-----------------------------------------------------------------------------
# Build GTRACT Tools
add_subdirectory(Common)
add_subdirectory(Cmdline)
