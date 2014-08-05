
##- project(BRAINSConstellationDetector)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(ITKHDF5
  ITKRegistrationMethodsv4
  ITKDistanceMap
  ITKImageFilterBase
  ITKImageFunction
  ITKImageFusion
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKRegistrationCommon
  ITKSmoothing
  ITKSpatialObjects
  ITKStatistics
  ITKThresholding
  ITKImageCompare
  ITKTransform
  ITKOptimizersv4
  ITKMetricsv4
  ITKRegionGrowing
  ITKTestKernel
  ITKVtkGlue)

#-----------------------------------------------------------------------------
# Output directories.
#
add_subdirectory(src)
add_subdirectory(landmarkStatistics)

#-----------------------------------------------------------------------------
## optional build.

# Build GUI option, the option is provided in src directory
if(BUILD_WITH_GUI OR BRAINS3_USE_QT OR BRAINSTools_USE_QT)
  add_subdirectory(gui)
endif()

# Build test option
if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()

