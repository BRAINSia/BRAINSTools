##- project(BRAINSABC)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(ITKAnisotropicSmoothing
  ITKDistanceMap
  ITKOptimizersv4
  ITKMetricsv4
  ITKRegistrationMethodsv4
  ITKRegionGrowing
  ITKBinaryMathematicalMorphology
  ITKConnectedComponents
  ITKCurvatureFlow
  ITKImageFilterBase
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKMathematicalMorphology
  ITKRegistrationCommon
  ITKSmoothing
  ITKStatistics
  ITKTestKernel
  ITKThresholding
  ITKTransform
  ITKImageCompare)

#-----------------------------------------------------------------------------
# Output directories.
#

###
add_subdirectory(brainseg)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()

