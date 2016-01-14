##- project(BRAINSABC)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(BRAINSABC_ITK
  ITKAnisotropicSmoothing
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

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()

