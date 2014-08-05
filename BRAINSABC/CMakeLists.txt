##- project(BRAINSABC)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(ITKAnisotropicSmoothing
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
  ITKTransform)

#-----------------------------------------------------------------------------
# Output directories.
#

###
add_subdirectory(brainseg)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()

