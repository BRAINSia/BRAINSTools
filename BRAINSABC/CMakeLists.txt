##- project(BRAINSABC)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKAnisotropicSmoothing
  ITKBinaryMathematicalMorphology
  ITKCommon
  ITKConnectedComponents
  ITKCurvatureFlow
  ITKIOImageBase
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
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Output directories.
#

###
add_subdirectory(brainseg)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()

