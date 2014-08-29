
##- project(BRAINSConstellationDetector)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(
  ITKThresholding
  ITKImageCompare
  ITKRegistrationCommon
  ITKSmoothing
  ITKConnectedComponents
  ITKMathematicalMorphology
  ITKBinaryMathematicalMorphology
  ITKDistanceMap
  ITKRegionGrowing
  ITKImageCompare
  ITKOptimizersv4
  ITKRegistrationMethodsv4
  ITKVTK
  ITKVtkGlue
  ITKTestKernel
  ITKConvolution
)

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

