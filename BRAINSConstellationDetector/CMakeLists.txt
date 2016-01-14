
##- project(BRAINSConstellationDetector)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(BRAINSConstellationDetector_ITK
  ITKBinaryMathematicalMorphology
  ITKCommon
  ITKConnectedComponents
  ITKConvolution
  ITKDistanceMap
  ITKHDF5
  ITKIOImageBase
  ITKIOTransformBase
  ITKImageCompare
  ITKImageFilterBase
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKLabelMap
  ITKMathematicalMorphology
  ITKOptimizersv4
  ITKRegionGrowing
  ITKRegistrationCommon
  ITKRegistrationMethodsv4
  ITKSmoothing
  ITKSpatialObjects
  ITKStatistics
  ITKTestKernel
  ITKThresholding
  ITKTransform
  ITKVTK
  ITKVtkGlue
)

# Build GUI option, the option is provided in src directory
if(BUILD_WITH_GUI OR BRAINS3_USE_QT OR BRAINSTools_USE_QT)
  set(LOCAL_VTKQT vtkGUISupportQt)
endif()

FindVTKUtil(
  vtkCommonColor
  vtkCommonComputationalGeometry
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMath
  vtkCommonMisc
  vtkCommonSystem
  vtkCommonTransforms
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersSources
  vtkImagingCore
  vtkInteractionStyle
  vtkRenderingCore
${LOCAL_VTKQT})

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
if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()

