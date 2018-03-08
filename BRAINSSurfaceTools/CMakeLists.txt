#-----------------------------------------------------------------------------
# Dependencies.
#

FindITKUtil(BRAINSSurfaceTools_ITK
  ITKAnisotropicSmoothing
  ITKBinaryMathematicalMorphology
  ITKCommon
  ITKConnectedComponents
  ITKDisplacementField
  ITKDistanceMap
  ITKIOGDCM
  ITKIOGE
  ITKIOImageBase
  ITKIOTransformBase
  ITKImageCompose
  ITKImageFilterBase
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKLabelMap
  ITKMathematicalMorphology
  ITKMesh
  ITKOptimizers
  ITKPath
  ITKQuadEdgeMesh
  ITKQuadEdgeMeshFiltering
  ITKRegionGrowing
  ITKSmoothing
  ITKSpatialObjects
  ITKStatistics
  ITKThresholding
  ITKTransform
  ITKVTK
  ITKVtkGlue
)
FindVTKUtil(
  vtkCommonComputationalGeometry
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMath
  vtkCommonMisc
  vtkCommonSystem
  vtkCommonTransforms
  vtkFiltersCore
  vtkFiltersExtraction
  vtkFiltersGeneral
  vtkFiltersGeometry
  vtkFiltersStatistics
  vtkIOCore
  vtkIOGeometry
  vtkIOImage
  vtkIOLegacy
  vtkIOXML
  vtkIOXMLParser
  vtkImagingCore
  vtkImagingFourier
  vtkalglib
  vtkRenderingCore

  vtkInteractionStyle
  vtkImagingSources
  )

set(vtkITK_FOLDER "Module-BRAINSSurfaceTools")
add_subdirectory(vtkITK)
include_directories(${CMAKE_CURRENT_LIST_DIR}/vtkITK ${CMAKE_CURRENT_BINARY_DIR}/vtkITK)
include_directories(${CMAKE_CURRENT_LIST_DIR}        ${CMAKE_CURRENT_BINARY_DIR})

#
# ITK
#
include_directories(${CMAKE_CURRENT_LIST_DIR}/BRAINSSurfaceCommon)

add_subdirectory(BRAINSAssignSurfaceFeatures)
add_subdirectory(BRAINSSurfaceFlattening)
add_subdirectory(BRAINSSurfaceGeneration)
add_subdirectory(BRAINSSurfaceRegister)
#add_subdirectory(BRAINSSurfaceStat)
add_subdirectory(GenusZeroImageFilter)
