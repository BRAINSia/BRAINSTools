#-----------------------------------------------------------------------------
# Dependencies.
#

FindITKUtil(
  ITKQuadEdgeMesh
  ITKOptimizers
  ITKQuadEdgeMeshFiltering
)
FindVTKUtil(vtkIOLegacy
  vtkIOXML
  vtkCommonDataModel
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersGeometry
  vtkImagingCore
  vtkFiltersExtraction
  )

#
# ITK
#
include_directories(${CMAKE_CURRENT_LIST_DIR}/BRAINSSurfaceCommon)

add_subdirectory(BRAINSAssignSurfaceFeatures)
add_subdirectory(BRAINSSurfaceFlattening)
add_subdirectory(BRAINSSurfaceGeneration)
add_subdirectory(BRAINSSurfaceRegister)
#add_subdirectory(BRAINSSurfaceStat)
