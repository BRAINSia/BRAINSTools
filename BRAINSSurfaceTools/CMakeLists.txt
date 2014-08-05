#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(ITKDeprecated
  ITKIOGDCM
  ITKIOGE
  ITKIOImageBase
  ITKImageGrid
  ITKImageIntensity
  ITKMesh
  ITKOptimizers
  ITKQuadEdgeMesh
  ITKQuadEdgeMeshFiltering
  ITKSpatialObjects
  ITKStatistics
  ITKTransform)

include_directories(${CMAKE_CURRENT_LIST_DIR}/BRAINSSurfaceCommon)

add_subdirectory(BRAINSAssignSurfaceFeatures)
add_subdirectory(BRAINSSurfaceFlattening)
add_subdirectory(BRAINSSurfaceGeneration)
add_subdirectory(BRAINSSurfaceRegister)
#add_subdirectory(BRAINSSurfaceStat)
