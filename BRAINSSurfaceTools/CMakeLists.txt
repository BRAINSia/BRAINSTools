#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
set(ITK_IO_MODULES_USED
  ITKIOBMP
  ITKIOBioRad
  ITKIOGIPL
  ITKIOJPEG
  ITKIOLSM
  ITKIOMeta
  ITKIONIFTI
  ITKIONRRD
  ITKIOPNG
  ITKIOStimulate
  ITKIOTIFF
  ITKIOTransformInsightLegacy
  ITKIOVTK
  )
find_package(ITK COMPONENTS
  ITKCommon
  ITKDeprecated
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
  ITKTransform
  ${ITK_IO_MODULES_USED}
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

include_directories(${CMAKE_CURRENT_LIST_DIR}/BRAINSSurfaceCommon)

add_subdirectory(BRAINSAssignSurfaceFeatures)
add_subdirectory(BRAINSSurfaceFlattening)
add_subdirectory(BRAINSSurfaceGeneration)
add_subdirectory(BRAINSSurfaceRegister)
#add_subdirectory(BRAINSSurfaceStat)
