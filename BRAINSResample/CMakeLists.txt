
##- project(BRAINSResample)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
set(ITK_IO_MODULES_USED
  ITKIOBMP
  ITKIOBioRad
  ITKIOGDCM
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
  ITKDisplacementField
  ITKDistanceMap
  ITKIOImageBase
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKImageSources
  ITKImageStatistics
  ITKThresholding
  ITKTransform
  ${ITK_IO_MODULES_USED}
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSResampleCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSResampleCommonLib.h
  )

set(ALL_PROGS_LIST
  BRAINSResample
  BRAINSResize
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

### Add the Testing Subdirectory.
if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
