
##- project(BRAINSConstellationDetector)

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
  ITKHDF5
  ITKIOImageBase
  ITKIOTransformBase
  ITKImageFilterBase
  ITKImageFunction
  ITKImageFusion
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKRegistrationCommon
  ITKSmoothing
  ITKSpatialObjects
  ITKStatistics
  ITKThresholding
  ITKTransform
  ITKVTK
  ITKVtkGlue
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

