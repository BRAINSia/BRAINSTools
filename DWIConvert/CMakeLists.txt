
project(DWIConvert)

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
  ITKIOPNG
  ITKIOStimulate
  ITKIOTIFF
  ITKIOTransformInsightLegacy
  ITKIOVTK
  )
find_package(ITK COMPONENTS
  ITKCommon
  ITKIODCMTK
  ITKIOImageBase
  ITKIONRRD
  ITKIORAW
  ITKImageIntensity
  ITKImageStatistics
  ITKTestKernel
  ${ITK_IO_MODULES_USED}
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

set(DWIConvertSupportLib_SRCS
  FSLToNrrd.cxx
  NrrdToFSL.cxx
  )

set(DWIConvertTest_SRC)
foreach(f ${DWIConvertSupportLib_SRCS})
  list(APPEND DWIConvertTest_SRC
    ${CMAKE_CURRENT_LIST_DIR}/${f})
endforeach()

add_library(DWIConvertSupportLib STATIC ${DWIConvertSupportLib_SRCS})
target_link_libraries(DWIConvertSupportLib ${ITK_LIBRARIES} )

#-----------------------------------------------------------------------------
# SEMMacroBuildCLI(
#     NAME ${MODULE_NAME}
#     LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
#     TARGET_LIBRARIES DWIConvertSupportLib
# )

# several files needed down in ExtenededTesting
StandardBRAINSBuildMacro(NAME DWIConvert EXECUTABLE_ONLY
  TARGET_LIBRARIES BRAINSCommonLib DWIConvertSupportLib)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
