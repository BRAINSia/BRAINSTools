
project(DWIConvert)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil( DWIConvert_ITK
  ITKImageIntensity
  ITKImageStatistics
  ITKTestKernel
  ITKIODCMTK
  ITKIORAW
  ITKDCMTK
)

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
target_link_libraries(DWIConvertSupportLib ${DWIConvert_ITK_LIBRARIES} BRAINSCommonLib)

#-----------------------------------------------------------------------------
# SEMMacroBuildCLI(
#     NAME ${MODULE_NAME}
#     LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
#     TARGET_LIBRARIES DWIConvertSupportLib
# )

# several files needed down in ExtenededTesting
StandardBRAINSBuildMacro(NAME DWIConvert EXECUTABLE_ONLY
  TARGET_LIBRARIES DWIConvertSupportLib)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
