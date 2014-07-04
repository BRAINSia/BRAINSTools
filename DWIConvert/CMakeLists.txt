
project(DWIConvert)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKCommon
  ITKIODCMTK
  ITKIOImageBase
  ITKIONRRD
  ITKIORAW
  ITKImageIntensity
  ITKImageStatistics
  ITKTestKernel
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

#
# TEMPORARY FIX -- added itkDCMTKFileReader.cxx locally until
# the next ITK blockpoint release.
if("${ITK_VERSION_MAJOR}" VERSION_GREATER "3"
    AND "${ITK_VERSION_MINOR}" VERSION_GREATER "5"
    AND "${ITK_VERSION_PATCH}" VERSION_GREATER "0")
message("ITK Version 4.6.1 and greater no longer require fix for itk::DCMTKFileReader
itkDCMTKFileReader.cxx & itkDCMTKFileReader.h can be removed from the repository,
and itkDCMTKFileReader.cxx can be removed from DWIConvertSupportLib_SRCS variable
in DWIConvert/CMakeLists.txt")
endif()

set(DWIConvertSupportLib_SRCS
  FSLToNrrd.cxx
  NrrdToFSL.cxx
  itkDCMTKFileReader.cxx
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
