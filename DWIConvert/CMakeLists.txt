
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
  ITKNrrdIO
)

#-----------------------------------------------------------------------------
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

set(DWIConvertSupportLib_SRCS
#  FSLToNrrd.cxx
#  NrrdToFSL.cxx
  DWIConvertUtils.cxx
  DWIConverterFactory.h
  DWIConverter.h
  DWIDICOMConverterBase.h
  DWIDICOMConverterBase.cxx
  GEDWIConverter.h
  GenericDWIConverter.h
  HitachiDWIConverter.h
  PhilipsDWIConverter.h
  SiemensDWIConverter.h
  StringContains.h
  )

add_library(DWIConvertSupportLib STATIC ${DWIConvertSupportLib_SRCS})
target_link_libraries(DWIConvertSupportLib ${DWIConvert_ITK_LIBRARIES} BRAINSCommonLib)

set(PROGS_LIST_1
  Convert4DImageTo3DSeries
  )
foreach(prog ${PROGS_LIST_1})
  add_executable(${prog} ${prog}.cxx)
  target_link_libraries(${prog} DWIConvertSupportLib)
endforeach()

#-----------------------------------------------------------------------------
# SEMMacroBuildCLI(
#     NAME ${MODULE_NAME}
#     LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
#     TARGET_LIBRARIES DWIConvertSupportLib
# )

# several files needed down in ExtenededTesting
StandardBRAINSBuildMacro(NAME DWIConvert EXECUTABLE_ONLY
  TARGET_LIBRARIES DWIConvertSupportLib)

set(DWIConvertTest_SRC)
foreach(f ${DWIConvertSupportLib_SRCS})
  list(APPEND DWIConvertTest_SRC
          ${CMAKE_CURRENT_LIST_DIR}/${f})
endforeach()

if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  add_subdirectory(TestSuite)
endif()
