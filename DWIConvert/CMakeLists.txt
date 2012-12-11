project (DWIConvert)

#-----------------------------------------------------------------------------
enable_testing()

find_package(ITK 4 REQUIRED)

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
target_link_libraries(DWIConvertSupportLib ${ITK_LIBRARIES}  )

set(MODULE_NAME DWIConvert)
#-----------------------------------------------------------------------------
# SEMMacroBuildCLI(
#     NAME ${MODULE_NAME}
#     LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
#     TARGET_LIBRARIES DWIConvertSupportLib
# )

# several files needed down in ExtenededTesting
StandardBRAINSBuildMacro(NAME DWIConvert
  TARGET_LIBRARIES DWIConvertSupportLib)

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  add_subdirectory(TestSuite)
endif()
