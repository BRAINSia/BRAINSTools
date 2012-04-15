project(BRAINSCommonLibProject)
#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib.h
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

set(BRAINSCommonLib_SRCS GenericTransformImage.cxx BRAINSFitHelper.cxx Slicer3LandmarkWeightIO.cxx Slicer3LandmarkIO.cxx itkOrthogonalize3DRotationMatrix.cxx BRAINSThreadControl.cxx)

## Always build BRAINSCommonLib as static
add_library(BRAINSCommonLib STATIC ${BRAINSCommonLib_SRCS})
target_link_libraries(BRAINSCommonLib ${ITK_LIBRARIES} )

install(TARGETS BRAINSCommonLib
  RUNTIME DESTINATION bin                 COMPONENT Development
  LIBRARY DESTINATION lib/BRAINSCommonLib COMPONENT Development
  ARCHIVE DESTINATION lib/BRAINSCommonLib COMPONENT Development
    )

file(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
install(FILES ${__files1} ${__files2} ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib.h DESTINATION include/BRAINSCommonLib COMPONENT Development)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLibConfig.cmake.in
               ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLibConfig.cmake
               @ONLY IMMEDIATE)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLibConfig.cmake DESTINATION lib/BRAINSCommonLib COMPONENT Development RENAME BRAINSCommonLibConfig.cmake)

set(BRAINSCommonLib_BUILDSCRIPTS_DIR_CONFIG ${CMAKE_CURRENT_BINARY_DIR}/BuildScripts)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/UseBRAINSCommonLib.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/UseBRAINSCommonLib.cmake @ONLY IMMEDIATE)

# install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/UseBRAINSCommonLib.cmake.in
#         DESTINATION lib/BRAINSCommonLib
#         COMPONENT Development
#         RENAME UseBRAINSCommonLib.cmake)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UseBRAINSCommonLib.cmake
  DESTINATION lib/BRAINSCommonLib
  COMPONENT Development)

#=======================================================================
# Install the common build scripts for BRAINS tools in a common location
#=======================================================================
file(GLOB BRAINSBUILD_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "BuildScripts/*")
foreach(BUILDFILE ${BRAINSBUILD_FILES})
  #message(STATUS "${BUILDFILE} ${CMAKE_CURRENT_BINARY_DIR}/${BUILDFILE}")
  configure_file(${BUILDFILE} ${CMAKE_CURRENT_BINARY_DIR}/${BUILDFILE} @ONLY IMMEDIATE)
endforeach()
install(FILES ${BRAINSBUILD_FILES} DESTINATION include/BRAINSCommonLib/BuildScripts COMPONENT Development)

#=======================================================================
# Install the common build scripts for BRAINS tools in a common location
#=======================================================================
file(GLOB_RECURSE BRAINSTEST_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/TestData" "*.md5")
set(TESTFILES_TO_INSTALL "")
foreach(TESTFILE ${BRAINSTEST_FILES})
  #message(STATUS "### ${TESTFILE} ${CMAKE_CURRENT_BINARY_DIR}/TestData/${TESTFILE}")
  configure_file("TestData/${TESTFILE}" ${CMAKE_CURRENT_BINARY_DIR}/TestData/${TESTFILE} @ONLY IMMEDIATE)
  list(APPEND TESTFILES_TO_INSTALL "TestData/${TESTFILE}")
endforeach()
install(FILES ${TESTFILES_TO_INSTALL} DESTINATION include/BRAINSCommonLib/TestData COMPONENT Development)

#CONFIGUREBRAINSORSLICERLIBRARY( BRAINSCommonLib "" ${BRAINSCommonLib_SRCS} "")

##HACK NEED BETTER TESTS add_directory( TestLargestForegroundFilledMaskImageFilter )
##HACK NEED BETTER TESTS add_directory( Test_FindCenterOfBrainFilter )
if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif(BUILD_TESTING)

