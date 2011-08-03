project(BRAINSABC)
set(LOCAL_PROJECT_NAME BRAINSABC)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

enable_testing()
include(CTest)

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/PreventInSourceBuilds.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/CMakeBuildMacros.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMMacroBuildCLI.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/CMakeBRAINS3BuildMacros.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/IJMacros.txt)

###
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

if(NOT ITK_FOUND)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

#-----------------------------------------------------------------------------
# Output directories.
#
#SETOPTIONALDEBUGIMAGEVIEWER()

###
add_subdirectory(brainseg)

if(1)
  add_subdirectory(TestSuite)
endif(1)

