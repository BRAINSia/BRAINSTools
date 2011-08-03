project(BRAINSConstellationDetectorTool)
set(LOCAL_PROJECT_NAME BRAINSConstellationDetectorTool)
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
    find_package(ITK 4 REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

#-----------------------------------------------------------------------------
# Output directories.
#
#SETOPTIONALDEBUGIMAGEVIEWER()


add_subdirectory(src)

#-----------------------------------------------------------------------------
## optional build.

# Build test option
if(BUILD_TESTING)
  add_subdirectory(TestSuite)
endif(BUILD_TESTING)

# Build GUI option, the option is provided in src directory
if(BUILD_WITH_GUI OR BRAINS3_USE_QT)
  add_subdirectory(gui)
endif()
