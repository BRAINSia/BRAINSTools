project(BRAINSInitializedControlPoints)
set(LOCAL_PROJECT_NAME BRAINSInitializedControlPoints)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

enable_testing()
include(Dart)
include(CPack)
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBuildMacros.cmake)
include(${BRAINS_CMAKE_HELPER_DIR}/PreventInSourceBuilds.cmake)

###
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

###
CHECKIFSLICER3BUILD()  ## Call the convenience macro

if(NOT ITK_FOUND)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

#-----------------------------------------------------------------------------
# Output directories.
#
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBRAINS3BuildMacros.cmake)

# No need to add with nothing to compile. add_subdirectory(BRAINSInitializedControlPoints_Common)
include(${BRAINS_CMAKE_HELPER_DIR}/IJMacros.txt)
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBRAINS3BuildMacros.cmake)
###


set(CLP BRAINSInitializedControlPoints)
set(${CLP}_SOURCE BRAINSInitializedControlPoints.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSInitializedControlPoints.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries(BRAINSInitializedControlPoints ${ITK_LIBRARIES} crypto)

