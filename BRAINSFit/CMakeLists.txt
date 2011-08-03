project(BRAINSFit)
set(LOCAL_PROJECT_NAME BRAINSFit)
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


configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSFitCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSFitCommonLib.h
  )

## ITKv4 requires that ITKReview is explicitly linked against
if(ITK_VERSION_MAJOR EQUAL 4 )
  set(ITKREVIEWLIB ITKReview)
endif(ITK_VERSION_MAJOR EQUAL 4 )

###
if(0) # Build against Slicer
  ## Include the Slicer macro for setting up default locations!
  SlicerMacroBuildCLI(
    NAME BRAINSFit
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
  )
else()
  SEMMacroBuildCLI(
    NAME BRAINSFit
    LOGO_HEADER ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h
    TARGET_LIBRARIES BRAINSCommonLib ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES}
    CLI_SHARED_LIBRARY_WRAPPER_CXX ${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx
    VERBOSE
    )
endif()

if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)
