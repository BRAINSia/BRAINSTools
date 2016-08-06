#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(BRAINSCommonLib_BUILDSCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts)
set(CMAKE_MODULE_PATH
  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}
  ${${PROJECT_NAME}_SOURCE_DIR}/CMake
  ${${PROJECT_NAME}_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )

#-----------------------------------------------------------------------------
# Version information
include(Version.cmake)

set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}")
if(DEFINED ${PROJECT_NAME}_VERSION_PATCH)
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.${${PROJECT_NAME}_VERSION_PATCH}")
  if(DEFINED ${PROJECT_NAME}_VERSION_TWEAK)
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.${${PROJECT_NAME}_VERSION_TWEAK}")
  endif()
endif()

if(DEFINED ${PROJECT_NAME}_VERSION_RC)
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}${${PROJECT_NAME}_VERSION_RC}")
endif()
if(DEFINED ${PROJECT_NAME}_VERSION_POST)
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.post${${PROJECT_NAME}_VERSION_POST}")
elseif(DEFINED ${PROJECT_NAME}_VERSION_DEV)
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.dev${${PROJECT_NAME}_VERSION_DEV}")
endif()

option( ${PROJECT_NAME}_BUILD_DISTRIBUTE "Remove '-g#####' from version. ( for official distribution only )" OFF )
mark_as_advanced( ${PROJECT_NAME}_BUILD_DISTRIBUTE )
if( NOT ${PROJECT_NAME}_BUILD_DISTRIBUTE )
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}-g${${PROJECT_NAME}_VERSION_HASH}")
endif()

message(STATUS "${PROJECT_NAME}_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR}")
message(STATUS "${PROJECT_NAME}_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR}")
message(STATUS "${PROJECT_NAME}_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH}")
message(STATUS "${PROJECT_NAME}_VERSION_TWEAK ${${PROJECT_NAME}_VERSION_TWEAK}")
message(STATUS "${PROJECT_NAME}_VERSION_RC    ${${PROJECT_NAME}_VERSION_RC}")
message(STATUS "${PROJECT_NAME}_VERSION_HASH  ${${PROJECT_NAME}_VERSION_HASH}")
message(STATUS "${PROJECT_NAME}_VERSION_POST  ${${PROJECT_NAME}_VERSION_POST}")
message(STATUS "${PROJECT_NAME}_VERSION_DEV   ${${PROJECT_NAME}_VERSION_DEV}")
message(STATUS "Building ${PROJECT_NAME} version \"${${PROJECT_NAME}_VERSION}\"")

include(FindITKUtil)
include(FindVTKUtil)
# #-----------------------------------------------------------------------------
# if(${PRIMARY_PROJECT_NAME}_REQUIRES_VTK)
# #  message("VTK_DIR:${VTK_DIR}")
#   find_package(VTK REQUIRED)
#   if(VTK_FOUND)
#     include(${VTK_USE_FILE})
#   endif()
# #  message("VTK_USE_FILE:${VTK_USE_FILE}")
# #  message("VTK_INCLUDE_DIRS:${VTK_INCLUDE_DIRS}")
#   include_directories(${VTK_INCLUDE_DIRS})
# endif()

#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)

if(USE_BRAINSABC)
  find_package(TBB REQUIRED)
  # set(VTK_SMP_IMPLEMENTATION_LIBRARIES ${TBB_LIBRARY})
  include_directories(${TBB_INCLUDE_DIRS})
endif()

if(USE_ANTS)
  ## Do a little sanity checking
  if( NOT (EXISTS "${ANTs_SOURCE_DIR}" AND IS_DIRECTORY "${ANTs_SOURCE_DIR}") )
    message(FATAL_ERROR "ANTs_SOURCE_DIR: '${ANTs_SOURCE_DIR}' does not exists")
  endif()
  # find ANTS includes
  include_directories(${BOOST_INCLUDE_DIR})
  include_directories(${ANTs_SOURCE_DIR}/Temporary)
  include_directories(${ANTs_SOURCE_DIR}/Tensor)
  include_directories(${ANTs_SOURCE_DIR}/Utilities)
  include_directories(${ANTs_SOURCE_DIR}/Examples)
  include_directories(${ANTs_SOURCE_DIR}/ImageRegistration)

  if( NOT (EXISTS "${ANTs_LIBRARY_DIR}" AND IS_DIRECTORY "${ANTs_LIBRARY_DIR}") )
    message(FATAL_ERROR "ANTs_LIBRARY_DIR: '${ANTs_LIBRARY_DIR}' does not exists")
  endif()

  link_directories(${BRAINSTools_LIBRARY_PATH} ${BRAINSTools_CLI_ARCHIVE_OUTPUT_DIRECTORY} ${ANTs_LIBRARY_DIR})
  set(ANTS_LIBS antsUtilities)
endif()

# Define the atlas subdirectory in one place
if(USE_ReferenceAtlas)
  set(ReferenceAtlas_XML_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  set(ATLAS_VERSION 20131115)
  set(ATLAS_NAME Atlas/Atlas_${ATLAS_VERSION})
  set(ATLAS_INSTALL_DIRECTORY ${ReferenceAtlas_XML_DIR}/${ATLAS_NAME})
endif()

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

# Note: Projects (e.g. Slicer) integrating BRAINSTools as a subtree that want
#       to disable BRAINSTools testing while managing their own test suite
#       also using the option "BUILD_TESTING" can explicitly set the
#       variable BRAINSTools_DISABLE_TESTING to 1.

#-----------------------------------------------------------------------------
# CTestCustom
#-----------------------------------------------------------------------------
if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  configure_file(
    CMake/CTestCustom.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake
    @ONLY)
endif()

# Some test are failing due to inadequate test construction, but
# the code seems to do the correct thing on real data.
# This is also for tests that are huge, and can not be running
# a long time.
option(ENABLE_EXTENDED_TESTING "Enable tests that are long running, or where the test itself is in error." OFF)
mark_as_advanced(ENABLE_EXTENDED_TESTING)

#Set the global max TIMEOUT for CTest jobs.  This is very large for the moment
#and should be revisted to reduce based on "LONG/SHORT" test times, set to 1 hr for now
set(CTEST_TEST_TIMEOUT 1800 CACHE STRING "Maximum seconds allowed before CTest will kill the test." FORCE)
set(DART_TESTING_TIMEOUT ${CTEST_TEST_TIMEOUT} CACHE STRING "Maximum seconds allowed before CTest will kill the test." FORCE)

## BRAINSTools_MAX_TEST_LEVEL adjusts how agressive the test suite is
## so that long running tests or incomplete tests can easily be
## silenced
## 1 - Run the absolute minimum very fast tests (These should always pass before any code commit)
## 3 - Run fast tests on continous builds (These need immediate attention if they begin to fail)
## 5 - Run moderate nightly tests (These need immediate attention if they begin to fail)
## 7 - Run long running extensive test that are a burden to normal development (perhaps test 1x per week)
## 8 - Run tests that fail due to incomplete test building, these are good ideas for test that we don't have time to make robust)
## 9 - Run silly tests that don't have much untility
set(BRAINSTools_MAX_TEST_LEVEL 4 CACHE STRING "Testing level for managing test burden")

#-----------------------------------------------------------------------
# Setup locations to find externally maintained test data.
#-----------------------------------------------------------------------
include(BRAINSToolsExternalData)

set(TestData_DIR ${CMAKE_CURRENT_SOURCE_DIR}/TestData)

#
# choose between using HDF5 or MAT format transform files
set(XFRM_EXT "h5" CACHE STRING "Choose the preferred transform file format")

#-----------------------------------------------------------------------------
# BRAINSCommonLib (Required)
#-----------------------------------------------------------------------------
include(CMakeBRAINS3BuildMacros)

add_subdirectory(BRAINSCommonLib)

set(BRAINSCommonLib_DIR ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib)
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib)

#-----------------------------------------------------------------------------
# Define list of module names
#-----------------------------------------------------------------------------
set(brains_modulenames
  BRAINSFit
  BRAINSLabelStats
  BRAINSResample
  BRAINSROIAuto
  GTRACT
  ImageCalculator
  BRAINSCut
  ## Temporarily Removed Need to update OpenCV BRAINSCut
  BRAINSLandmarkInitializer
  BRAINSSnapShotWriter
  BRAINSDemonWarp ## NOTE: This is off by default, but is valid for both ITKv3/4
                  ##       This builds just fine with ITKv3/4, but test cases need
                  ##       further review before trusting it.
  ##TODO: KENT:  This is broken with latest builds,  I think something in ITKv4 changed slightly, or VTK6 compatibility -->BRAINSSurfaceTools
  BRAINSSurfaceTools
  ICCDEF
  BRAINSContinuousClass
  BRAINSPosteriorToContinuousClass
  BRAINSMush
  BRAINSMultiModeSegment
  BRAINSInitializedControlPoints
  BRAINSTransformConvert
  BRAINSTalairach
  BRAINSConstellationDetector
  BRAINSABC
  ConvertBetweenFileFormats
  DWIConvert
  BRAINSCreateLabelMapFromProbabilityMaps
  BRAINSMultiSTAPLE
  BRAINSStripRotation
  AutoWorkup
  BRAINSDWICleanup
  ReferenceAtlas
  BRAINSSuperResolution
  )

if(USE_DebugImageViewer)
  list(APPEND brains_modulenames
    DebugImageViewer)
endif()

## HACK: This is needed to get DWIConvert to build in installed tree
## KENT: Please remove this line and make DWIConvert build by fixing ITK install of DCMTK
include_directories(${ITK_INSTALL_PREFIX}/install)

#-----------------------------------------------------------------------------
# Add module sub-directory if USE_<MODULENAME> is both defined and true
#-----------------------------------------------------------------------------
set(BRAINSToolsModules "")
foreach(modulename ${brains_modulenames})
  # message("DEFINED USE_${modulename} AND ${USE_${modulename}}")
  if(DEFINED USE_${modulename} AND USE_${modulename})
  #  message("Adding ${modulename}")
    add_subdirectory(${modulename})
    list(APPEND BRAINSToolsModules ${modulename})
  #else()
  #  message("USE_${modulename} = ${USE_${modulename}}")
  endif()
endforeach()

if(USE_ITKMatlabIO)
  add_subdirectory(ITKMatlabIO)
endif()

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
