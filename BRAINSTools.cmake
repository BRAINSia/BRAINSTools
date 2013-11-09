include(${CMAKE_CURRENT_LIST_DIR}/Common.cmake)

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
find_package(ITK REQUIRED)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)

#-----------------------------------------------------------------------------
if(USE_AutoWorkup OR USE_GTRACT OR USE_BRAINSTalairach OR USE_BRAINSSurfaceTools)
  find_package(VTK REQUIRED)
  if(VTK_FOUND)
    include(${VTK_USE_FILE})
  endif()
endif()

if(USE_ANTS)
  # find ANTS includes
  message("ANTs_SOURCE_DIR=${ANTs_SOURCE_DIR}")
  include_directories(${BOOST_INCLUDE_DIR})
  include_directories(${ANTs_SOURCE_DIR}/Temporary)
  include_directories(${ANTs_SOURCE_DIR}/Tensor)
  include_directories(${ANTs_SOURCE_DIR}/Utilities)
  include_directories(${ANTs_SOURCE_DIR}/Examples)
  include_directories(${ANTs_SOURCE_DIR}/ImageRegistration)
  link_directories(${BRAINSTools_LIBRARY_PATH} ${BRAINSTools_CLI_ARCHIVE_OUTPUT_DIRECTORY} ${ANTs_LIBRARY_DIR})
  set(ANTS_LIBS antsUtilities)
endif()

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

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
  BRAINSResample
  BRAINSROIAuto
  GTRACT
  ImageCalculator
  BRAINSCut
  BRAINSLandmarkInitializer
  BRAINSSnapShotWriter
  BRAINSDemonWarp ## NOTE: This is off by default, but is valid for both ITKv3/4
                  ##       This builds just fine with ITKv3/4, but test cases need
                  ##       further review before trusting it.
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
  AutoWorkup
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

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
