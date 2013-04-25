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
find_package(VTK REQUIRED)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
endif()

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

if(USE_ANTS)
  # find ANTS includes
  message( STATUS "XXXXXXXXXXXXX ${BOOST_INCLUDE_DIR} XXXXXXX")
  include_directories(${BOOST_INCLUDE_DIR})
  include_directories(${ANTS_SOURCE_DIR}/Temporary)
  include_directories(${ANTS_SOURCE_DIR}/Tensor)
  include_directories(${ANTS_SOURCE_DIR}/Utilities)
  include_directories(${ANTS_SOURCE_DIR}/Examples)
  include_directories(${ANTS_SOURCE_DIR}/ImageRegistration)
  link_directories(${BRAINSTools_LIBRARY_PATH})
  set(ANTS_LIBS ${ANTS_LIBS} antsUtilities)
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
  AutoWorkup
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
foreach(modulename ${brains_modulenames})
  # message("DEFINED USE_${modulename} AND ${USE_${modulename}}")
  if(DEFINED USE_${modulename} AND USE_${modulename})
  #  message("Adding ${modulename}")
    add_subdirectory(${modulename})
  #else()
  #  message("USE_${modulename} = ${USE_${modulename}}")
  endif()
endforeach()

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
