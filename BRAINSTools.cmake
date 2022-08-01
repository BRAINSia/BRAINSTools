if(BUILD_COVERAGE)
  if(BUILD_OPTIMIZED)
    message(FATAL_ERROR "Can not build optimized when building for coverage, debug information needed")
  endif()
  if(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    message(FATAL_ERROR "BUILD_COVERAGE Requires Debug build")
  endif()
  message(INFO " ${CMAKE_CXX_COMPILER_ID}")
  if( CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options( -g -O0 --coverage )
    add_link_options(--coverage)
  else()
    message(FATAL_ERROR "COVERAGE Requires GNU or Clang compilers.")
  endif()
    message(INFO " CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
    message(INFO " CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
    message(INFO " CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
    message(INFO " CMAKE_MODULE_LINKER_FLAGS: ${CMAKE_MODULE_LINKER_FLAGS}")
endif()


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
#-----------------------------------------------------------------------------
if(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)
  include(FindVTKUtil)
  #  message("VTK_DIR:${VTK_DIR}")
  find_package(VTK 9.1 REQUIRED NO_MODULE)
  #  message("VTK_INCLUDE_DIRS:${VTK_INCLUDE_DIRS}")
  include_directories(${VTK_INCLUDE_DIRS})
endif()

#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)

if(USE_BRAINSABC)
  if( NOT USE_AutoWorkup )
     message(FATAL_ERROR "BRAINSABC requires USE_AutoWorkup to be ON: ${USE_BRAINSABC} != ${USE_AutoWorkup}")
  endif()
  if( APPLE )
    set( TBB_MIN_VERSION "2019.1") ## Actually 2019.0.11002 is needed for when OSX MIN version < 10.12
  else()
    set( TBB_MIN_VERSION "2017.0")
  endif()
  find_package(TBB ${TBB_MIN_VERSION} REQUIRED
#               COMPONENTS tbb tbbmalloc
               NO_MODULE PATHS ${TBB_DIR} )
#message(FATAL_ERROR "${TBB_DIR}")

  # set(VTK_SMP_IMPLEMENTATION_LIBRARIES ${tbb_LIBRARY})
  include_directories(${tbb_INCLUDE_DIRS})
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

# Note: Projects (e.g. Slicer) integrating BRAINSTools as a subtree that want
#       to disable BRAINSTools testing while managing their own test suite
#       also using the option "${LOCAL_PROJECT_NAME}_BUILD_TESTING" can explicitly set the
#       variable BRAINSTools_DISABLE_TESTING to 1.

#-----------------------------------------------------------------------------
# CTestCustom
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
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
# long tests should explicitly set the "TIMEOUT" property
# https://stackoverflow.com/questions/45009595/how-to-overwrite-ctest-default-timeout-1500-in-cmakelists-txt
#HACK THIS DOES NOT CHANGE TEST TIMEOUTS set(CTEST_TEST_TIMEOUT 60 CACHE STRING "Maximum seconds allowed before CTest will kill the test." FORCE)

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(BRAINSCommonLib_BUILDSCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts)
list(INSERT CMAKE_MODULE_PATH 0
  ${BRAINSCommonLib_BUILDSCRIPTS_DIR}
  ${${PROJECT_NAME}_SOURCE_DIR}/CMake
  ${${PROJECT_NAME}_BINARY_DIR}/CMake
  )

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
  ReferenceAtlas
  BRAINSFit
  BRAINSResample
  BRAINSABC
  BRAINSROIAuto
  DWIConvert
  ConvertBetweenFileFormats
  BRAINSConstellationDetector
  ImageCalculator
  GTRACT
  BRAINSLabelStats
  BRAINSLandmarkInitializer
  BRAINSSnapShotWriter
  BRAINSDWICleanup
  BRAINSPosteriorToContinuousClass
  BRAINSMush
  BRAINSMultiModeSegment
  BRAINSInitializedControlPoints
  BRAINSTransformConvert
  BRAINSCreateLabelMapFromProbabilityMaps
  BRAINSMultiSTAPLE
  BRAINSStripRotation
  BRAINSSuperResolution
  BRAINSDeface
  BRAINSIntensityNormalize
)

if(BUILD_ARCHIVE)
   list(APPEND brains_modulenames
  BRAINSTalairach
)

endif()

if(USE_DebugImageViewer)
  list(APPEND brains_modulenames DebugImageViewer)
endif()

## HACK: This is needed to get DWIConvert to build in installed tree
## KENT: Please remove this line and make DWIConvert build by fixing ITK install of DCMTK
# --include_directories(${ITK_INSTALL_PREFIX}/install)

#-----------------------------------------------------------------------------
# Add module sub-directory if USE_<MODULENAME> is both defined and true
#-----------------------------------------------------------------------------
set(BRAINSToolsModules "")
foreach(modulename ${brains_modulenames})
  # message("DEFINED USE_${modulename} AND ${USE_${modulename}}")
  if(DEFINED USE_${modulename} AND USE_${modulename})
    if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${modulename})
      #message("++++++++++++++++++++++++++++++++++++++Adding ${modulename}")
      add_subdirectory(${modulename})
      #message("--------------------------------------USE_${modulename} = ${USE_${modulename}}")
    else()
      if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ARCHIVE/${modulename})
        add_subdirectory( ARCHIVE/${modulename} )
      else()
        message(FATAL_ERROR "ERROR: Missing directory for module '${modulename}' and 'ARCHIVE/${modulename}'")
      endif()
    endif()
    list(APPEND BRAINSToolsModules ${modulename})
  else()
  endif()
endforeach()

if(USE_ITKMatlabIO)
  add_subdirectory(ITKMatlabIO)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cpack_brainstools.cmake)
# Name of data management target
if(${CMAKE_PROJECT_NAME} STREQUAL "BRAINSTools")
  ExternalData_Add_Target( ${BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET} SHOW_PROGRESS ON )
endif()
