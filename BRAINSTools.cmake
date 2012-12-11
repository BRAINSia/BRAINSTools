
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

if(USE_ANTS)
  # find ANTS includes
  include_directories(${ANTS_SOURCE_DIR}/Temporary)
  include_directories(${ANTS_SOURCE_DIR}/Utilities)
  include_directories(${ANTS_SOURCE_DIR}/Examples)
  include_directories(${ANTS_SOURCE_DIR}/ImageRegistration)
  link_directories(${BRAINSStandAlone_LIBRARY_PATH})
  set(ANTS_LIBS ${ANTS_LIBS} antsUtilities)
endif()

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

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
  BRAINSImageConvert
  BRAINSMush
  BRAINSMultiModeSegment
  BRAINSInitializedControlPoints
  BRAINSTransformConvert
  BRAINSConstellationDetector
  BRAINSABC
  ConvertBetweenFileFormats
  DWIConvert
  BRAINSCreateLabelMapFromProbabilityMaps
  )

if(USE_DebugImageViewer)
  list(APPEND brains_modulenames
    DebugImageViewer)
endif()

#-----------------------------------------------------------------------------
# Add module sub-directory if USE_<MODULENAME> is both defined and true
#-----------------------------------------------------------------------------
foreach(modulename ${brains_modulenames})
  # message("DEFINED USE_${modulename} AND ${USE_${modulename}}")
  if(DEFINED USE_${modulename} AND USE_${modulename})
    add_subdirectory(${modulename})
  endif()
endforeach()

ExternalData_Add_Target( ${PROJECT_NAME}FetchData )  # Name of data management target
