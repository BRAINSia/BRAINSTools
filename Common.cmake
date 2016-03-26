#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)


include(ExternalProjectDependency)

include(CMakeDependentOption)

option(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT "Install development support include and libraries for external packages." OFF)
mark_as_advanced(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)

option(${LOCAL_PROJECT_NAME}_USE_QT "Find and use Qt with VTK to build GUI Tools" OFF)
mark_as_advanced(${LOCAL_PROJECT_NAME}_USE_QT)

set(USE_ITKv4 ON)
set(ITK_VERSION_MAJOR 4 CACHE STRING "Choose the expected ITK major version to build BRAINS only version 4 allowed.")
# Set the possible values of ITK major version for cmake-gui
set_property(CACHE ITK_VERSION_MAJOR PROPERTY STRINGS "4")
set(expected_ITK_VERSION_MAJOR ${ITK_VERSION_MAJOR})
if(${ITK_VERSION_MAJOR} VERSION_LESS ${expected_ITK_VERSION_MAJOR})
  # Note: Since ITKv3 doesn't include a ITKConfigVersion.cmake file, let's check the version
  #       explicitly instead of passing the version as an argument to find_package() command.
  message(FATAL_ERROR "Could not find a configuration file for package \"ITK\" that is compatible "
                      "with requested version \"${expected_ITK_VERSION_MAJOR}\".\n"
                      "The following configuration files were considered but not accepted:\n"
                      "  ${ITK_CONFIG}, version: ${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}\n")
endif()


#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set(CTEST_CONFIGURATION_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of test." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
  #  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(${ITK_VERSION_MAJOR} STREQUAL "3")
  message(FATAL_ERROR "ITKv3 is no longer supported")
endif()

#-----------------------------------------------------------------------------
# Build option(s)
#-----------------------------------------------------------------------------
option(USE_AutoWorkup                     "Build AutoWorkup"                     ON)
option(USE_ANTS                           "Build ANTS"                           ON)
option(USE_ReferenceAtlas                 "Build the Reference Atlas"            ON)

option(USE_BRAINSFit                      "Build BRAINSFit"                      ON)
option(USE_BRAINSResample                 "Build BRAINSResample"                 ON)
option(USE_BRAINSROIAuto                  "Build BRAINSROIAuto"                  ON)
option(USE_DWIConvert                     "Build DWIConvert"                     ON)
option(USE_BRAINSLabelStats               "Build BRAINSLabelStats"               ON)
option(USE_BRAINSStripRotation            "Build BRAINSStripRotation"            ON)
option(USE_BRAINSTransformConvert         "Build BRAINSTransformConvert"         ON)
option(USE_BRAINSConstellationDetector    "Build BRAINSConstellationDetector"    ON)
option(USE_BRAINSInitializedControlPoints "Build BRAINSInitializedControlPoints" ON)
option(USE_BRAINSLandmarkInitializer      "Build BRAINSLandmarkInitializer"      ON)
option(USE_ImageCalculator                "Build ImageCalculator"                ON)
option(USE_ConvertBetweenFileFormats      "Build ConvertBetweenFileFormats"      ON)
option(USE_BRAINSDWICleanup               "Build BRAINSDWICleanup"               ON)
option(USE_BRAINSCreateLabelMapFromProbabilityMaps "Build BRAINSCreateLabelMapFromProbabilityMaps" OFF)
option(USE_BRAINSSnapShotWriter           "Build BRAINSSnapShotWriter"           ON)
if(${CMAKE_CXX_STANDARD} STREQUAL "11" OR ${CMAKE_CXX_STANDARD} STREQUAL "14")
  option(USE_BRAINSABC                      "Build BRAINSABC"                      ON)
else()
  option(USE_BRAINSABC                      "Build BRAINSABC"                      OFF)
endif()


## These are no longer needed on a day to day basis
if(NOT BUILD_FOR_DASHBOARD)
  set(BUILD_FOR_DASHBOARD OFF)
endif()
option(USE_BRAINSCut                      "Build BRAINSCut"                      ${BUILD_FOR_DASHBOARD})
option(USE_BRAINSMultiSTAPLE              "Build BRAINSMultiSTAPLE"              ${BUILD_FOR_DASHBOARD})
option(USE_BRAINSDemonWarp                "Build BRAINSDemonWarp "               ${BUILD_FOR_DASHBOARD})
option(USE_GTRACT                         "Build GTRACT"                         ${BUILD_FOR_DASHBOARD})
option(USE_BRAINSMush                     "Build BRAINSMush"                     ${BUILD_FOR_DASHBOARD})
option(USE_BRAINSMultiModeSegment         "Build BRAINSMultiModeSegment"         ${BUILD_FOR_DASHBOARD})
option(USE_BRAINSMultiSTAPLE              "Build BRAINSMultiSTAPLE"              ${BUILD_FOR_DASHBOARD})

## These are not yet ready for prime time.
option(USE_BRAINSTalairach                "Build BRAINSTalairach"                ${BUILD_FOR_DASHBOARD})
option(USE_BRAINSContinuousClass          "Build BRAINSContinuousClass "         OFF)
option(USE_BRAINSSurfaceTools             "Build BRAINSSurfaceTools     "        ${BUILD_FOR_DASHBOARD})
option(USE_ICCDEF                         "Build ICCDEF     "                    OFF)
option(USE_BRAINSPosteriorToContinuousClass             "Build BRAINSPosteriorToContinuousClass" OFF)
option(USE_DebugImageViewer "Build DebugImageViewer" OFF)
option(BRAINS_DEBUG_IMAGE_WRITE "Enable writing out intermediate image results" OFF)
option(USE_TBB "Build TBB as an internal module. This feature is still experimental and unsupported" OFF)
mark_as_advanced(USE_TBB)

if(Slicer_BUILD_BRAINSTOOLS OR USE_AutoWorkup OR USE_GTRACT OR USE_BRAINSTalairach OR USE_BRAINSSurfaceTools OR USE_BRAINSConstellationDetector OR USE_BRAINSDemonWarp OR USE_ConvertBetweenFileFormats )
  set(BRAINSTools_REQUIRES_VTK ON)
endif()

if(USE_ICCDEF OR ITK_USE_FFTWD OR ITK_USE_FFTWF)
  set(${PROJECT_NAME}_BUILD_FFTWF_SUPPORT ON)
endif()

if(${LOCAL_PROJECT_NAME}_USE_QT)
  if(NOT QT4_FOUND)
    find_package(Qt4 4.8 COMPONENTS QtCore QtGui QtNetwork QtXml REQUIRED)
    include(${QT_USE_FILE})
  endif()
endif()

#-----------------------------------------------------------------------------
# CMake Function(s) and Macro(s)
#-----------------------------------------------------------------------------
if(CMAKE_VERSION VERSION_LESS 2.8.3)
  include(Pre283CMakeParseArguments)
else()
  include(CMakeParseArguments)
endif()

#-----------------------------------------------------------------------------
# Platform check
#-----------------------------------------------------------------------------
set(PLATFORM_CHECK true)
if(PLATFORM_CHECK)
  # See CMake/Modules/Platform/Darwin.cmake)
  #   6.x == Mac OSX 10.2 (Jaguar)
  #   7.x == Mac OSX 10.3 (Panther)
  #   8.x == Mac OSX 10.4 (Tiger)
  #   9.x == Mac OSX 10.5 (Leopard)
  #  10.x == Mac OSX 10.6 (Snow Leopard)
  if (DARWIN_MAJOR_VERSION LESS "9")
    message(FATAL_ERROR "Only Mac OSX >= 10.5 are supported !")
  endif()
endif()


#-----------------------------------------------------------------------------
if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_INSTALL_LIBRARY_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_ARCHIVE_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_RUNTIME_DESTINATION bin)

#-------------------------------------------------------------------------
SETIFEMPTY(BRAINSTools_CLI_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SETIFEMPTY(BRAINSTools_CLI_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
SETIFEMPTY(BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#-------------------------------------------------------------------------
SETIFEMPTY(BRAINSTools_CLI_INSTALL_LIBRARY_DESTINATION ${CMAKE_INSTALL_LIBRARY_DESTINATION})
SETIFEMPTY(BRAINSTools_CLI_INSTALL_ARCHIVE_DESTINATION ${CMAKE_INSTALL_ARCHIVE_DESTINATION})
SETIFEMPTY(BRAINSTools_CLI_INSTALL_RUNTIME_DESTINATION ${CMAKE_INSTALL_RUNTIME_DESTINATION})

#-------------------------------------------------------------------------
# Augment compiler flags
#-------------------------------------------------------------------------
include(ITKSetStandardCompilerFlags)
if(ITK_LEGACY_REMOVE)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS} " )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS} " )
  else() # Release, or anything else
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS} " )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS} " )
  endif()
endif()
