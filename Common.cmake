#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

  mark_as_superbuild(VARS CMAKE_BUILD_TYPE ALL_PROJECTS)
endif()

#-----------------------------------------------------------------------------
if(APPLE)
#-----------------------------------------------------------------------------
# Platform check
#-----------------------------------------------------------------------------
  # See CMake/Modules/Platform/Darwin.cmake)
  #   6.x == Mac OSX 10.2 (Jaguar)
  #   7.x == Mac OSX 10.3 (Panther)
  #   8.x == Mac OSX 10.4 (Tiger)
  #   9.x == Mac OSX 10.5 (Leopard)
  #  10.x == Mac OSX 10.6 (Snow Leopard)
  #  11.x == Mac OSX 10.7 (Lion)
  #  12.x == Mac OSX 10.8 (Mountain Lion)
  if (DARWIN_MAJOR_VERSION LESS "13")
    message(FATAL_ERROR "Only Mac OSX >= 10.9 are supported !")
  endif()

  ## RPATH-RPATH-RPATH
  ## https://cmake.org/Wiki/CMake_RPATH_handling
  ## Always full RPATH
  # In many cases you will want to make sure that the required libraries are
  # always found independent from LD_LIBRARY_PATH and the install location. Then
  # you can use these settings:

  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif("${isSystemDir}" STREQUAL "-1")
  ## RPATH-RPATH-RPATH
endif()

#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)
include(itkCheckSourceTree)

include(CMakeDependentOption)
include(CMakeParseArguments)

#------------------------------------------------------------------------------
if(NOT Slicer_BUILD_BRAINSTOOLS)
   set(BUILD_SHARED_LIBS OFF) ## Build everything static for non-slicer builds
endif()
#------------------------------------------------------------------------------

# Enable this option to avoid unnecessary re-compilation associated with command line module
set(GENERATECLP_USE_MD5 ON)

#-----------------------------------------------------------------------------
# Build option(s)
#-----------------------------------------------------------------------------

option(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK "Determine if tools depending on VTK need to be built." ON)
mark_as_advanced(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)

cmake_dependent_option(${LOCAL_PROJECT_NAME}_USE_QT
      "Find and use Qt with VTK to build GUI Tools" OFF
      "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)

# bt_option: Convenience macro allowing to set an option and call mark_as_superbuild.
macro(bt_option name)
  option(${name} ${ARGN})
  mark_as_superbuild(VARS     ${name}:BOOL
                     PROJECTS ${LOCAL_PROJECT_NAME}
  )
endmacro()

bt_option(USE_ANTS                           "Build ANTS"                           ON)

bt_option(USE_BRAINSFit                      "Build BRAINSFit"                      ON)
bt_option(USE_BRAINSResample                 "Build BRAINSResample"                 ON)
bt_option(USE_BRAINSROIAuto                  "Build BRAINSROIAuto"                  ON)
bt_option(USE_DWIConvert                     "Build DWIConvert"                     ON)
bt_option(USE_BRAINSLabelStats               "Build BRAINSLabelStats"               ON)
bt_option(USE_BRAINSStripRotation            "Build BRAINSStripRotation"            ON)
bt_option(USE_BRAINSTransformConvert         "Build BRAINSTransformConvert"         ON)
bt_option(USE_BRAINSConstellationDetector    "Build BRAINSConstellationDetector"    ON)
bt_option(USE_BRAINSInitializedControlPoints "Build BRAINSInitializedControlPoints" ON)
bt_option(USE_BRAINSLandmarkInitializer      "Build BRAINSLandmarkInitializer"      ON)
bt_option(USE_ImageCalculator                "Build ImageCalculator"                ON)
bt_option(USE_ConvertBetweenFileFormats      "Build ConvertBetweenFileFormats"      ON)
bt_option(USE_BRAINSDWICleanup               "Build BRAINSDWICleanup"               ON)
bt_option(USE_BRAINSSnapShotWriter           "Build BRAINSSnapShotWriter"           ON)

cmake_dependent_option(USE_BRAINSConstellationDetectorGUI
              "Build BRAINSConstellationDetectorGUI" OFF
              "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_BRAINSConstellationDetectorGUI:BOOL
                   PROJECTS ${LOCAL_PROJECT_NAME} )

bt_option(USE_AutoWorkup                     "Build AutoWorkup"                     ON)
bt_option(USE_ReferenceAtlas                 "Build the Reference Atlas"            ON)
cmake_dependent_option(USE_BRAINSABC
               "Build BRAINSABC" ON
               "USE_AutoWorkup;USE_ReferenceAtlas" OFF)
mark_as_superbuild(VARS USE_BRAINSABC:BOOL
                   PROJECTS ${LOCAL_PROJECT_NAME} )

## These are no longer needed on a day to day basis
if(NOT BUILD_FOR_DASHBOARD)
  set(BUILD_FOR_DASHBOARD OFF)
endif()
bt_option(USE_BRAINSCreateLabelMapFromProbabilityMaps "Build BRAINSCreateLabelMapFromProbabilityMaps" ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSSuperResolution          "Build BRAINSSuperResolution"          ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSCut                      "Build BRAINSCut"                      ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSMultiSTAPLE              "Build BRAINSMultiSTAPLE"              ${BUILD_FOR_DASHBOARD})

cmake_dependent_option(USE_BRAINSDemonWarp "Build BRAINSDemonWarp " ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_BRAINSDemonWarp:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
cmake_dependent_option(USE_GTRACT "Build GTRACT" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_GTRACT:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
cmake_dependent_option(USE_ITKMatlabIO "Build ITKMatlabIO" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_ITKMatlabIO:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
bt_option(USE_BRAINSMush                     "Build BRAINSMush"                     ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSMultiModeSegment         "Build BRAINSMultiModeSegment"         ${BUILD_FOR_DASHBOARD})

## These are not yet ready for prime time.
## TODO: Move to ARCHIVE directory
bt_option(BUILD_ARCHIVE                      "Build old tools from archive"        OFF)
cmake_dependent_option(USE_BRAINSSurfaceTools "Build BRAINSSurfaceTools" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_BRAINSSurfaceTools:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
bt_option(USE_BRAINSPosteriorToContinuousClass             "Build BRAINSPosteriorToContinuousClass" OFF)
cmake_dependent_option(USE_DebugImageViewer "Build DebugImageViewer" OFF "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_DebugImageViewer:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
bt_option(BRAINS_DEBUG_IMAGE_WRITE "Enable writing out intermediate image results" OFF)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Items that are archived.  May be usefult for compiler testing,
## but probably wont be useful for research work.
cmake_dependent_option(USE_BRAINSContinuousClass          "Build BRAINSContinuousClass" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK;BUILD_ARCHIVE" ON)
mark_as_superbuild(VARS USE_BRAINSContinuousClass:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )

cmake_dependent_option(USE_ICCDEF "Build ICCDEF" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK;BUILD_ARCHIVE" OFF)
mark_as_superbuild(VARS USE_ICCDEF:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )

cmake_dependent_option(USE_BRAINSRefacer "BRAINSRefacer is still under development." ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK;BUILD_ARCHIVE" ON)
mark_as_superbuild(VARS USE_BRAINSRefacer:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )

cmake_dependent_option(USE_BRAINSTalairach "Build BRAINSTalairach is under development" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK;BUILD_ARCHIVE" ON)
mark_as_superbuild(VARS USE_BRAINSTalairach:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )

cmake_dependent_option(USE_DicomSignature "Build DicomSignature is under development" ${BUILD_FOR_DASHBOARD} "BUILD_ARCHIVE" ON)
mark_as_superbuild(VARS USE_DicomSignature:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#if(NOT ${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)
#  message("NOTE: Following toolkits are dependent to VTK:
#      -GTRACT
#      -BRAINSDemonWarp
#      -BRAINSTalairach
#      -BRAINSSurfaceTools
#      -DebugImageViewer
#      -ITKMatlabIO
#      -BRAINSConstellationDetectorGUI
#      First you need to set ${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK to ON to be able to choose above application for build.")
#endif()

if(${LOCAL_PROJECT_NAME}_USE_QT) #//TODO:  BRAINSTools only indirectly needs QT!,
  set(${LOCAL_PROJECT_NAME}_REQUIRED_QT_MODULES
    Core Widgets
    Multimedia
    Network OpenGL
    PrintSupport # Required by "Annotations" module
    UiTools #no dll
    Xml XmlPatterns
    Svg Sql
    )
  find_package(Qt5 COMPONENTS Core QUIET)
  if(Qt5_VERSION VERSION_LESS "5.6.0")
    list(APPEND ${LOCAL_PROJECT_NAME}_REQUIRED_QT_MODULES
      WebKit
      )
  else()
    list(APPEND ${LOCAL_PROJECT_NAME}_REQUIRED_QT_MODULES
      WebEngine
      WebEngineWidgets
      WebChannel
      )
  endif()
  if(${LOCAL_PROJECT_NAME}_BUILD_EXTENSIONMANAGER_SUPPORT)
    list(APPEND ${LOCAL_PROJECT_NAME}_REQUIRED_QT_MODULES Script)
  endif()
  if(${LOCAL_PROJECT_NAME}_BUILD_I18N_SUPPORT)
    list(APPEND ${LOCAL_PROJECT_NAME}_REQUIRED_QT_MODULES LinguistTools) # no dll
  endif()
  if(BUILD_TESTING)
    list(APPEND ${LOCAL_PROJECT_NAME}_REQUIRED_QT_MODULES Test)
  endif()
  find_package(Qt5 COMPONENTS ${${LOCAL_PROJECT_NAME}_REQUIRED_QT_MODULES})
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
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS}" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS}")
  else() # Release, or anything else
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS}" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS}" )
  endif()
endif()

#-----------------------------------------------------------------------------
# Add needed flag for gnu on linux like enviroments to build static common libs
# suitable for linking with shared object libs.
#message("CMAKE_SYSTEM_PROCESSOR:${CMAKE_SYSTEM_PROCESSOR}")
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  # message("Adding fPIC")
  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
  if(NOT "${CMAKE_C_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  endif()
endif()

option(BUILD_OPTIMIZED "Set compiler flags for native host building." OFF)
mark_as_advanced(BUILD_OPTIMIZED)

if(BUILD_OPTIMIZED)
  include(CheckCXXCompilerFlag)
  include(CheckCCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" CAN_BUILD_CXX_OPTIMIZED)
  CHECK_C_COMPILER_FLAG("-march=native" CAN_BUILD_C_OPTIMIZED)
  if(CAN_BUILD_CXX_OPTIMIZED AND CAN_BUILD_C_OPTIMIZED)
    string(FIND CMAKE_CXX_FLAGS "-march=native" PREVIOUS_CXX_OPTIMIZED_SET)
    if(PREVIOUS_CXX_OPTIMIZED_SET LESS 0)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
    string(FIND CMAKE_C_FLAGS "-march=native" PREVIOUS_C_OPTIMIZED_SET)
    if(PREVIOUS_C_OPTIMIZED_SET LESS 0 )
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    endif()
  else()
    message("WARNING: Requested optimized build, but -march=native flag not supported by"
            "${CMAKE_CXX_COMPILER}"
            "${CMAKE_C_COMPILER}")
  endif()
endif()
mark_as_superbuild(VARS BUILD_OPTIMIZED:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )

#-------------------------------------------------------------------------
if(NOT DEFINED BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET)
  set(BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET "BRAINSToolsFetchData")
endif()
mark_as_superbuild(VARS BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET:STRING PROJECTS ${LOCAL_PROJECT_NAME} )

#------------------------------------------------------------------------------
# Configure and build
#------------------------------------------------------------------------------
## BRAINSTools_MAX_TEST_LEVEL adjusts how agressive the test suite is
## so that long running tests or incomplete tests can easily be
## silenced
## 1 - Run the absolute minimum very fast tests (These should always pass before any code commit)
## 3 - Run fast tests on continous builds (These need immediate attention if they begin to fail)
## 5 - Run moderate nightly tests (These need immediate attention if they begin to fail)
## 7 - Run long running extensive test that are a burden to normal development (perhaps test 1x per week)
## 8 - Run tests that fail due to incomplete test building, these are good ideas for test that we don't have time to make robust)
## 9 - Run silly tests that don't have much untility
set(BRAINSTools_MAX_TEST_LEVEL 3 CACHE STRING "Testing level for managing test burden")
