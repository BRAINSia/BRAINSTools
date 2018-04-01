#-----------------------------------------------------------------------------
include(ExternalProject)
include(ExternalProjectDependency) #<-- Must be after project() calls
include(ExternalProjectGenerateProjectDescription)
#-----------------------------------------------------------------------------
if(APPLE)
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

  # Note: By setting CMAKE_OSX_* variables before any enable_language() or project() calls,
  #       we ensure that the bitness, and C++ standard library will be properly detected.
  include(BlockSetCMakeOSXVariables)
  mark_as_superbuild(
    VARS CMAKE_OSX_ARCHITECTURES:STRING CMAKE_OSX_SYSROOT:PATH CMAKE_OSX_DEPLOYMENT_TARGET:STRING
    ALL_PROJECTS
    )
endif()

#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)
include(itkCheckSourceTree)

include(CMakeDependentOption)
include(CMakeParseArguments)

if(NOT Slicer_BUILD_BRAINSTOOLS)
   set(BUILD_SHARED_LIBS OFF) ## Build everything static for non-slicer builds
endif()
#------------------------------------------------------------------------------
#if(Slicer_BUILD_BRAINSTOOLS OR USE_AutoWorkup OR USE_GTRACT OR USE_BRAINSTalairach OR USE_BRAINSSurfaceTools OR USE_BRAINSConstellationDetector OR USE_BRAINSDemonWarp OR USE_ConvertBetweenFileFormats )

## VTK is not easy to build on all platforms
if(Slicer_BUILD_BRAINSTOOLS)
  option(${PRIMARY_PROJECT_NAME}_REQUIRES_VTK "Determine if tools depending on VTK need to be built." ON)
else()
  option(${PRIMARY_PROJECT_NAME}_REQUIRES_VTK "Determine if tools depending on VTK need to be built." OFF)
  # Enable this option to avoid unnecessary re-compilation associated with command line module
  set(GENERATECLP_USE_MD5 ON)
endif()
mark_as_advanced(${PRIMARY_PROJECT_NAME}_REQUIRES_VTK)

option(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT "Install development support include and libraries for external packages." OFF)
mark_as_advanced(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)

CMAKE_DEPENDENT_OPTION(${LOCAL_PROJECT_NAME}_USE_QT "Find and use Qt with VTK to build GUI Tools" OFF "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK" OFF)

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
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()
if(NOT CMAKE_CONFIGURATION_TYPES)
  mark_as_superbuild(VARS CMAKE_BUILD_TYPE ALL_PROJECTS)
endif()

if(${ITK_VERSION_MAJOR} STREQUAL "3")
  message(FATAL_ERROR "ITKv3 is no longer supported")
endif()

#-----------------------------------------------------------------------------
# Build option(s)
#-----------------------------------------------------------------------------

# bt_option: Convenience macro allowing to set an option and call mark_as_superbuild.
macro(bt_option name)
  option(${name} ${ARGN})
  mark_as_superbuild(${name})
endmacro()

bt_option(USE_AutoWorkup                     "Build AutoWorkup"                     ON)
bt_option(USE_ReferenceAtlas                 "Build the Reference Atlas"            ON)

bt_option(USE_ANTS                           "Build ANTS"                           ON)
bt_option(BUILD_ARCHIVE                      "Build old tools from archive"        OFF)

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
bt_option(USE_BRAINSCreateLabelMapFromProbabilityMaps "Build BRAINSCreateLabelMapFromProbabilityMaps" OFF)
bt_option(USE_BRAINSSnapShotWriter           "Build BRAINSSnapShotWriter"           ON)
bt_option(USE_BRAINSSuperResolution          "Build BRAINSSuperResolution"          OFF)

cmake_dependent_option(USE_BRAINSConstellationDetectorGUI "Build BRAINSConstellationDetectorGUI" OFF "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK" OFF)
mark_as_superbuild(USE_BRAINSConstellationDetectorGUI)

cmake_dependent_option(USE_BRAINSABC "Build BRAINSABC" ON "USE_AutoWorkup;USE_ReferenceAtlas" OFF)
mark_as_superbuild(USE_BRAINSABC)

## These are no longer needed on a day to day basis
if(NOT BUILD_FOR_DASHBOARD)
  set(BUILD_FOR_DASHBOARD OFF)
endif()
bt_option(USE_BRAINSCut                      "Build BRAINSCut"                      ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSMultiSTAPLE              "Build BRAINSMultiSTAPLE"              ${BUILD_FOR_DASHBOARD})
cmake_dependent_option(USE_BRAINSDemonWarp "Build BRAINSDemonWarp " ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK" OFF)
mark_as_superbuild(USE_BRAINSDemonWarp)
cmake_dependent_option(USE_GTRACT "Build GTRACT" ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK" OFF)
mark_as_superbuild(USE_GTRACT)
cmake_dependent_option(USE_ITKMatlabIO "Build ITKMatlabIO" ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK" OFF)
mark_as_superbuild(USE_ITKMatlabIO)
bt_option(USE_BRAINSMush                     "Build BRAINSMush"                     ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSMultiModeSegment         "Build BRAINSMultiModeSegment"         ${BUILD_FOR_DASHBOARD})

## These are not yet ready for prime time.
cmake_dependent_option(USE_BRAINSSurfaceTools "Build BRAINSSurfaceTools" ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK" OFF)
mark_as_superbuild(USE_BRAINSSurfaceTools)
bt_option(USE_BRAINSPosteriorToContinuousClass             "Build BRAINSPosteriorToContinuousClass" OFF)
cmake_dependent_option(USE_DebugImageViewer "Build DebugImageViewer" OFF "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK" OFF)
mark_as_superbuild(USE_DebugImageViewer)
bt_option(BRAINS_DEBUG_IMAGE_WRITE "Enable writing out intermediate image results" OFF)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Items that are archived.  May be usefult for compiler testing,
## but probably wont be useful for research work.
cmake_dependent_option(USE_BRAINSContinuousClass          "Build BRAINSContinuousClass" ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK;BUILD_ARCHIVE" ON)
mark_as_superbuild(USE_BRAINSContinuousClass)

cmake_dependent_option(USE_ICCDEF "Build ICCDEF" ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK;BUILD_ARCHIVE" OFF)
mark_as_superbuild(USE_ICCDEF)

cmake_dependent_option(USE_BRAINSRefacer "BRAINSRefacer is still under development." ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK;BUILD_ARCHIVE" ON)
mark_as_superbuild(USE_BRAINSRefacer)

cmake_dependent_option(USE_BRAINSTalairach "Build BRAINSTalairach is under development" ${BUILD_FOR_DASHBOARD} "${PRIMARY_PROJECT_NAME}_REQUIRES_VTK;BUILD_ARCHIVE" ON)
mark_as_superbuild(USE_BRAINSTalairach)

cmake_dependent_option(USE_DicomSignature "Build DicomSignature is under development" ${BUILD_FOR_DASHBOARD} "BUILD_ARCHIVE" ON)
mark_as_superbuild(USE_DicomSignature)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#if(NOT ${PRIMARY_PROJECT_NAME}_REQUIRES_VTK)
#  message("NOTE: Following toolkits are dependent to VTK:
#      -GTRACT
#      -BRAINSDemonWarp
#      -BRAINSTalairach
#      -BRAINSSurfaceTools
#      -DebugImageViewer
#      -ITKMatlabIO
#      -BRAINSConstellationDetectorGUI
#      First you need to set ${PRIMARY_PROJECT_NAME}_REQUIRES_VTK to ON to be able to choose above application for build.")
#endif()

if(${LOCAL_PROJECT_NAME}_USE_QT)
  if(NOT QT4_FOUND)
    find_package(Qt4 4.8 COMPONENTS QtCore QtGui QtNetwork QtXml REQUIRED)
    include(${QT_USE_FILE})
  endif()
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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS} " )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS} " )
  else() # Release, or anything else
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS} " )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS} " )
  endif()
endif()

#-------------------------------------------------------------------------
if(NOT DEFINED BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET)
  set(BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET "BRAINSToolsFetchData")
endif()
mark_as_superbuild(BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET:STRING)
