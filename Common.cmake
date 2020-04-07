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
endif()
#-----------------------------------------------------------------------------
# Set a default external project build type if none was specified
set(EXTERNAL_PROJECT_BUILD_TYPE "Release" CACHE STRING "Default build type for support libraries")
set_property(CACHE EXTERNAL_PROJECT_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

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
  #  13.x == Mac OSX 10.9 (Yosemite)
  #  14.x == Mac OSX 10.10 (El Capitan)
  #  15.x == Mac OSX 10.12 (Sierra)    # Sept 2016 -- Improve C++11 support by default, for TBB
  #  17.x == Mac OSX 10.13 (High Sierra)
  #  18.x == Mac OSX 10.14 (Mojave)
  if (DARWIN_MAJOR_VERSION LESS "13")  #https://en.wikipedia.org/wiki/Darwin_(operating_system)
    message(FATAL_ERROR "Only Mac OSX >= 10.13 are supported !")
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
  set(CMAKE_MACOSX_RPATH 0)
  mark_as_superbuild(
    VARS
      CMAKE_SKIP_BUILD_RPATH:STRING
      CMAKE_BUILD_WITH_INSTALL_RPATH:STRING
      CMAKE_INSTALL_RPATH:STRING
      CMAKE_INSTALL_RPATH_USE_LINK_PATH:STRING
      CMAKE_OSX_ARCHITECTURES:STRING
      CMAKE_OSX_SYSROOT:PATH
      CMAKE_OSX_DEPLOYMENT_TARGET:STRING
      CMAKE_MACOSX_RPATH:BOOL
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

# -- if(NOT Slicer_BUILD_BRAINSTOOLS AND BUILD_TESTING)
# --   message(FATAL_ERROR "USE ${LOCAL_PROJECT_NAME}_BUILD_TESTING for controlling brainstools testing")
# -- endif()
option(${LOCAL_PROJECT_NAME}_BUILD_TESTING "Control if BRAINSTools builds testing" ON)

option(${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT "Build Dicom Support"    OFF)
mark_as_advanced(${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT)
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
bt_option(USE_DWIConvert                     "Build DWIConvert"                     ${${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT})
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

if(USE_DWIConvert AND NOT ${${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT})
  message(FATAL_ERROR "USE_DWIConvert requires ${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT:BOOL=ON")
endif()
## These are no longer needed on a day to day basis
if(NOT BUILD_FOR_DASHBOARD)
  set(BUILD_FOR_DASHBOARD OFF)
endif()
bt_option(USE_BRAINSCreateLabelMapFromProbabilityMaps "Build BRAINSCreateLabelMapFromProbabilityMaps" ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSSuperResolution          "Build BRAINSSuperResolution"          ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSMultiSTAPLE              "Build BRAINSMultiSTAPLE"              ${BUILD_FOR_DASHBOARD})

message(STATUS "BUILD_FOR_DASHBOARD: ${BUILD_FOR_DASHBOARD}")
message(STATUS "${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT: ${${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT}")
cmake_dependent_option(USE_GTRACT "Build GTRACT" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" ${${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT})
mark_as_superbuild(VARS USE_GTRACT:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
cmake_dependent_option(USE_ITKMatlabIO "Build ITKMatlabIO" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_ITKMatlabIO:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
bt_option(USE_BRAINSMush                     "Build BRAINSMush"                     ${BUILD_FOR_DASHBOARD})
bt_option(USE_BRAINSMultiModeSegment         "Build BRAINSMultiModeSegment"         ${BUILD_FOR_DASHBOARD})

## These are not yet ready for prime time.
## INFO: Move to ARCHIVE directory
bt_option(BUILD_ARCHIVE                                    "Build old tools from archive"        OFF)
bt_option(USE_BRAINSPosteriorToContinuousClass             "Build BRAINSPosteriorToContinuousClass" OFF)

cmake_dependent_option(USE_DebugImageViewer "Build DebugImageViewer" OFF "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK" OFF)
mark_as_superbuild(VARS USE_DebugImageViewer:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )
bt_option(BRAINS_DEBUG_IMAGE_WRITE "Enable writing out intermediate image results" OFF)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Items that are archived.  May be usefult for compiler testing,
## but probably wont be useful for research work.

cmake_dependent_option(USE_BRAINSTalairach "Build BRAINSTalairach is under development" ${BUILD_FOR_DASHBOARD} "${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK;BUILD_ARCHIVE" OFF)
mark_as_superbuild(VARS USE_BRAINSTalairach:BOOL PROJECTS ${LOCAL_PROJECT_NAME} )

#if(NOT ${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK)
#  message("NOTE: Following toolkits are dependent to VTK:
#      -GTRACT
#      -BRAINSTalairach
#      -DebugImageViewer
#      -ITKMatlabIO
#      -BRAINSConstellationDetectorGUI
#      First you need to set ${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_VTK to ON to be able to choose above application for build.")
#endif()

if(${LOCAL_PROJECT_NAME}_USE_QT) #//INFO:  BRAINSTools only indirectly needs QT!,
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
  if(${LOCAL_PROJECT_NAME}_BUILD_TESTING)
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
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-${CMAKE_BUILD_TYPE}-EP${EXTERNAL_PROJECT_BUILD_TYPE}-build/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-${CMAKE_BUILD_TYPE}-EP${EXTERNAL_PROJECT_BUILD_TYPE}-build/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-${CMAKE_BUILD_TYPE}-EP${EXTERNAL_PROJECT_BUILD_TYPE}-build/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-${CMAKE_BUILD_TYPE}-EP${EXTERNAL_PROJECT_BUILD_TYPE}-build/bin)
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

string(APPEND CMAKE_C_FLAGS  " ${ITK_REQUIRED_C_FLAGS}")
string(APPEND CMAKE_CXX_FLAGS " ${ITK_REQUIRED_CXX_FLAGS}")
string(APPEND CMAKE_EXE_LINKER_FLAGS " ${ITK_REQUIRED_LINK_FLAGS}")
string(APPEND CMAKE_SHARED_LINKER_FLAGS " ${ITK_REQUIRED_LINK_FLAGS}")
string(APPEND CMAKE_MODULE_LINKER_FLAGS " ${ITK_REQUIRED_LINK_FLAGS}")

#-----------------------------------------------------------------------------
# Add needed flag for gnu on linux like enviroments to build static common libs
# suitable for linking with shared object libs.
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "-fPIC")
    string(APPEND CMAKE_CXX_FLAGS " -fPIC")
    string(APPEND ep_CMAKE_CXX_FLAGS " -fPIC")
  endif()
  if(NOT "${CMAKE_C_FLAGS}" MATCHES "-fPIC")
    string(APPEND CMAKE_C_FLAGS " -fPIC")
    string(APPEND ep_CMAKE_C_FLAGS " -fPIC")
  endif()
endif()

option(BUILD_COVERAGE "Build ${LOCAL_PROJECT_NAME} for coverage analysis" OFF)
mark_as_advanced(BUILD_COVERAGE)

# ------ duplicate from ITK
function(check_compiler_optimization_flags c_optimization_flags_var cxx_optimization_flags_var)
  set(${c_optimization_flags_var} "" PARENT_SCOPE)
  set(${cxx_optimization_flags_var} "" PARENT_SCOPE)

  if(MSVC)
    check_avx_flags(InstructionSetOptimizationFlags)
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "4")
       list(APPEND InstructionSetOptimizationFlags
            /arch:SSE /arch:SSE2)
    endif()
  elseif(NOT EMSCRIPTEN)
    if (${CMAKE_C_COMPILER} MATCHES "icc.*$")
      set(USING_INTEL_ICC_COMPILER TRUE)
    endif()
    if (${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
      set(USING_INTEL_ICC_COMPILER TRUE)
    endif()
    if(USING_INTEL_ICC_COMPILER)
      set(InstructionSetOptimizationFlags "")
    else()
      set(InstructionSetOptimizationFlags "")
    endif ()

    # Check this list on C compiler only
    set(c_flags "")

    # Check this list on C++ compiler only
    set(cxx_flags "")

    # Check this list on both C and C++ compilers
    set(InstructionSetOptimizationFlags
       # https://gcc.gnu.org/onlinedocs/gcc-4.8.0/gcc/i386-and-x86_002d64-Options.html
       # NOTE the corei7 release date was 2008
       -mtune=native # Tune the code for the computer used compile ITK, but allow running on generic cpu archetectures
       ## use corei7-avx for ARGON cluster builds
       -march=native # Use ABI settings to support corei7 (circa 2008 ABI feature sets, corei7-avx circa 2013)
       )
  endif()
  set(c_and_cxx_flags ${InstructionSetOptimizationFlags})

  check_c_compiler_flags(    CMAKE_C_WARNING_FLAGS ${c_and_cxx_flags} ${c_flags})
  check_cxx_compiler_flags(CMAKE_CXX_WARNING_FLAGS ${c_and_cxx_flags} ${cxx_flags})

  set(${c_optimization_flags_var} "${CMAKE_C_WARNING_FLAGS}" PARENT_SCOPE)
  set(${cxx_optimization_flags_var} "${CMAKE_CXX_WARNING_FLAGS}" PARENT_SCOPE)
endfunction()
if(NOT BRAINSToools_C_OPTIMIZATION_FLAGS OR NOT BRAINSToools_CXX_OPTIMIZATION_FLAGS ) # Only check once if not explicitly set on command line
  #-----------------------------------------------------------------------------
  #Check the set of warning flags the compiler supports
  check_compiler_optimization_flags(C_OPTIMIZATION_FLAGS CXX_OPTIMIZATION_FLAGS)
endif()
if(NOT BRAINSToools_C_OPTIMIZATION_FLAGS) #Not set on cmake command line option -DBRAINSToools_C_OPTIMIZATION_FLAGS:STRING=""
  set(BRAINSToools_C_OPTIMIZATION_FLAGS ${C_OPTIMIZATION_FLAGS} CACHE STRING "ITK C Compiler ABI/Optimization flags, Use '-march=native' to maximize performance, but break portabilitly.")
else()
  set(BRAINSToools_C_OPTIMIZATION_FLAGS ${BRAINSToools_C_OPTIMIZATION_FLAGS} CACHE STRING "ITK C Compiler ABI/Optimization flags, Use '-march=native' to maximize performance, but break portabilitly.")
endif()
if(NOT BRAINSToools_CXX_OPTIMIZATION_FLAGS) #Not set on cmake command line option -DBRAINSToools_CXX_OPTIMIZATION_FLAGS:STRING=""
  set(BRAINSToools_CXX_OPTIMIZATION_FLAGS ${CXX_OPTIMIZATION_FLAGS} CACHE STRING "ITK CXX Compiler ABI/Optimization flags, Use '-march=native' to maximize performance, but break portabilitly.")
else()
  set(BRAINSToools_CXX_OPTIMIZATION_FLAGS ${BRAINSToools_CXX_OPTIMIZATION_FLAGS} CACHE STRING "ITK CXX Compiler ABI/Optimization flags, Use '-march=native' to maximize performance, but break portabilitly.")
endif()
mark_as_advanced(BRAINSToools_CXX_OPTIMIZATION_FLAGS)
mark_as_advanced(BRAINSToools_C_OPTIMIZATION_FLAGS)
unset(C_OPTIMIZATION_FLAGS)
unset(CXX_OPTIMIZATION_FLAGS)

option(BUILD_OPTIMIZED "Set compiler flags for native host building." OFF)
mark_as_advanced(BUILD_OPTIMIZED)

if(BUILD_OPTIMIZED AND NOT BUILD_COVERAGE)
  if(CAN_BUILD_CXX_OPTIMIZED AND CAN_BUILD_C_OPTIMIZED)
    string(APPEND CMAKE_CXX_FLAGS ${BRAINSToools_CXX_OPTIMIZATION_FLAGS})
    string(APPEND CMAKE_C_FLAGS ${BRAINSToools_CXX_OPTIMIZATION_FLAGS})
    string(REPLACE " " ";" CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${${PROJECT_NAME}_C_OPTIMIZATION_FLAGS} ${${PROJECT_NAME}_C_WARNING_FLAGS}")
    list(REMOVE_DUPLICATES CMAKE_C_FLAGS)
    string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

    string(REPLACE " " ";" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${PROJECT_NAME}_CXX_OPTIMIZATION_FLAGS} ${${PROJECT_NAME}_CXX_WARNING_FLAGS}")
    list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS)
    string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    message("WARNING: Requested optimized build, but -march=native flag not supported by"
            "${CMAKE_CXX_COMPILER}"
            "${CMAKE_C_COMPILER}")
  endif()
endif()

if(NOT Slicer_BUILD_BRAINSTOOLS)

## Let Slicer use it
mark_as_superbuild(
   VARS
#--    CMAKE_BUILD_TYPE:STRING
#--    CMAKE_CXX_FLAGS_DEBUG:STRING
#--    CMAKE_CXX_FLAGS_MINSIZEREL:STRING
#--    CMAKE_CXX_FLAGS_RELEASE:STRING
#--    CMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING
#--    CMAKE_C_FLAGS_DEBUG:STRING
#--    CMAKE_C_FLAGS_MINSIZEREL:STRING
#--    CMAKE_C_FLAGS_RELEASE:STRING
#--    CMAKE_C_FLAGS_RELWITHDEBINFO:STRING
#--    CMAKE_EXE_LINKER_FLAGS:STRING
#--    CMAKE_EXE_LINKER_FLAGS_DEBUG:STRING
#--    CMAKE_EXE_LINKER_FLAGS_MINSIZEREL:STRING
#--    CMAKE_EXE_LINKER_FLAGS_RELEASE:STRING
#--    CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO:STRING
#--    CMAKE_EXTRA_GENERATOR:STRING
#--    CMAKE_GENERATOR:STRING
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH
    CMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH
    CMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH
    CMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH
#--    CMAKE_MODULE_PATH:PATH
#--    CMAKE_SHARED_LINKER_FLAGS:STRING
#--    CMAKE_SHARED_LINKER_FLAGS_DEBUG:STRING
#--    CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL:STRING
#--    CMAKE_SHARED_LINKER_FLAGS_RELEASE:STRING
#--    CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO:STRING
#--    CMAKE_SKIP_RPATH:BOOL
#--    CTEST_NEW_FORMAT:BOOL
#--    MEMORYCHECK_COMMAND:PATH
#--    MEMORYCHECK_COMMAND_OPTIONS:STRING
  ALL_PROJECTS
)
endif()

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
