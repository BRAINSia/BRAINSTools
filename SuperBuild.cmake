#-----------------------------------------------------------------------------
enable_language(C)
enable_language(CXX)
#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)

#-----------------------------------------------------------------------------
include(SlicerMacroGetOperatingSystemArchitectureBitness)

#-----------------------------------------------------------------------------
# Where should the superbuild source files be downloaded to?
# By keeping this outside of the build tree, you can share one
# set of external source trees for multiple build trees
#-----------------------------------------------------------------------------
# set( SOURCE_DOWNLOAD_CACHE ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
set( SOURCE_DOWNLOAD_CACHE ${CMAKE_CURRENT_BINARY_DIR} )

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)
#-----------------------------------------------------------------------------
# CTestCustom
#-----------------------------------------------------------------------------
if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  configure_file(
    CMake/CTestCustom.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake
    @ONLY)
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

#-----------------------------------------------------------------------------
# Git protocole option
#-----------------------------------------------------------------------------
option(${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT ${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

CMAKE_DEPENDENT_OPTION(${CMAKE_PROJECT_NAME}_USE_CTKAPPLAUNCHER "CTKAppLauncher used with python" ON
  "NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_python" OFF)

find_package(Git REQUIRED)

# I don't know who removed the Find_Package for QT, but it needs to be here
# in order to build VTK if ${LOCAL_PROJECT_NAME}_USE_QT is set.
if(${LOCAL_PROJECT_NAME}_USE_QT)
find_package(Qt4 REQUIRED)
endif()


#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------
include(ExternalProject)
include(SlicerMacroEmptyExternalProject)
include(SlicerMacroCheckExternalProjectDependency)

# Compute -G arg for configuring external projects with the same CMake generator:
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()


# With CMake 2.8.9 or later, the UPDATE_COMMAND is required for updates to occur.
# For earlier versions, we nullify the update state to prevent updates and
# undesirable rebuild.
option(FORCE_EXTERNAL_BUILDS "Force rebuilding of external project (if they are updated)" ON)
if(CMAKE_VERSION VERSION_LESS 2.8.9 OR NOT FORCE_EXTERNAL_BUILDS)
  set(cmakeversion_external_update UPDATE_COMMAND)
  set(cmakeversion_external_update_value "" )
else()
  set(cmakeversion_external_update LOG_UPDATE )
  set(cmakeversion_external_update_value 1)
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
# Superbuild option(s)
#-----------------------------------------------------------------------------
option(BUILD_STYLE_UTILS "Build uncrustify, cppcheck, & KWStyle" OFF)
CMAKE_DEPENDENT_OPTION(
  USE_SYSTEM_Uncrustify "Use system Uncrustify program" OFF
  "BUILD_STYLE_UTILS" OFF
  )
CMAKE_DEPENDENT_OPTION(
  USE_SYSTEM_KWStyle "Use system KWStyle program" OFF
  "BUILD_STYLE_UTILS" OFF
  )
CMAKE_DEPENDENT_OPTION(
  USE_SYSTEM_Cppcheck "Use system Cppcheck program" OFF
  "BUILD_STYLE_UTILS" OFF
  )

set(EXTERNAL_PROJECT_BUILD_TYPE "Release" CACHE STRING "Default build type for support libraries")

option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
option(USE_SYSTEM_DCMTK "Build using an externally defined version of DCMTK" OFF)
option(${PROJECT_NAME}_BUILD_DICOM_SUPPORT "Build Dicom Support" ON)

#------------------------------------------------------------------------------
# ${LOCAL_PROJECT_NAME} dependency list
#------------------------------------------------------------------------------

set(${LOCAL_PROJECT_NAME}_DEPENDENCIES DCMTK ITKv4 SlicerExecutionModel)

list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES teem)
if(USE_TBB) ## Needed for BRAINSABC
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES TBB)
endif()
#list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES Boost)

if(BUILD_STYLE_UTILS)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES Cppcheck KWStyle Uncrustify)
endif()

if(USE_BRAINSABC OR USE_BRAINSCut)
  #list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES qhull)
endif()

if(${PRIMARY_PROJECT_NAME}_REQUIRES_VTK)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES VTK)
endif()

if(USE_BRAINSCut)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES OpenCV)
endif()

if(USE_ANTS)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES ANTs)
endif()

if(USE_BRAINSSuperResolution)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES RTK)
endif()

if(USE_BRAINSSuperResolution)
mark_as_superbuild(
  VARS
  RTK_DIR:PATH
ALL_PROJECTS
)
endif()
#-----------------------------------------------------------------------------
# Common external projects CMake variables
#-----------------------------------------------------------------------------
mark_as_superbuild(
  VARS
    MAKECOMMAND:STRING
    CMAKE_SKIP_RPATH:BOOL
    BUILD_SHARED_LIBS:BOOL
    CMAKE_MODULE_PATH:PATH
    CMAKE_BUILD_TYPE:STRING
    # BUILD_SHARED_LIBS:BOOL
    CMAKE_INCLUDE_DIRECTORIES_BEFORE:BOOL
    CMAKE_CXX_COMPILER:PATH
    CMAKE_CXX_FLAGS:STRING
    CMAKE_CXX_FLAGS_DEBUG:STRING
    CMAKE_CXX_FLAGS_MINSIZEREL:STRING
    CMAKE_CXX_FLAGS_RELEASE:STRING
    CMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_C_COMPILER:PATH
    CMAKE_C_FLAGS:STRING
    CMAKE_C_FLAGS_DEBUG:STRING
    CMAKE_C_FLAGS_MINSIZEREL:STRING
    CMAKE_C_FLAGS_RELEASE:STRING
    CMAKE_C_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_SHARED_LINKER_FLAGS:STRING
    CMAKE_GENERATOR:STRING
    CMAKE_EXTRA_GENERATOR:STRING
    CMAKE_INSTALL_PREFIX:PATH
    CMAKE_EXE_LINKER_FLAGS:STRING
    CMAKE_EXE_LINKER_FLAGS_DEBUG:STRING
    CMAKE_EXE_LINKER_FLAGS_MINSIZEREL:STRING
    CMAKE_EXE_LINKER_FLAGS_RELEASE:STRING
    CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_MODULE_LINKER_FLAGS:STRING
    CMAKE_MODULE_LINKER_FLAGS_DEBUG:STRING
    CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL:STRING
    CMAKE_MODULE_LINKER_FLAGS_RELEASE:STRING
    CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_SHARED_LINKER_FLAGS:STRING
    CMAKE_SHARED_LINKER_FLAGS_DEBUG:STRING
    CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL:STRING
    CMAKE_SHARED_LINKER_FLAGS_RELEASE:STRING
    CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH
    CMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH
    CMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH
    CTEST_NEW_FORMAT:BOOL
    MEMORYCHECK_COMMAND_OPTIONS:STRING
    MEMORYCHECK_COMMAND:PATH
    CMAKE_SHARED_LINKER_FLAGS:STRING
    CMAKE_EXE_LINKER_FLAGS:STRING
    CMAKE_MODULE_LINKER_FLAGS:STRING
    SITE:STRING
    BUILDNAME:STRING
    ${PROJECT_NAME}_BUILD_DICOM_SUPPORT:BOOL
    SlicerExecutionModel_DIR:PATH
    SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_LIBRARY_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
    ${PRIMARY_PROJECT_NAME}_REQUIRES_VTK:BOOL
  ALL_PROJECTS
  )

if(${LOCAL_PROJECT_NAME}_USE_QT)
  mark_as_superbuild(
    VARS
      ${PRIMARY_PROJECT_NAME}_USE_QT:BOOL
      QT_QMAKE_EXECUTABLE:PATH
      QT_MOC_EXECUTABLE:PATH
      QT_UIC_EXECUTABLE:PATH
    ALL_PROJECTS
    )
endif()

if(USE_BRAINSCut)
  mark_as_superbuild(
    VARS
     OpenCV_DIR:PATH
    ALL_PROJECTS
  )
endif()

set(extProjName ${LOCAL_PROJECT_NAME})
set(proj        ${LOCAL_PROJECT_NAME})

#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------

set(ep_common_c_flags "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_INIT} ${ADDITIONAL_C_FLAGS}")
set(ep_common_cxx_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_INIT} ${ADDITIONAL_CXX_FLAGS}")

#-----------------------------------------------------------------------------
# Set CMake OSX variable to pass down the external project
#-----------------------------------------------------------------------------
set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
if(APPLE)
  list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
    -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
endif()

set(${LOCAL_PROJECT_NAME}_CLI_RUNTIME_DESTINATION  bin)
set(${LOCAL_PROJECT_NAME}_CLI_LIBRARY_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION  bin)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION  lib)

#-----------------------------------------------------------------------------
# Add external project CMake args
#-----------------------------------------------------------------------------
mark_as_superbuild(
  VARS
  BUILD_EXAMPLES:BOOL
  BUILD_TESTING:BOOL
  ITK_VERSION_MAJOR:STRING
  ITK_DIR:PATH
  VTK_DIR:PATH
  SlicerExecutionModel_DIR:PATH
  ${LOCAL_PROJECT_NAME}_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
  ${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
  ${LOCAL_PROJECT_NAME}_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
  ${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION:PATH
  ${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
  ${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION:PATH

  BOOST_INCLUDE_DIR:PATH
  BRAINS_DEBUG_IMAGE_WRITE:BOOL
  INSTALL_RUNTIME_DESTINATION:STRING
  INSTALL_LIBRARY_DESTINATION:STRING
  INSTALL_ARCHIVE_DESTINATION:STRING

  CMAKE_CXX_STANDARD:STRING

  TBB_ROOT:PATH
  TBB_BUILD_PREFIX:STRING
  TBB_BUILD_DIR:PATH
ALL_PROJECTS
)

#
# By default we want to build BRAINSTools stuff using the CMAKE_BUILD_TYPE of
# the top level build, but build the support libraries in Release.
# So make one list of parameters to pass to BRAINSTools when we build it and
# another for all the prerequisite libraries
#
# since we use a macro to build the list of arguments, it's easier to modify the
# list after it's built than try and conditionally change just the build type in the macro.


#------------------------------------------------------------------------------
# Calling this macro last will ensure all prior calls to 'mark_as_superbuild' are
# considered when updating the variable '${proj}_EP_ARGS' passed to the main project
# below.
ExternalProject_Include_Dependencies(${proj} DEPENDS_VAR ${PRIMARY_PROJECT_NAME}_DEPENDENCIES)

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
set(proj ${LOCAL_PROJECT_NAME})
ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  DEPENDS ${${LOCAL_PROJECT_NAME}_DEPENDENCIES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${LOCAL_PROJECT_NAME}-build
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    --no-warn-unused-cli    # HACK Only expected variables should be passed down.
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DBRAINSTools_MAX_TEST_LEVEL:STRING=${BRAINSTools_MAX_TEST_LEVEL}
    -DBRAINS_DEBUG_IMAGE_WRITE:BOOL=${BRAINS_DEBUG_IMAGE_WRITE}
    -DBUILD_STYLE_UTILS:BOOL=${BUILD_STYLE_UTILS}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${CMAKE_OSX_DEPLOYMENT_TARGET}
    -DCMAKE_OSX_SYSROOT:STRING=${CMAKE_OSX_SYSROOT}
    -DEXTERNAL_PROJECT_BUILD_TYPE:STRING=${EXTERNAL_PROJECT_BUILD_TYPE}
    -DFORCE_EXTERNAL_BUILDS:BOOL=${FORCE_EXTERNAL_BUILDS}
    -DITK_VERSION_MAJOR:STRING=${ITK_VERSION_MAJOR}
    -DSuperBuild_BRAINSTools_BUILD_DICOM_SUPPORT:BOOL=${SuperBuild_BRAINSTools_BUILD_DICOM_SUPPORT}
    -DSuperBuild_BRAINSTools_USE_CTKAPPLAUNCHER:BOOL=${SuperBuild_BRAINSTools_USE_CTKAPPLAUNCHER}
    -DSuperBuild_BRAINSTools_USE_GIT_PROTOCOL:BOOL=${SuperBuild_BRAINSTools_USE_GIT_PROTOCOL}
    -DUSE_SYSTEM_DCMTK:BOOL=${USE_SYSTEM_DCMTK}
    -DUSE_SYSTEM_ITK:BOOL=${USE_SYSTEM_ITK}
    -DUSE_SYSTEM_SlicerExecutionModel:BOOL=${USE_SYSTEM_SlicerExecutionModel}
    -DUSE_SYSTEM_VTK:BOOL=${USE_SYSTEM_VTK}
    -DVTK_GIT_REPOSITORY:STRING=${VTK_GIT_REPOSITORY}
    -D${LOCAL_PROJECT_NAME}_SUPERBUILD:BOOL=OFF
    -DANTs_SOURCE_DIR:PATH=${ANTs_SOURCE_DIR}
    -DANTs_LIBRARY_DIR:PATH=${ANTs_LIBRARY_DIR}
    -DBRAINSTools_LIBRARY_PATH:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
    -DTBB_ROOT:PATH=${TBB_ROOT}
    -DTBB_BUILD_PREFIX:STRING=${TBB_BUILD_PREFIX}
    -DTBB_BUILD_DIR:PATH=${TBB_BUILD_DIR}
    -DRTK_DIR:PATH=${RTK_DIR}
  INSTALL_COMMAND ""
  )

# This custom external project step forces the build and later
# steps to run whenever a top level build is done...
#
# BUILD_ALWAYS flag is available in CMake 3.1 that allows force build
# of external projects without this workaround. Remove this workaround
# and use the CMake flag instead, when BRAINSTools's required minimum CMake
# version will be at least 3.1.
#
if(CMAKE_CONFIGURATION_TYPES)
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${CMAKE_CFG_INTDIR}/${proj}-build")
else()
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build")
endif()
ExternalProject_Add_Step(${proj} forcebuild
  COMMAND ${CMAKE_COMMAND} -E remove ${BUILD_STAMP_FILE}
  COMMENT "Forcing build step for '${proj}'"
  DEPENDEES build
  ALWAYS 1
  )
