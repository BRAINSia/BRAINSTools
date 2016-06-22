if(DEFINED CMAKE_CXX_STANDARD AND (NOT ${CMAKE_CXX_STANDARD} EQUAL "98" ))
  cmake_minimum_required(VERSION 3.2.0) ## 3.2.0 needed for AppleClang CXX11 finding
else()
  cmake_minimum_required(VERSION 2.8.9) ## Needed for slicer compliance
endif()

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 OLD) #OLD project() behavior for setting VERSION variable
endif()
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD) #Prefer files from the CMake module directory when including from there
endif()
foreach(p
  CMP0054 # CMake 3.1 Only interpret ``if()`` arguments as variables or keywords when unquoted.
  CMP0042 # ``MACOSX_RPATH`` is enabled by default.
  CMP0025 # AppleClang vs. regular Clang
  CMP0063 # Honor visibility properties for all target types
  )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
include(ExternalProject)
include(ExternalProjectDependency)
#-----------------------------------------------------------------------------
if(APPLE)
  # Note: By setting CMAKE_OSX_* variables before any enable_language() or project() calls,
  #       we ensure that the bitness, and C++ standard library will be properly detected.
  include(BlockSetCMakeOSXVariables)
  mark_as_superbuild(
    VARS CMAKE_OSX_ARCHITECTURES:STRING CMAKE_OSX_SYSROOT:PATH CMAKE_OSX_DEPLOYMENT_TARGET:STRING
    ALL_PROJECTS
    )
endif()

set(LOCAL_PROJECT_NAME BRAINSTools)
set(PRIMARY_PROJECT_NAME ${LOCAL_PROJECT_NAME})

## NOTE THERE SHOULD BE NO PROJECT STATEMENT HERE!
## This file acts as a simple switch to initiate
## two completely independant CMake build environments.

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#                   Phase I:  ${LOCAL_PROJECT_NAME}_SUPERBUILD is set to ON, and the
#                             supporting packages defined in "SuperBuild.cmake"
#                             are built.  The last package in "SuperBuild.cmake"
#                             to be built is a recursive call to this
#                             file with ${LOCAL_PROJECT_NAME}_SUPERBUILD explicitly
#                             set to "OFF" to initiate Phase II
#
#                   Phase II: Build the ${LOCAL_PROJECT_NAME}, referencing the support
#                             packages built in Phase I.
#-----------------------------------------------------------------------------
option(${LOCAL_PROJECT_NAME}_SUPERBUILD "Build ${LOCAL_PROJECT_NAME} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${LOCAL_PROJECT_NAME}_SUPERBUILD)

if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  project(${PRIMARY_PROJECT_NAME} C CXX)
else()
  project(${LOCAL_PROJECT_NAME} C CXX)
endif()

## For compatibility with Slicer
if(NOT CMAKE_CXX_STANDARD)
   set(CMAKE_CXX_STANDARD 98) # Supported values are ``98``, ``11`` and ``14``.
endif()

enable_language(C)
enable_language(CXX)

#-----------------------------------------------------------------------------
# Common build features for both the superbuild and the main build
#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)
include(itkCheckSourceTree)

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  project(SuperBuild_${LOCAL_PROJECT_NAME} C CXX)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")

