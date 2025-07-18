cmake_minimum_required(VERSION 3.16.3...3.20.3)
cmake_policy(VERSION 3.16.3...3.20.3)

#-----------------------------------------------------------------------------
# Enable C++17
#-----------------------------------------------------------------------------
#####
##  Set the default target properties for ITK
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17) # Supported values are 17, 20, and 23.
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
set(VALID_CXX_STANDARDS "17" "20" "23")
if(NOT CMAKE_CXX_STANDARD IN_LIST VALID_CXX_STANDARDS )
   MESSAGE(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD} not in know standards list\n ${VALID_CXX_STANDARDS}  for BRAINSTools version 5 and greater.")
endif()

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_LIST_DIR}/CMake)

#-----------------------------------------------------------------------------
if(APPLE)
  # Note: By setting CMAKE_OSX_* variables before any enable_language() or project() calls,
  #       we ensure that the bitness, and C++ standard library will be properly detected.
  include(InitializeOSXVariables)
  set(required_deployment_target "11.0")
  if("x${CMAKE_OSX_DEPLOYMENT_TARGET}x" STREQUAL "xx")
    set(CMAKE_OSX_DEPLOYMENT_TARGET ${required_deployment_target})
  endif()
  if(CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS ${required_deployment_target})
    message(FATAL_ERROR "CMAKE_OSX_DEPLOYMENT_TARGET ${CMAKE_OSX_DEPLOYMENT_TARGET} must be ${required_deployment_target} or greater.")
  endif()
endif()

#-----------------------------------------------------------------------------
set(LOCAL_PROJECT_NAME BRAINSTools) # <-- the primary product endpoint (Often matches the superbuild name)
project(${LOCAL_PROJECT_NAME}
   LANGUAGES  C CXX
   DESCRIPTION "BRAINSTools see github at https://github.com/BRAINSia/BRAINSTools" #<-- cmake version 3.9 or greater
   VERSION 5.8.0
)
# VERSION 5.8.0 Update to fix bug in computing ComputeRobustMinMaxMean
# VERSION 5.7.0 Update ITK, VTK, ANTs, improve landmarksConstellationAligner with bug fix for finding RP point
# VERSION 5.6.0 Update ITK, VTK, ANTs, Remove DCMTK external, and improve landmarksConstellationAligner
# VERSION 5.5.0 The version that requires C++14 by default
# VERSION 5.4.0 The version used to create the final PREDICTHD_BIDS_DEFACE dataset.
# VERSION 5.3.2 Many tools placed in ARCHIVE

#-- The project() command stores the version number and its components in variables
#-- PROJECT_VERSION, <PROJECT-NAME>_VERSION
#-- PROJECT_VERSION_MAJOR, <PROJECT-NAME>_VERSION_MAJOR
#-- PROJECT_VERSION_MINOR, <PROJECT-NAME>_VERSION_MINOR
#-- PROJECT_VERSION_PATCH, <PROJECT-NAME>_VERSION_PATCH
#-- PROJECT_VERSION_TWEAK, <PROJECT-NAME>_VERSION_TWEAK
#-----------------------------------------------------------------------------
set(SUPERBUILD_TOPLEVEL_PROJECT ${LOCAL_PROJECT_NAME})
set(EXTERNAL_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild)
include(ExternalProject)
include(ExternalProjectDependency) #<-- Must be after SUPERBUILD_TOPLEVEL_PROJECT and EXTERNAL_PROJECT_DIR
include(ExternalProjectGenerateProjectDescription)
if( IS_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/.git")
  #-----------------------------------------------------------------------------
  # Version strings need to be set outside of project so that
  # we can get version information from Version.cmake (which rquires
  # that the project() command was already run.
  #-----------------------------------------------------------------------------
  # Version information
  include(Version.cmake) ## Do introspection for fine grained versioning, and tag matching
endif()
#-----------------------------------------------------------------------------
#-------------------------------------------------------------------------
if(NOT DEFINED BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET)
  set(BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET "BRAINSToolsFetchData")
endif()

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
option(${LOCAL_PROJECT_NAME}_PKGBUILD "Build ${LOCAL_PROJECT_NAME} SUPPORT_FILES ONLY (i.e. do not build BRAINSTools)." ON)
mark_as_advanced(${LOCAL_PROJECT_NAME}_PKGBUILD)

#-----------------------------------------------------------------------------
# Common build features for both the superbuild and the main build
#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake) #<-- All feature options for top superbuild and inner product build


#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake") #<-- Harness for managing top superbuild, and finally delegate inner product build
  return()
elseif(${LOCAL_PROJECT_NAME}_PKGBUILD )
#-----------------------------------------------------------------------------
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake") #<-- inner product build (a separate build environment from top level build)
  return()
else()
  message(STATUS "BRAINSTools package not built because ${LOCAL_PROJECT_NAME}_PKGBUILD=${${LOCAL_PROJECT_NAME}_PKGBUILD}")
  return()
endif()

if(CMAKE_EXPORT_COMPILE_COMMANDS)
  message(STATUS "Exporting compiler flags") #silence cmake warnings when building jason db of compile options
endif()
message(FATAL_ERROR "You should never reach this point !")
