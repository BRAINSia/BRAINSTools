cmake_minimum_required(VERSION 3.9.5 FATAL_ERROR)

#-----------------------------------------------------------------------------
# Enable C++11
#-----------------------------------------------------------------------------
#####
##  Set the default target properties for ITK
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11) # Supported values are ``11``, ``14``, and ``17``.
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
set(VALID_CXX_STANDARDS "11" "14" "17")
if(NOT CMAKE_CXX_STANDARD IN_LIST VALID_CXX_STANDARDS )
   MESSAGE(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD} not in know standards list\n ${VALID_CXX_STANDARDS}  for BRAINSTools version 5 and greater.")
endif()

foreach(p
    ## Only policies introduced after the cmake_minimum_required
    ## version need to explicitly be set to NEW.
    CMP0070 #3.10.0 Define ``file(GENERATE)`` behavior for relative paths.
    CMP0071 #3.10.0 Let ``AUTOMOC`` and ``AUTOUIC`` process ``GENERATED`` files.
    CMP0067
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)

set(LOCAL_PROJECT_NAME BRAINSTools)
set(PRIMARY_PROJECT_NAME ${LOCAL_PROJECT_NAME})

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
  set(PROJECT_NAME SuperBuild_${LOCAL_PROJECT_NAME} )
  project(${PROJECT_NAME} LANGUAGES C CXX)
else()
  set(PROJECT_NAME ${LOCAL_PROJECT_NAME} )
  # Version information
  include(Version.cmake)
  set(_project_version "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MAJOR}")
  project(${PROJECT_NAME} LANGUAGES  C CXX VERSION ${_project_version})
endif()

#-----------------------------------------------------------------------------
# Common build features for both the superbuild and the main build
#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
#-----------------------------------------------------------------------------
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
  return()
endif()

if(CMAKE_EXPORT_COMPILE_COMMANDS)
  message(STATUS "Exporting compiler flags") #silence cmake warnings when building jason db of compile options
endif()
message(FATAL_ERROR "You should never reach this point !")
