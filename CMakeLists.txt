cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

## NOTE THERE SHOULD BE NO PROJECT STATEMENT HERE!
## This file acts as a simple switch to initiate
## two completely independant CMake build environments.

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#                   Phase I:  BRAINSTools_SUPERBUILD is set to ON, and the
#                             supporting packages defined in "SuperBuild.cmake"
#                             are built.  The last package in "SuperBuild.cmake"
#                             to be built is a recursive call to this
#                             file with BRAINSTools_SUPERBUILD explicitly
#                             set to "OFF" to initiate Phase II
#
#                   Phase II: Build the BRAINSTools, referencing the support
#                             packages built in Phase I.
#-----------------------------------------------------------------------------
option(BRAINSTools_SUPERBUILD "Build BRAINSTools and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(BRAINSTools_SUPERBUILD)

option(BRAINSTools_USE_QT "Find and use Qt with VTK to build GUI Tools" OFF)
mark_as_advanced(BRAINSTools_USE_QT)

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(BRAINSTools_SUPERBUILD)
  project(SuperBuild_BRAINSTools)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  project(BRAINSTools)
  include("${CMAKE_CURRENT_SOURCE_DIR}/BRAINSTools.cmake")
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")

