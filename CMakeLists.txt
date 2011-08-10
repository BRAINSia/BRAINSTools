project(BRAINSTools)
cmake_minimum_required(VERSION 2.8)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#-----------------------------------------------------------------------------
option(BRAINSTools_SUPERBUILD "Build Slicer and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(BRAINSTools_SUPERBUILD)
#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(BRAINSTools_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/BRAINSTools.cmake")
  return()
endif()

