#-----------------------------------------------------------------------------
enable_language(C)
enable_language(CXX)

#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
# Git protocole option
#-----------------------------------------------------------------------------
option(${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT ${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

find_package(Git REQUIRED)

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

#option(BRAINSTools_BUILD_WITH_ITKv4 "Build BRAINS against ITKv4" ON)
#mark_as_advanced(BRAINSTools_BUILD_WITH_ITKv4)
# TODO To explicitly support building against ITKv3, External_ITKv3.cmake should be added.
set(BRAINSTools_BUILD_WITH_ITKv4 ON CACHE BOOL "Build BRAINS against ITKv4")
  
option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)

# The following variable could be set at configuration time to customize where the 
# CLI executables and libraries are either built or installed.
#BRAINSTools_CLI_LIBRARY_OUTPUT_DIRECTORY
#BRAINSTools_CLI_ARCHIVE_OUTPUT_DIRECTORY
#BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY
#BRAINSTools_CLI_INSTALL_LIBRARY_DESTINATION
#BRAINSTools_CLI_INSTALL_ARCHIVE_DESTINATION
#BRAINSTools_CLI_INSTALL_RUNTIME_DESTINATION

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

set(ep_common_compiler_args
  -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_MODULE_LINKER_FLAGS}
  )

#------------------------------------------------------------------------------
# BRAINS dependency list
#------------------------------------------------------------------------------
if(BRAINSTools_BUILD_WITH_ITKv4)
  set(ITK_EXTERNAL_NAME "ITKv4")
else()
  set(ITK_EXTERNAL_NAME "ITKv3")
endif()

set(BRAINSTools_DEPENDENCIES VTK ${ITK_EXTERNAL_NAME} SlicerExecutionModel)

if(BUILD_STYLE_UTILS)
  list(APPEND BRAINSTools_DEPENDENCIES Cppcheck KWStyle Uncrustify)
endif()

if(USE_BRAINSABC) # OR USE_BRAINSCut)
  list(APPEND BRAINSTools_DEPENDENCIES ReferenceAtlas)
endif()

#-----------------------------------------------------------------------------
# Include external projects
#-----------------------------------------------------------------------------
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS)

# This variable will contain the list of CMake variable specific to each external project
# that should passed to ${CMAKE_PROJECT_NAME}.
# The item of this list should have the following form: <EP_VAR>:<TYPE>
# where '<EP_VAR>' is an external project variable and TYPE is either BOOL, STRING, PATH or FILEPATH.
# TODO Variable appended to this list will be automatically exported in BRAINSToolsConfig.cmake, 
# prefix 'BRAINSTools_' will be prepended if it applied.
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS)

SlicerMacroCheckExternalProjectDependency(BRAINSTools)

#-----------------------------------------------------------------------------
# Additionnal superbuild args
#-----------------------------------------------------------------------------

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS

  BRAINSTools_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
  BRAINSTools_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
  BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
  BRAINSTools_CLI_INSTALL_LIBRARY_DESTINATION:PATH
  BRAINSTools_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
  BRAINSTools_CLI_INSTALL_RUNTIME_DESTINATION:PATH

  BUILD_TESTING:BOOL
  BRAINSTools_BUILD_WITH_ITKv4:BOOL
  USE_GTRACT:BOOL
  USE_BRAINSFit:BOOL
  USE_BRAINSCommonLib:BOOL
  USE_BRAINSABC:BOOL
  USE_BRAINSCUT:BOOL
  USE_BRAINSMush:BOOL
  USE_BRAINSROIAuto:BOOL
  USE_BRAINSResample:BOOL
  USE_BRAINSConstellationDetector:BOOL
  USE_BRAINSDemonWarp:BOOL
  USE_BRAINSMultiModeSegment:BOOL
  USE_BRAINSInitializedControlPoints:BOOL
  USE_BRAINSTransformConvert:BOOL
  USE_ImageCalculator:BOOL
  )

#-----------------------------------------------------------------------------
# Expand superbuild external project args
#-----------------------------------------------------------------------------
SET(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES)
FOREACH(arg ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS})
  STRING(REPLACE ":" ";" varname_and_vartype ${arg})
  SET(target_info_list ${target_info_list})
  LIST(GET varname_and_vartype 0 _varname)
  LIST(GET varname_and_vartype 1 _vartype)
  LIST(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS -D${_varname}:${_vartype}=${${_varname}})
  LIST(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES ${_varname})
ENDFOREACH()
STRING(REPLACE ";" "^" ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES "${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES}")

# MESSAGE("CMake external project args:")
# FOREACH(arg ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})
#   MESSAGE("  ${arg}")
# ENDFOREACH()

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

#------------------------------------------------------------------------------
# Configure and build
#------------------------------------------------------------------------------
set(proj BRAINSTools)
ExternalProject_Add(${proj}
  DEPENDS ${BRAINSTools_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR BRAINSTools-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
    ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS}
    -DBRAINSTools_SUPERBUILD:BOOL=OFF
    -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
  INSTALL_COMMAND ""
  )

