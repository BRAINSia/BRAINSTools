# This file contains content that is only used for not BRAINSTools superbuild
# configuration

#-----------------------------------------------------------------------------
include(SlicerMacroGetOperatingSystemArchitectureBitness)

#-----------------------------------------------------------------------------
# Where should the superbuild source files be downloaded to?
# By keeping this outside of the build tree, you can share one
# set of external source trees for multiple build trees
#-----------------------------------------------------------------------------
# set( SOURCE_DOWNLOAD_CACHE ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
set( SOURCE_DOWNLOAD_CACHE ${CMAKE_CURRENT_BINARY_DIR} ) #<-- Note same as default

#-----------------------------------------------------------------------------
# Git protocol option
#-----------------------------------------------------------------------------
option(${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "https")
if(NOT ${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL)
  set(git_protocol "https")
endif()

find_package(Git REQUIRED)

cmake_dependent_option(${CMAKE_PROJECT_NAME}_USE_CTKAPPLAUNCHER "CTKAppLauncher used with python" ON
  "NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_python" OFF)

# Compute -G arg for configuring external projects with the same CMake generator:
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()


set(cmakeversion_external_update LOG_UPDATE )
set(cmakeversion_external_update_value 1)

#-----------------------------------------------------------------------------
# Superbuild option(s)
#-----------------------------------------------------------------------------
option(BUILD_STYLE_UTILS "Build uncrustify, cppcheck, & KWStyle" OFF)
cmake_dependent_option(
  USE_SYSTEM_Uncrustify "Use system Uncrustify program" OFF
  "BUILD_STYLE_UTILS" OFF
  )
cmake_dependent_option(
  USE_SYSTEM_KWStyle "Use system KWStyle program" OFF
  "BUILD_STYLE_UTILS" OFF
  )
cmake_dependent_option(
  USE_SYSTEM_Cppcheck "Use system Cppcheck program" OFF
  "BUILD_STYLE_UTILS" OFF
  )

option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
option(USE_SYSTEM_zlib "build using the system version of zlib" OFF)
#option(USE_SYSTEM_DCMTK "Build using an externally defined version of DCMTK" OFF)

#------------------------------------------------------------------------------
# ${LOCAL_PROJECT_NAME} dependency list
#------------------------------------------------------------------------------

#if(${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT)
#  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES DCMTK)
#endif()
list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES ITKv5)
list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES SlicerExecutionModel)
# Teem is not needed list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES teem)
#list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES Boost)

if(BUILD_STYLE_UTILS)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES Cppcheck KWStyle Uncrustify)
endif()


if(USE_BRAINSABC)
  if(${SUPERBUILD_TOPLEVEL_PROJECT}_REQUIRES_TBB)
    list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES TBB)
  else()
    message(FATAL_ERROR "USE_BRAINSABC Requires TBB")
  endif()
endif()

if(${LOCAL_PROJECT_NAME}_REQUIRES_VTK)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES VTK)
endif()

if(USE_ANTS)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES ANTs)
endif()

## RTK now part of ITK remote module if(USE_BRAINSSuperResolution)
## RTK now part of ITK remote module   list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES RTK)
## RTK now part of ITK remote module endif()

set(extProjName ${LOCAL_PROJECT_NAME})
set(proj        ${LOCAL_PROJECT_NAME})

#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------

set(${LOCAL_PROJECT_NAME}_CLI_RUNTIME_DESTINATION  bin)
set(${LOCAL_PROJECT_NAME}_CLI_LIBRARY_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION  bin)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION  lib)

set(${LOCAL_PROJECT_NAME}_INSTALL_LIB_DIR ${${LOCAL_PROJECT_NAME}_CLI_LIBRARY_DESTINATION} )

#-----------------------------------------------------------------------------
# Common external projects CMake variables
#-----------------------------------------------------------------------------

if(${LOCAL_PROJECT_NAME}_USE_QT)
  mark_as_superbuild(
    VARS
      ${LOCAL_PROJECT_NAME}_USE_QT:BOOL
      QT_QMAKE_EXECUTABLE:PATH
      QT_MOC_EXECUTABLE:PATH
      QT_UIC_EXECUTABLE:PATH
    PROJECTS VTK ${LOCAL_PROJECT_NAME}
    )
endif()


#
# By default we want to build BRAINSTools stuff using the CMAKE_BUILD_TYPE of
# the top level build, but build the support libraries in Release.
# So make one list of parameters to pass to BRAINSTools when we build it and
# another for all the prerequisite libraries
#
# since we use a macro to build the list of arguments, it's easier to modify the
# list after it's built than try and conditionally change just the build type in the macro.
set(BRAINSTools_LIBRARY_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# These variables are written to the cache used to initialize this external project
mark_as_superbuild( # Local project
  VARS
  #    ${LOCAL_PROJECT_NAME}_SUPERBUILD:BOOL
    BUILD_OPTIMIZED:BOOL
    BUILD_COVERAGE:BOOL
    BUILD_STYLE_UTILS:BOOL

    BRAINSTools_ExternalData_DATA_MANAGEMENT_TARGET:STRING
    BRAINSTools_LIBRARY_PATH:PATH
    BRAINSTools_MAX_TEST_LEVEL:STRING
    BRAINS_DEBUG_IMAGE_WRITE:BOOL
    BRAINSTools_USE_CTKAPPLAUNCHER:BOOL
    BRAINSTools_USE_GIT_PROTOCOL:BOOL

    #BOOST_INCLUDE_DIR:PATH
    ${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT:BOOL
    #    USE_SYSTEM_DCMTK:BOOL
    USE_SYSTEM_ITK:BOOL
    ITK_DIR:PATH

    ${LOCAL_PROJECT_NAME}_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
    ${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
    ${LOCAL_PROJECT_NAME}_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION:PATH
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION:PATH

    USE_SYSTEM_SlicerExecutionModel:BOOL
    SlicerExecutionModel_DIR:PATH
    SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_LIBRARY_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
  PROJECTS ${LOCAL_PROJECT_NAME}
)
if(${LOCAL_PROJECT_NAME}_REQUIRES_VTK)
  mark_as_superbuild( # Local project
    VARS
      ${LOCAL_PROJECT_NAME}_REQUIRES_VTK:BOOL
      USE_SYSTEM_VTK:BOOL
      VTK_DIR:PATH
      VTK_GIT_REPOSITORY:STRING
    PROJECTS ${LOCAL_PROJECT_NAME}
  )
endif()


#------------------------------------------------------------------------------
# Calling this macro last will ensure all prior calls to 'mark_as_superbuild' are
# considered when updating the variable '${LOCAL_PROJECT_NAME}_EP_ARGS' passed to the main project
# below.
ExternalProject_Include_Dependencies( ${LOCAL_PROJECT_NAME}
   PROJECT_VAR proj
   EP_ARGS_VAR MYBRAINSTools_EP_ARGS ## Get mark_as_superbuild items for writing to cache file BRAINSTools-prefix/tmp/BRAINSTools-cache-Debug.cmake
   DEPENDS_VAR ${LOCAL_PROJECT_NAME}_DEPENDENCIES
)
#message(FATAL_ERROR "MYBRAINSTools_EP_ARGS:${MYBRAINSTools_EP_ARGS}:")
#------------------------------------------------------------------------------
# Write values to a file for demonstrating the config options
#set(WRITE_BRAINSTOOLS_ARGS "${MYBRAINSTools_EP_ARGS}")
#separate_arguments(WRITE_BRAINSTOOLS_ARGS)
#list(REMOVE_DUPLICATES WRITE_BRAINSTOOLS_ARGS)
#string(REPLACE ";" " " WRITE_BRAINSTOOLS_ARGS "${WRITE_BRAINSTOOLS_ARGS}")
#string(REPLACE "CMAKE_CACHE_ARGS" "" WRITE_BRAINSTOOLS_ARGS "${WRITE_BRAINSTOOLS_ARGS}")
#string(REPLACE "LIST_SEPARATOR.*" "" WRITE_BRAINSTOOLS_ARGS "${WRITE_BRAINSTOOLS_ARGS}")

#set( WRITE_BRAINSTOOLS_ARGS " cmake ${WRITE_BRAINSTOOLS_ARGS}")
#file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/BRAINSToolsArgs.sh ${WRITE_BRAINSTOOLS_ARGS})
#message(FATAL_ERROR "${cmd_string}")
#message(FATAL_ERROR "\n${WRITE_BRAINSTOOLS_ARGS}\n")
# message(FATAL_ERROR "${MYBRAINSTools_EP_ARGS}")

#------------------------------------------------------------------------------
# Configure and build ${PROJECT_NAME}
#------------------------------------------------------------------------------
ExternalProject_Add(${LOCAL_PROJECT_NAME}
  DEPENDS ${${LOCAL_PROJECT_NAME}_DEPENDENCIES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        BINARY_DIR ${LOCAL_PROJECT_NAME}-${CMAKE_BUILD_TYPE}-EP${EXTERNAL_PROJECT_BUILD_TYPE}-build
        DOWNLOAD_COMMAND ""
        UPDATE_COMMAND ""
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS
        --no-warn-unused-cli    # HACK Only expected variables should be passed down.

        CMAKE_CACHE_ARGS
        -D${LOCAL_PROJECT_NAME}_SUPERBUILD:BOOL=OFF  # This must be here
        -DBUILD_TESTING:BOOL=${${LOCAL_PROJECT_NAME}_BUILD_TESTING}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DBUILD_EXAMPLES:BOOL=${${LOCAL_PROJECT_NAME}_BUILD_TESTING}
        -DCMAKE_INSTALL_PREFIX_SET:BOOL=TRUE
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
        ${MYBRAINSTools_EP_ARGS}  # This appends the "mark_as_superbuild" variables to the written out cache for clion use
        #INSTALL_COMMAND ""
        LOG_CONFIGURE ON
        LOG_BUILD ON
        )

if(CMAKE_CONFIGURATION_TYPES)
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-prefix/src/${LOCAL_PROJECT_NAME}-stamp/${CMAKE_CFG_INTDIR}/${LOCAL_PROJECT_NAME}-build")
else()
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-prefix/src/${LOCAL_PROJECT_NAME}-stamp/${LOCAL_PROJECT_NAME}-build")
endif()
