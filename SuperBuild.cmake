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
# CTestCustom
#-----------------------------------------------------------------------------
if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
  configure_file(
    CMake/CTestCustom.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake
    @ONLY)
endif()

#-----------------------------------------------------------------------------
# Git protocol option
#-----------------------------------------------------------------------------
option(${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
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

#-----------------------------------------------------------------------------
# Set a default external project build type if none was specified
set(EXTERNAL_PROJECT_BUILD_TYPE "Release" CACHE STRING "Default build type for support libraries")
set_property(CACHE EXTERNAL_PROJECT_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
option(USE_SYSTEM_zlib "build using the system version of zlib" OFF)
option(USE_SYSTEM_DCMTK "Build using an externally defined version of DCMTK" OFF)

#------------------------------------------------------------------------------
# ${LOCAL_PROJECT_NAME} dependency list
#------------------------------------------------------------------------------

if(${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT)
list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES DCMTK)
endif()
list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES ITKv5)
list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES SlicerExecutionModel)
list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES teem)
#list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES Boost)

if(BUILD_STYLE_UTILS)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES Cppcheck KWStyle Uncrustify)
endif()

if(USE_BRAINSABC OR USE_BRAINSCut)
  #list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES qhull)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES TBB)
endif()

if(${LOCAL_PROJECT_NAME}_REQUIRES_VTK)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES VTK)
endif()

if(USE_BRAINSCut )
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES OpenCV)
endif()

if(USE_ANTS)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES ANTs)
endif()

if(USE_BRAINSSuperResolution)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES RTK)
endif()

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

set(extProjName ${LOCAL_PROJECT_NAME})
set(proj        ${LOCAL_PROJECT_NAME})

#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------

set(ep_common_c_flags "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_INIT} ${ADDITIONAL_C_FLAGS} ${BRAINSTools_C_OPTIMIZATION_FLAGS} ${BRAINSTools_C_WARNING_FLAGS}")
set(ep_common_cxx_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_INIT} ${ADDITIONAL_CXX_FLAGS} ${BRAINSTools_CXX_OPTIMIZATION_FLAGS} ${BRAINSTools_CXX_WARNING_FLAGS}")

set(${LOCAL_PROJECT_NAME}_CLI_RUNTIME_DESTINATION  bin)
set(${LOCAL_PROJECT_NAME}_CLI_LIBRARY_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION  bin)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION  lib)

set(${LOCAL_PROJECT_NAME}_INSTALL_LIB_DIR ${${LOCAL_PROJECT_NAME}_CLI_LIBRARY_DESTINATION} )

#-----------------------------------------------------------------------------
# Add external project CMake args
#-----------------------------------------------------------------------------
mark_as_superbuild(
  VARS
    BUILD_EXAMPLES:BOOL
    BUILD_TESTING:BOOL
    BUILD_SHARED_LIBS:BOOL

    MAKECOMMAND:STRING

    INSTALL_RUNTIME_DESTINATION:STRING
    INSTALL_LIBRARY_DESTINATION:STRING
    INSTALL_ARCHIVE_DESTINATION:STRING

    SITE:STRING
    BUILDNAME:STRING
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
set(BRAINSTools_LIBRARY_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
mark_as_superbuild(
  VARS
    USE_SYSTEM_SlicerExecutionModel:BOOL
#    SlicerExecutionModel_DIR:PATH
#    ITK_DIR:PATH
    VTK_DIR:PATH
    #BOOST_INCLUDE_DIR:PATH

    BRAINSTools_LIBRARY_PATH:PATH
    BRAINSTools_MAX_TEST_LEVEL:STRING

    ${LOCAL_PROJECT_NAME}_REQUIRES_VTK:BOOL
    BUILD_STYLE_UTILS:BOOL
    ${SUPERBUILD_TOPLEVEL_PROJECT}_BUILD_DICOM_SUPPORT:BOOL
    BRAINS_DEBUG_IMAGE_WRITE:BOOL

    BRAINSTools_USE_CTKAPPLAUNCHER:BOOL
    BRAINSTools_USE_GIT_PROTOCOL:BOOL
    EXTERNAL_PROJECT_BUILD_TYPE:STRING

    USE_SYSTEM_DCMTK:BOOL
    USE_SYSTEM_ITK:BOOL
    USE_SYSTEM_VTK:BOOL
    VTK_GIT_REPOSITORY:STRING

    ${LOCAL_PROJECT_NAME}_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
    ${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
    ${LOCAL_PROJECT_NAME}_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION:PATH
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION:PATH

    SlicerExecutionModel_DIR:PATH
    SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_LIBRARY_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
  PROJECTS ${LOCAL_PROJECT_NAME}
)


#------------------------------------------------------------------------------
# Calling this macro last will ensure all prior calls to 'mark_as_superbuild' are
# considered when updating the variable '${LOCAL_PROJECT_NAME}_EP_ARGS' passed to the main project
# below.
ExternalProject_Include_Dependencies( ${LOCAL_PROJECT_NAME}
   PROJECT_VAR proj
   EP_ARGS_VAR MYBRAINSTools_EP_ARGS
   DEPENDS_VAR ${LOCAL_PROJECT_NAME}_DEPENDENCIES
)

#------------------------------------------------------------------------------
# Write values to a file for demonstrating the config options
set(WRITE_BRAINSTOOLS_ARGS "${MYBRAINSTools_EP_ARGS}")
separate_arguments(WRITE_BRAINSTOOLS_ARGS)
list(REMOVE_DUPLICATES WRITE_BRAINSTOOLS_ARGS)
string(REPLACE ";" " " WRITE_BRAINSTOOLS_ARGS "${WRITE_BRAINSTOOLS_ARGS}")
string(REPLACE "CMAKE_CACHE_ARGS" "" WRITE_BRAINSTOOLS_ARGS "${WRITE_BRAINSTOOLS_ARGS}")
string(REPLACE "LIST_SEPARATOR.*" "" WRITE_BRAINSTOOLS_ARGS "${WRITE_BRAINSTOOLS_ARGS}")

set( WRITE_BRAINSTOOLS_ARGS " cmake -DBRAINSTools_SUPERBUILD:BOOL=OFF ${WRITE_BRAINSTOOLS_ARGS}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/BRAINSToolsArgs.sh ${WRITE_BRAINSTOOLS_ARGS})
#message(FATAL_ERROR "${cmd_string}")
#message(FATAL_ERROR "\n${WRITE_BRAINSTOOLS_ARGS}\n")

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
  ${MYBRAINSTools_EP_ARGS}  # All superbuild options should be passed by mark_as_superbuild
  CMAKE_CACHE_ARGS
    -D${LOCAL_PROJECT_NAME}_SUPERBUILD:BOOL=OFF #<-- Critical override
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DCMAKE_INCLUDE_DIRECTORIES_BEFORE:BOOL=OFF
  INSTALL_COMMAND ""
)

if(CMAKE_CONFIGURATION_TYPES)
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-prefix/src/${LOCAL_PROJECT_NAME}-stamp/${CMAKE_CFG_INTDIR}/${LOCAL_PROJECT_NAME}-build")
else()
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${LOCAL_PROJECT_NAME}-prefix/src/${LOCAL_PROJECT_NAME}-stamp/${LOCAL_PROJECT_NAME}-build")
endif()
