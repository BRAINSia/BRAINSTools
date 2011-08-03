find_package(Git REQUIRED)
include(ExternalProject)
enable_language(C)
enable_language(CXX)

if(NOT SETIFEMPTY)
macro(SETIFEMPTY)
  set(KEY ${ARGV0})
  set(VALUE ${ARGV1})
  if(NOT ${KEY})
    set(${ARGV})
  endif(NOT ${KEY})
endmacro(SETIFEMPTY KEY VALUE)
endif(NOT SETIFEMPTY)
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/SuperBuild
  ${CMAKE_BINARY_DIR}/CMake
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake #  CMake directory
  ${CMAKE_CURRENT_SOURCE_DIR}/src/CMake # CMake directory
  ${CMAKE_MODULE_PATH}
  )

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
# Prerequisites
#------------------------------------------------------------------------------
#
# BRAINS4 Addition: install to the common library
# directory, so that all libs/include etc ends up
# in one common tree
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Where all the prerequisite libraries go" FORCE)
set(${CMAKE_PROJECT_NAME}_BUILD_TESTING ON CACHE BOOL "Turn on Testing for BRAINS")

# Compute -G arg for configuring external projects with the same CMake generator:
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()


#-------------------------------------------------------------------------
# augment compiler flags
#-------------------------------------------------------------------------
include(CompilerFlagSettings)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS}" )
else() # Release, or anything else
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS}" )
endif()

#------------------------------------------------------------------------------
# Conditionnaly include ExternalProject Target
#------------------------------------------------------------------------------

set(ep_common_args
  --no-warn-unused-cli
  -DMAKECOMMAND:STRING=${MAKECOMMAND}
  -DCMAKE_SKIP_RPATH:BOOL=ON
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
  -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
  -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=${BUILD_TESTING}
  -DCMAKE_GENERATOR:STRING=${CMAKE_GENERATOR}
  -DCMAKE_EXTRA_GENERATOR:STRING=${CMAKE_EXTRA_GENERATOR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  -DCMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH=${CMAKE_BUNDLE_OUTPUT_DIRECTORY}
  -DCTEST_NEW_FORMAT:BOOL=ON
  -DMEMORYCHECK_COMMAND_OPTIONS:STRING=${MEMORYCHECK_COMMAND_OPTIONS}
  -DMEMORYCHECK_COMMAND:PATH=${MEMORYCHECK_COMMAND}
  -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_MODULE_LINKER_FLAGS}
  -DSITE:STRING=${SITE}
  -DBUILDNAME:STRING=${BUILDNAME}
)


#------------------------------------------------------------------------------
# Determine if building stand-alone, or using external versions of ITKv4
# and SEM (i.e. for tight integration with Slicer)
#------------------------------------------------------------------------------
include(SlicerMacroEmptyExternalProject)
option(USE_SYSTEM_ITKv4             "Build using an externally defined version of ITKv4" OFF)
if(USE_SYSTEM_ITKv4)
  find_package(ITK 4 REQUIRED)
  include(${ITK_USE_FILE})
  SlicerMacroEmptyExternalProject("ITKv4" "")
else()
  include(External_ITKv4)
endif()

option(USE_SYSTEM_SEM               "Build using an externally defined version of SEM"  OFF)
if(USE_SYSTEM_SEM)
  find_package(SlicerExecutionModel NO_MODULE REQUIRED GenerateCLP)
  if(GenerateCLP_DIR)
    include(${GenerateCLP_USE_FILE})
  else(GenerateCLP_DIR)
    message(FATAL_ERROR "Can't build without GenerateCLP. Please set GenerateCLP_DIR")
  endif(GenerateCLP_DIR)
  SlicerMacroEmptyExternalProject("SlicerExecutionModel" "")
else()
  set(SlicerExecutionModel_DEPENDENCIES ITKv4)
  include(External_SlicerExecutionModel)
endif()

include(External_BRAINSTools)

set(BRAINSStandAlone_DEPENDENCIES ${BRAINSTools_PROVIDES})
#------------------------------------------------------------------------------
# Configure and build
#------------------------------------------------------------------------------
set(proj BRAINSStandAlone)
ExternalProject_Add(${proj}
  DEPENDS ${BRAINSStandAlone_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR BRAINSStandAlone-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${ep_common_args}
    -DBRAINSStandAlone_SUPERBUILD:BOOL=OFF
    -DADDITIONAL_CXX_FLAGS:STRING=${ADDITIONAL_CXX_FLAGS}
    -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
    # ITK
    -DITK_DIR:PATH=${ITK_DIR}
    # VTK
    -DVTK_DIR:PATH=${VTK_DIR}
    # SlicerExecutionModel_DIR
    -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
  INSTALL_COMMAND ""
  )
