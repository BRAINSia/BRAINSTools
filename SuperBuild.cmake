#
# Stuff that needs to be set if you're building with Slicer
if(NOT INTEGRATE_WITH_SLICER)
  set(Slicer_VTK_GIT_REPOSITORY
    "github.com/Slicer/VTK.git" CACHE STRING "repository from which to get VTK" FORCE)
  mark_as_advanced(Slicer_VTK_GIT_REPOSITORY)
  # By default, use a specific SHA1 associated with branch slicer-4.0-gamma on github.com/Slicer/VTK.git
  set(Slicer_VTK_GIT_TAG
    "9330805c45444eb5b740bf401aec50f4c32f3cab" CACHE STRING "VTK git tag to use" FORCE)
  mark_as_advanced(Slicer_VTK_GIT_TAG)
  option(Slicer_USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)

  set(git_protocol "git")
  if(NOT Slicer_USE_GIT_PROTOCOL)
    set(git_protocol "http")
  endif()

endif(NOT INTEGRATE_WITH_SLICER)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


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
option(USE_SYSTEM_ITK             "Build using an externally defined version of ITKv4" OFF)
if(USE_SYSTEM_ITK)
  find_package(ITK REQUIRED)
  include(${ITK_USE_FILE})
  message(STATUS "ITKv4 EmptyProject Set")
  SlicerMacroEmptyExternalProject("ITKv4" "")
else()
  include(External_ITKv4)
endif()

option(USE_SYSTEM_VTK             "Build using an externally defined version of VTK" OFF)
if(USE_SYSTEM_VTK)
  find_package(VTK  REQUIRED)
  include(${VTK_USE_FILE})
  SlicerMacroEmptyExternalProject("VTK" "")
  #link_libraries( vtkGraphics vtkImaging vtkIO vtkFiltering vtkCommon )
else()
  include(External_VTK)
endif()


option(USE_SYSTEM_SEM               "Build using an externally defined version of SEM"  OFF)
if(USE_SYSTEM_SEM)
  find_package(SlicerExecutionModel NO_MODULE REQUIRED GenerateCLP)
  if(GenerateCLP_DIR)
    include(${GenerateCLP_USE_FILE})
  else(GenerateCLP_DIR)
    message(FATAL_ERROR "Can't build without GenerateCLP. Please set GenerateCLP_DIR")
  endif(GenerateCLP_DIR)
  SlicerMacroEmptyExternalProject("SlicerExecutionModel" "ITKv4")
else()
  set(SlicerExecutionModel_DEPENDENCIES ITKv4)
  include(External_SlicerExecutionModel)
endif()

option(USE_GTRACT                      "Build GTRACT"                      ON)
option(USE_BRAINSCommonLib             "Build BRAINSCommonLib"             ON)
option(USE_BRAINSFit                   "Build BRAINSFit"                   ON)
option(USE_BRAINSABC                   "Build BRAINSABC"                   ON)
option(USE_BRAINSROIAuto               "Build BRAINSROIAuto"               ON)
option(USE_BRAINSDemonWarp             "Build BRAINSDemonWarp"             ON)
option(USE_BRAINSResample              "Build BRAINSResample"              ON)
option(USE_BRAINSConstellationDetector "Build BRAINSConstellationDetector" ON)
option(USE_BRAINSMush                  "Build BRAINSMush"                  ON)
option(USE_BRAINSMultiModeSegment      "Build BRAINSMultiModeSegment"      ON)
option(USE_BRAINSInitializedControlPoints      "Build BRAINSInitializedControlPoints"      ON)
option(USE_BRAINSTransformConvert       "Build BRAINSTransformConvert"     ON)
#option(USE_BRAINSCut                   "Build BRAINSCut"                   OFF)

set(BRAINSTools_DEPENDENCIES ITKv4 SlicerExecutionModel VTK)
if(USE_BRAINSABC) # OR USE_BRAINSCut)
  # Define the atlas subdirectory in one place
  set(${CMAKE_PROJECT_NAME}_RUNTIME_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/bin)
  include(External_ReferenceAtlas)
  list(APPEND  BRAINSTools_DEPENDENCIES ${ReferenceAtlas_DEPEND})
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
    ${ep_common_args}
    -DBRAINSTools_SUPERBUILD:BOOL=OFF
    -DADDITIONAL_CXX_FLAGS:STRING=${ADDITIONAL_CXX_FLAGS}
    -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
    #################### Propogate the SEM Environment
    ## -- This could be some other variable to indicate a slicer build
    -DSEM_PLUGINS_BIN_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/bin
    -DSEM_PLUGINS_LIB_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/lib
    -DSEM_INSTALL_PLUGINS_BIN_DIR:STRING=/bin
    -DSEM_INSTALL_PLUGINS_LIB_DIR:STRING=/lib
    #################### Propogate the SEM Environment
    # ITK
    -DITK_DIR:PATH=${ITK_DIR}
    # VTK
    -DVTK_DIR:PATH=${VTK_DIR}
    # SlicerExecutionModel_DIR
    -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
    -DINTEGRATE_WITH_SLICER:BOOL=${INTEGRATE_WITH_SLICER}
    -DSlicer_SOURCE_DIR:PATH=${Slicer_SOURCE_DIR}
    -DBUILD_TESTING:BOOL=ON
    -DUSE_GTRACT:BOOL=${USE_GTRACT}
    -DUSE_BRAINSFit:BOOL=${USE_BRAINSFit}
    -DUSE_BRAINSCommonLib:BOOL=${USE_BRAINSCommonLib}
    -DUSE_BRAINSABC:BOOL=${USE_BRAINSABC}
    -DUSE_BRAINSMush:BOOL=${USE_BRAINSMush}
    -DUSE_BRAINSROIAuto:BOOL=${USE_BRAINSROIAuto}
    -DUSE_BRAINSResample:BOOL=${USE_BRAINSResample}
    -DUSE_BRAINSConstellationDetector:BOOL=${USE_BRAINSConstellationDetector}
    -DUSE_BRAINSDemonWarp:BOOL=${USE_BRAINSDemonWarp}
    -DUSE_BRAINSMultiModeSegment:BOOL=${USE_BRAINSMultiModeSegment}
    -DUSE_BRAINSInitializedControlPoints:BOOL=${USE_BRAINSInitializedControlPoints}
    -DUSE_BRAINSTransformConvert:BOOL=${USE_BRAINSTransformConvert}
    -D${CMAKE_PROJECT_NAME}_USE_ITK4:BOOL=ON
    ${VTK_BUILD_FLAGS}
    ${OpenCV_BUILD_FLAGS}
    ${QT_BUILD_FLAGS}
  INSTALL_COMMAND ""
  )
