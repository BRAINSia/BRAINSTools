
#-----------------------------------------------------------------------------
if(INTEGRATE_WITH_SLICER)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY        "${Slicer_PLUGINS_BIN_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY        "${Slicer_PLUGINS_LIB_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY        "${Slicer_PLUGINS_LIB_DIR}")
set(CMAKE_INSTALL_RUNTIME_DESTINATION     "${Slicer_INSTALL_PLUGINS_BIN_DIR}")
set(CMAKE_INSTALL_LIBRARY_DESTINATION     "${Slicer_INSTALL_PLUGINS_LIB_DIR}")
set(CMAKE_INSTALL_ARCHIVE_DESTINATION     "${Slicer_INSTALL_PLUGINS_LIB_DIR}")
else()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY        "${SEM_PLUGINS_BIN_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY        "${SEM_PLUGINS_LIB_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY        "${SEM_PLUGINS_LIB_DIR}")
set(CMAKE_INSTALL_RUNTIME_DESTINATION     "${SEM_INSTALL_PLUGINS_BIN_DIR}")
set(CMAKE_INSTALL_LIBRARY_DESTINATION     "${SEM_INSTALL_PLUGINS_LIB_DIR}")
set(CMAKE_INSTALL_ARCHIVE_DESTINATION     "${SEM_INSTALL_PLUGINS_LIB_DIR}")
set(CMAKE_BUNDLE_OUTPUT_DIRECTORY         "${SEM_INSTALL_PLUGINS_BIN_DIR}")
endif()
#message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY        ${SEM_PLUGINS_BIN_DIR}")
#message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY        ${SEM_PLUGINS_LIB_DIR}")
#message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY        ${SEM_PLUGINS_LIB_DIR}")
#message(STATUS "CMAKE_INSTALL_RUNTIME_DESTINATION     ${SEM_INSTALL_PLUGINS_BIN_DIR}")
#message(STATUS "CMAKE_INSTALL_LIBRARY_DESTINATION     ${SEM_INSTALL_PLUGINS_LIB_DIR}")
#message(STATUS "CMAKE_INSTALL_ARCHIVE_DESTINATION     ${SEM_INSTALL_PLUGINS_LIB_DIR}")
#message(STATUS "CMAKE_BUNDLE_OUTPUT_DIRECTORY         ${SEM_INSTALL_PLUGINS_BIN_DIR}")

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel NO_MODULE REQUIRED GenerateCLP)
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)
include(${SlicerExecutionModel_USE_FILE})
if(GenerateCLP_DIR)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_DIR)
  message(FATAL_ERROR "Can't build without GenerateCLP. Please set GenerateCLP_DIR")
endif(GenerateCLP_DIR)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

enable_testing()
include(CTest)

###
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

# Add needed flag for gnu on linux like enviroments to build static common libs
# suitable for linking with shared object libs.
if(NOT APPLE AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
  if(NOT "${CMAKE_C_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  endif()
endif()

##-----------------------------------------------------------------------
## Setup locaitons to find externally maintained test data.
##-----------------------------------------------------------------------
list(APPEND ExternalData_URL_TEMPLATES
  # Local data store populated by the ITK pre-commit hook
  "file:///${CMAKE_SOURCE_DIR}/.ExternalData/%(algo)/%(hash)"
  # Data published by Iowa Psychiatry web interface
  "http://www.psychiatry.uiowa.edu/users/brainstestdata/ctestdata/%(algo)/%(hash)"
  # Data published by MIDAS
  "http://midas.kitware.com/api/rest/midas.bitstream.by.hash?hash=%(hash)&algorithm=%(algo)"
  # Data published by developers using git-gerrit-push.
  "http://www.itk.org/files/ExternalData/%(algo)/%(hash)"
)

# Tell ExternalData commands to transform raw files to content links.
# TODO: Condition this feature on presence of our pre-commit hook.
set(ExternalData_LINK_CONTENT MD5)

set(TestData_DIR ${CMAKE_CURRENT_SOURCE_DIR}/TestData)


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## NOTE:  BRAINSCommonLib is REQUIRED.
# BRAINSCommonLib
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
set(BRAINSCommonLib_BUILDSCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts)

include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/PreventInSourceBuilds.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/CMakeBuildMacros.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/CMakeBRAINS3BuildMacros.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/IJMacros.txt)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/ExternalData.cmake)

add_subdirectory(BRAINSCommonLib)
set(BRAINSCommonLib_DIR ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib)

#find_package(BRAINSCommonLib NO_MODULE REQUIRED)
#include(${BRAINSCommonLib_USE_FILE})

#-----------------------------------------------------------------------------
# BRAINSFit
#-----------------------------------------------------------------------------
if(USE_BRAINSFit)
  add_subdirectory(BRAINSFit)
endif()

#-----------------------------------------------------------------------------
# BRAINSResample
#-----------------------------------------------------------------------------
if(USE_BRAINSResample)
  add_subdirectory(BRAINSResample)
endif()

#-----------------------------------------------------------------------------
# BRAINSDemonWarp
#-----------------------------------------------------------------------------
if(USE_BRAINSDemonWarp)
  add_subdirectory(BRAINSDemonWarp)
endif()

#-----------------------------------------------------------------------------
# BRAINSROIAuto
#-----------------------------------------------------------------------------
if(USE_BRAINSROIAuto)
  add_subdirectory(BRAINSROIAuto)
endif()

#-----------------------------------------------------------------------------
# GTRACT
#-----------------------------------------------------------------------------
if(USE_GTRACT)
  add_subdirectory(GTRACT)
endif()


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
if(ITK_VERSION_MAJOR GREATER 3)  ## Tools that only work with ITKv4
#-----------------------------------------------------------------------------
# BRAINSConstellationDetector
#-----------------------------------------------------------------------------
if(USE_BRAINSConstellationDetector)
  add_subdirectory(BRAINSConstellationDetector)
endif()


# Define the atlas subdirectory in one place
#set(${CMAKE_PROJECT_NAME}_RUNTIME_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/bin)

#-----------------------------------------------------------------------------
# BRAINSABC
#-----------------------------------------------------------------------------
if(USE_BRAINSABC)
  add_subdirectory(BRAINSABC)
endif()

#-----------------------------------------------------------------------------
# BRAINSCut
#-----------------------------------------------------------------------------
if(USE_BRAINSCut)
  add_subdirectory(BRAINSCut)
endif()

#-----------------------------------------------------------------------------
# BRAINSMush
#-----------------------------------------------------------------------------
if(USE_BRAINSMush)
  add_subdirectory(BRAINSMush)
endif()

#-----------------------------------------------------------------------------
# BRAINSMultiModeSegment
#-----------------------------------------------------------------------------
if(USE_BRAINSMultiModeSegment)
  add_subdirectory(BRAINSMultiModeSegment)
endif()

#-----------------------------------------------------------------------------
# BRAINSInitializedControlPoints
#-----------------------------------------------------------------------------
if(USE_BRAINSInitializedControlPoints)
  add_subdirectory(BRAINSInitializedControlPoints)
endif()

#-----------------------------------------------------------------------------
# BRAINSTransformConvert
#-----------------------------------------------------------------------------
if(USE_BRAINSTransformConvert)
  add_subdirectory(BRAINSTransformConvert)
endif()

endif(ITK_VERSION_MAJOR GREATER 3)

