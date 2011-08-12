cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

enable_testing()
include(CTest)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel NO_MODULE REQUIRED GenerateCLP)
if(GenerateCLP_DIR)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_DIR)
  message(FATAL_ERROR "Can't build without GenerateCLP. Please set GenerateCLP_DIR")
endif(GenerateCLP_DIR)
include(${SlicerExecutionModel_MACROS}/SEMMacroBuildCLI.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts/PreventInSourceBuilds.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts/CMakeBuildMacros.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts/CMakeBRAINS3BuildMacros.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts/IJMacros.txt)

###
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

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

include(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib/BuildScripts/ExternalData.cmake)
set(TestData_DIR ${CMAKE_CURRENT_SOURCE_DIR}/TestData)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## NOTE:  BRAINSCommonLib is REQUIRED.
# BRAINSCommonLib
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
add_subdirectory(BRAINSCommonLib)
set(BRAINSCommonLib_DIR ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib)

find_package(BRAINSCommonLib NO_MODULE REQUIRED)
include(${BRAINSCommonLib_USE_FILE})

#-----------------------------------------------------------------------------
# BRAINSFit
#-----------------------------------------------------------------------------
if(USE_BRAINSFit)
  add_subdirectory(BRAINSFit)
endif()

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
# GTRACT
#-----------------------------------------------------------------------------
if(USE_GTRACT)
  add_subdirectory(GTRACT)
endif()

#  BuildExtPackage(BRAINSCut "BRAINSCommonLib;${OpenCV_DEPEND}" )

