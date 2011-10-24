
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_SOURCE_DIR}/CMake
  ${${PROJECT_NAME}_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )

#-----------------------------------------------------------------------------
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel NO_MODULE REQUIRED GenerateCLP)
include(${SlicerExecutionModel_USE_FILE})

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

#-----------------------------------------------------------------------------
# TODO Should be moved in a subdirectory
#link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

#-----------------------------------------------------------------------------
# Add needed flag for gnu on linux like enviroments to build static common libs
# suitable for linking with shared object libs.
#if(NOT APPLE AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
#  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "-fPIC")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
#  endif()
#  if(NOT "${CMAKE_C_FLAGS}" MATCHES "-fPIC")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
#  endif()
#endif()

##-----------------------------------------------------------------------
## Setup locations to find externally maintained test data.
##-----------------------------------------------------------------------
list(APPEND ExternalData_URL_TEMPLATES
  # Local data store populated by the ITK pre-commit hook
  "file:///${${PROJECT_NAME}_SOURCE_DIR}/.ExternalData/%(algo)/%(hash)"
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

include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/CMakeBuildMacros.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/CMakeBRAINS3BuildMacros.cmake)
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/IJMacros.txt)

set(ExternalData_SOURCE_ROOT ${${PROJECT_NAME}_SOURCE_DIR})
include(${BRAINSCommonLib_BUILDSCRIPTS_DIR}/ExternalData.cmake)

add_subdirectory(BRAINSCommonLib)
set(BRAINSCommonLib_DIR ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/BRAINSCommonLib ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCommonLib)

#find_package(BRAINSCommonLib NO_MODULE REQUIRED)
#include(${BRAINSCommonLib_USE_FILE})

#-----------------------------------------------------------------------------
# Define list of module names
#-----------------------------------------------------------------------------
set(brains_modulenames
  BRAINSFit
  BRAINSResample
  BRAINSROIAuto
  GTRACT
  ImageCalculator
  )

## Tools that only work with ITKv4
if(ITK_VERSION_MAJOR GREATER 3)
  list(APPEND brains_modulenames
    BRAINSConstellationDetector
    BRAINSABC
    BRAINSMush
    BRAINSMultiModeSegment
    BRAINSInitializedControlPoints
    BRAINSTransformConvert
    # Not yet ready BRAINSCut
    )
else(ITK_VERSION_MAJOR GREATER 3)
endif(ITK_VERSION_MAJOR GREATER 3)
## HACK:  Need to have it available just for testing puposes.
  list(APPEND brains_modulenames
    BRAINSDemonWarp  # This is only working in ITKv3,  ITKv4 does not work correctly, and we are moving to a new program.
  )

#-----------------------------------------------------------------------------
# Add module sub-directory if USE_<MODULENAME> is both defined and true
#-----------------------------------------------------------------------------
foreach(modulename ${brains_modulenames})
  if(DEFINED USE_${modulename} AND USE_${modulename})
    add_subdirectory(${modulename})
  endif()
endforeach()

