#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
# #-----------------------------------------------------------------------------
if(POLICY CMP0016)
  cmake_policy(SET CMP0016 NEW)
endif()
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

#
#-----------------------------------------------------------------------------
# Build the optional DEBUGIMAGEVIEWER
if(NOT SETOPTIONALDEBUGIMAGEVIEWER)
macro(SETOPTIONALDEBUGIMAGEVIEWER)
if(BRAINS_BUILD)
  option(USE_DEBUG_IMAGE_VIEWER "Use the DEBUG_IMAGE_VIEWER for debugging" ON)
else(BRAINS_BUILD)
  option(USE_DEBUG_IMAGE_VIEWER "Use the DEBUG_IMAGE_VIEWER for debugging" OFF)
endif(BRAINS_BUILD)

mark_as_advanced(USE_DEBUG_IMAGE_VIEWER)
set(OPTIONAL_DEBUG_LINK_LIBRARIES) ## Set it to empty as the default
if( USE_DEBUG_IMAGE_VIEWER )
   if(NOT KWWidgets_SOURCE_DIR)
     find_package(KWWidgets REQUIRED)
     include(${KWWidgets_USE_FILE})
   endif(NOT KWWidgets_SOURCE_DIR)
   add_definitions(-DUSE_DEBUG_IMAGE_VIEWER)
   find_path(DEBUG_IMAGE_VIEWER_INCLUDE_DIR DebugImageViewerClient.h ${CMAKE_INSTALL_PREFIX}/include)
   include_directories(${DEBUG_IMAGE_VIEWER_INCLUDE_DIR})
   set(OPTIONAL_DEBUG_LINK_LIBRARIES ${KWWidgets_LIBRARIES})
endif( USE_DEBUG_IMAGE_VIEWER )
endmacro(SETOPTIONALDEBUGIMAGEVIEWER)
endif(NOT SETOPTIONALDEBUGIMAGEVIEWER)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## A macro to create CLP programs for Slicer or BRAINS3

if(NOT StandardBRAINSBuildMacro)
macro(StandardBRAINSBuildMacro)
  set(options
  )
  set(oneValueArgs
    NAME
  )
  set(multiValueArgs
    ADDITIONAL_SRCS
    TARGET_LIBRARIES
    # LINK_DIRECTORIES
    # INCLUDE_DIRECTORIES
  )
  CMAKE_PARSE_ARGUMENTS(BRAINS_SEM
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
###
# if(NOT INTEGRATE_WITH_SLICER)
#   set(EXECUTABLE_ONLY_FLAG EXECUTABLE_ONLY)
# endif(NOT INTEGRATE_WITH_SLICER)

SEMMacroBuildCLI(
  NAME "${BRAINS_SEM_NAME}"
#  ${EXECUTABLE_ONLY_FLAG}
    ADDITIONAL_SRCS "${BRAINS_SEM_ADDITIONAL_SRCS}"
    LOGO_HEADER "${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h"
    TARGET_LIBRARIES "${BRAINS_SEM_TARGET_LIBRARIES}"
    CLI_SHARED_LIBRARY_WRAPPER_CXX "${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx"
    VERBOSE
    RUNTIME_OUTPUT_DIRECTORY    "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    LIBRARY_OUTPUT_DIRECTORY    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY    "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    INSTALL_RUNTIME_DESTINATION "${CMAKE_INSTALL_RUNTIME_DESTINATION}"
    INSTALL_LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBRARY_DESTINATION}"
    INSTALL_ARCHIVE_DESTINATION "${CMAKE_INSTALL_ARCHIVE_DESTINATION}"
)
endmacro(StandardBRAINSBuildMacro TOOLNAME)
endif(NOT StandardBRAINSBuildMacro)

###############################################################################
###############################################################################
## MakeTestDriverFromSEMTool
## For tools made with the slicer execution model,
## This macro will build a test driver that adds the
## --compare
## --compareIntensityTolerance
## --compareRadiusTolerance
## --compareNumberOfPixelsTolerance
## to the SEM tools.
macro(MakeTestDriverFromSEMTool SEMToolName SEMToolTestSourceName)
  set(SEMToolLibName        ${SEMToolName}Lib)

  if(ITK_VERSION_MAJOR LESS 4)
    ## BackPort files from ITKv4 need to be pushed to ITKv3 for backwards compatibility
    include_directories(${BRAINSTools_SOURCE_DIR}/BRAINSCommonLib/itkV3TestKernel/include)
  endif()

  set(CMAKE_TESTDRIVER_BEFORE_TESTMAIN "#include \"itkTestDriverBeforeTest.inc\"")
  set(CMAKE_TESTDRIVER_AFTER_TESTMAIN "#include \"itkTestDriverAfterTest.inc\"")

  create_test_sourcelist(${SEMToolName}   ${SEMToolName}TestDriver.cxx ${SEMToolTestSourceName}
     EXTRA_INCLUDE itkTestDriverIncludeRequiredIOFactories.h
     FUNCTION  ProcessArgumentsAndRegisterRequiredFactories
     )
  add_executable(${SEMToolName}TestDriver ${SEMToolName}TestDriver.cxx ${SEMToolTestSourceName})
  target_link_libraries(${SEMToolName}TestDriver ${SEMToolLibName} ${ITKTestKernel_LIBRARIES})
  set_target_properties(${SEMToolName}TestDriver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endmacro(MakeTestDriverFromSEMTool SEMToolName)

