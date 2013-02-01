#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
# #-----------------------------------------------------------------------------
if(POLICY CMP0016)
  cmake_policy(SET CMP0016 NEW)
endif()
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## A macro to create CLP programs for Slicer or BRAINS3

if(NOT StandardBRAINSBuildMacro)
  include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)
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
  if(NOT INTEGRATE_WITH_SLICER AND NOT BUILD_SHARED_LIBS)
    set(BRAINSTools_BUILD_STATIC_LIBRARY BUILD_STATIC_LIBRARY)
  endif()

  if( NOT INTEGRATE_WITH_SLICER AND NOT BUILD_TESTING )
    ## The testing framework requires the shared lib infrastructure.
    ## to avoid building the executables twice.
    set(BRAINSTools_EXECUTABLE_ONLY EXECUTABLE_ONLY)
  endif()

    SEMMacroBuildCLI(
      NAME "${BRAINS_SEM_NAME}"
        ${BRAINSTools_EXECUTABLE_ONLY}
        ${BRAINSTools_BUILD_STATIC_LIBRARY}
        ADDITIONAL_SRCS "${BRAINS_SEM_ADDITIONAL_SRCS}"
        LOGO_HEADER "${BRAINSCommonLib_BUILDSCRIPTS_DIR}/BRAINSLogo.h"
        TARGET_LIBRARIES "${BRAINS_SEM_TARGET_LIBRARIES}"
        CLI_LIBRARY_WRAPPER_CXX "${BRAINSCommonLib_BUILDSCRIPTS_DIR}/SEMCommanLineSharedLibraryWrapper.cxx"
        RUNTIME_OUTPUT_DIRECTORY ${BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${BRAINSTools_CLI_LIBRARY_OUTPUT_DIRECTORY}
        ARCHIVE_OUTPUT_DIRECTORY ${BRAINSTools_CLI_ARCHIVE_OUTPUT_DIRECTORY}
        INSTALL_RUNTIME_DESTINATION ${BRAINSTools_CLI_INSTALL_RUNTIME_DESTINATION}
        INSTALL_LIBRARY_DESTINATION ${BRAINSTools_CLI_INSTALL_LIBRARY_DESTINATION}
        INSTALL_ARCHIVE_DESTINATION ${BRAINSTools_CLI_INSTALL_ARCHIVE_DESTINATION}
    #    VERBOSE
    )
  endmacro()
endif()

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
endmacro()

# DebugImageViewer Macro
if(USE_DebugImageViewer)

  if(NOT VTK_FOUND) #See if it has already been found because this is a nested project
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
  endif()

  macro(DebugImageViewerLibAdditions LibVariable)
    list(APPEND ${LibVariable} ${VTK_LIBRARIES})
  endmacro()
else()
  macro(DebugImageViewerLibAdditions LibVariable)
  endmacro()

endif()
