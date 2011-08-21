#  This file contains common Macros that are useful for controlling
#  hierarchial build where common components may be included in more
#  than one subdirectory, but we only want a single copy of it built.


###############################################################################
###############################################################################
## A simple macro to set variables ONLY if it has not been set
## This is needed when stand-alone packages are combined into
## a larger package, and the desired behavior is that all the
## binary results end up in the combined directory.
if(NOT SETIFEMPTY)
macro(SETIFEMPTY)
  set(KEY ${ARGV0})
  set(VALUE ${ARGV1})
  if(NOT ${KEY})
    set(${ARGV})
  endif(NOT ${KEY})
endmacro(SETIFEMPTY KEY VALUE)
endif(NOT SETIFEMPTY)

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
endmacro(MakeTestDriverFromSEMTool SEMToolName)

