
##- project(BRAINSSuperResolution)

#-----------------------------------------------------------------------------
# Dependencies.
#
#-----------------------------------------------------------------------------
# Output directories.
#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/BRAINSSuperResolutionCommonLib.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/BRAINSSuperResolutionCommonLib.h
  )

set(ALL_PROGS_LIST
  GenerateEdgeMapImage
  #BRAINSSuperResolution
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

### Add the Testing Subdirectory.
#if(BUILD_TESTING AND NOT BRAINSTools_DISABLE_TESTING)
#  add_subdirectory(TestSuite)
#endif()
